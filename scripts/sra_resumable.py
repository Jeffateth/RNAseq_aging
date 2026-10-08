#!/usr/bin/env python3
"""Resumable, Wi-Fi-tolerant SRA download (prefetch) and optional FASTQ conversion (fasterq-dump).

Standard library only. Designed to run detached in the background (see download_raw_data.ipynb),
but also usable directly:

    python scripts/sra_resumable.py SRR21355398 --approve-gb 1        # one test run
    python scripts/sra_resumable.py --all --approve-gb 55             # all 27 runs
    python scripts/sra_resumable.py --all --approve-gb 55 --convert   # + FASTQ conversion
    python scripts/sra_resumable.py --status

Resumability, in layers:
  1. Truth is the filesystem: a run counts as done only if its .sra passes `vdb-validate`.
     Re-running the script skips finished runs and continues the rest.
  2. Per-run retry loop with backoff; before each attempt it waits until NCBI is reachable.
  3. Stall watchdog: if the partial download stops growing for --stall-sec, prefetch is killed and retried
     (prefetch resumes partial files in recent SRA Toolkit versions; check the log to confirm).
  4. state.json is written atomically, and a lock file stops two workers from running at once.
Raw .sra files are never modified or deleted by this script.
"""
import argparse, csv, json, os, shutil, signal, subprocess, sys, time, urllib.request
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "raw_data"
SRA_DIR, FASTQ_DIR, TMP_DIR, LOG_DIR = RAW / "sra", RAW / "fastq", RAW / "fastq_tmp", RAW / "logs"
STATE, LOCK = RAW / "state.json", RAW / "worker.lock"
SRR_LIST, META = ROOT / "metadata/srr_list.txt", ROOT / "metadata/sra_metadata.csv"
PROBE_URL = "https://ftp.ncbi.nlm.nih.gov/"
GB = 1024 ** 3
CONVERT_SPACE_FACTOR = 7          # rough uncompressed FASTQ + temp space needed relative to .sra size

os.environ["PATH"] = "/opt/homebrew/bin:/usr/local/bin:" + os.environ.get("PATH", "")
_child = None


def log(msg):
    line = f"[{datetime.now():%Y-%m-%d %H:%M:%S}] {msg}"
    print(line, flush=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with open(LOG_DIR / "sra_resumable.log", "a") as f:
        f.write(line + "\n")


def online(timeout=8):
    try:
        urllib.request.urlopen(PROBE_URL, timeout=timeout).close()
        return True
    except Exception:
        return False


def wait_online():
    if online():
        return
    log("network unreachable - waiting for connection (checks every 15 s)")
    while not online():
        time.sleep(15)
    log("network is back")


def load_state():
    return json.loads(STATE.read_text()) if STATE.exists() else {}


def save_state(st):
    RAW.mkdir(parents=True, exist_ok=True)
    tmp = STATE.with_suffix(".tmp")
    tmp.write_text(json.dumps(st, indent=1))
    tmp.replace(STATE)


def set_state(st, srr, **kw):
    st.setdefault(srr, {}).update(kw, updated=datetime.now().isoformat(timespec="seconds"))
    save_state(st)


def sizes_mb():
    with open(META) as f:
        return {r["SRR"]: float(r["sra_size_MB"]) for r in csv.DictReader(f)}


def dir_bytes(p):
    return sum(f.stat().st_size for f in Path(p).rglob("*") if f.is_file()) if Path(p).exists() else 0


def sra_path(srr):
    return SRA_DIR / srr / f"{srr}.sra"


def validate(srr):
    p = sra_path(srr)
    if not p.exists() or p.stat().st_size == 0:
        return False
    r = subprocess.run(["vdb-validate", str(p)], capture_output=True, text=True)
    out = (r.stdout + r.stderr).lower()
    return r.returncode == 0 and "inconsistent" not in out and "corrupt" not in out


def fastq_done(srr):
    return (FASTQ_DIR / f"{srr}_1.fastq.gz").exists() and (FASTQ_DIR / f"{srr}_2.fastq.gz").exists()


def run_watched(cmd, watch_dir, stall_sec):
    """Run cmd; kill it if watch_dir stops growing for stall_sec. Returns exit code (or -1 if stalled)."""
    global _child
    _child = subprocess.Popen(cmd, stdout=open(LOG_DIR / "tool_output.log", "a"), stderr=subprocess.STDOUT)
    last, last_change = dir_bytes(watch_dir), time.time()
    while _child.poll() is None:
        time.sleep(10)
        cur = dir_bytes(watch_dir)
        if cur != last:
            last, last_change = cur, time.time()
        elif time.time() - last_change > stall_sec:
            log(f"no progress for {stall_sec}s - killing stalled process")
            _child.kill()
            _child.wait()
            _child = None
            return -1
    rc = _child.returncode
    _child = None
    return rc


def backoff(attempt):
    t = min(30 * attempt, 600)
    log(f"retrying in {t}s")
    time.sleep(t)


def download(srr, st, max_attempts, stall_sec):
    if validate(srr):
        set_state(st, srr, download="done")
        log(f"{srr}: already downloaded and valid, skipping")
        return True
    for attempt in range(1, max_attempts + 1):
        wait_online()
        set_state(st, srr, download="running", attempts=attempt)
        log(f"{srr}: prefetch attempt {attempt}/{max_attempts}")
        rc = run_watched(["prefetch", srr, "--output-directory", str(SRA_DIR), "--max-size", "100G"],
                         SRA_DIR / srr, stall_sec)
        if rc == 0 and validate(srr):
            set_state(st, srr, download="done", error="")
            log(f"{srr}: OK ({sra_path(srr).stat().st_size / GB:.2f} GB, validated)")
            return True
        if rc == 0:
            log(f"{srr}: file failed validation - removing it and retrying")
            shutil.rmtree(SRA_DIR / srr, ignore_errors=True)
        set_state(st, srr, download="retrying", error=f"exit code {rc}")
        backoff(attempt)
    set_state(st, srr, download="failed")
    log(f"{srr}: FAILED after {max_attempts} attempts (rerun to continue)")
    return False


def convert(srr, st, size_mb, threads, max_attempts):
    if fastq_done(srr):
        set_state(st, srr, convert="done")
        log(f"{srr}: FASTQ already present, skipping")
        return True
    need = size_mb * 1024 ** 2 * CONVERT_SPACE_FACTOR + 5 * GB
    free = shutil.disk_usage(RAW).free
    if free < need:
        log(f"{srr}: not enough free disk ({free / GB:.0f} GB < ~{need / GB:.0f} GB needed) - skipping")
        set_state(st, srr, convert="skipped_no_space")
        return False
    for attempt in range(1, max_attempts + 1):
        out = TMP_DIR / srr
        shutil.rmtree(out, ignore_errors=True)
        (out / "tmp").mkdir(parents=True)
        set_state(st, srr, convert="running")
        log(f"{srr}: fasterq-dump attempt {attempt}/{max_attempts}")
        rc = run_watched(["fasterq-dump", str(sra_path(srr)), "--split-3", "--threads", str(threads),
                          "--temp", str(out / "tmp"), "--outdir", str(out)], out, 1800)
        files = sorted(out.glob("*.fastq"))
        if rc == 0 and len(files) == 2:
            gz = shutil.which("pigz") or "gzip"
            subprocess.run([gz, "-f", *map(str, files)], check=True)
            FASTQ_DIR.mkdir(exist_ok=True)
            for f in sorted(out.glob("*.fastq.gz")):
                f.replace(FASTQ_DIR / f.name)       # move only when complete -> no half-written finals
            shutil.rmtree(out, ignore_errors=True)
            set_state(st, srr, convert="done")
            log(f"{srr}: FASTQ written to {FASTQ_DIR}")
            return True
        log(f"{srr}: conversion failed (exit {rc}, {len(files)} fastq files)")
        shutil.rmtree(out, ignore_errors=True)
        backoff(attempt)
    set_state(st, srr, convert="failed")
    return False


def print_status():
    st, sz = load_state(), sizes_mb()
    print(f"{'run':13}{'size GB':>8}  {'download':10}{'convert':10}{'on disk GB':>11}  note")
    for srr in SRR_LIST.read_text().split():
        s = st.get(srr, {})
        # cheap: trust the state file written after validation (full re-validation is in notebook section 8)
        d = "done" if (s.get("download") == "done" and sra_path(srr).exists()) else s.get("download", "pending")
        part = dir_bytes(SRA_DIR / srr) / GB
        print(f"{srr:13}{sz[srr] / 1024:8.2f}  {d:10}{('done' if fastq_done(srr) else s.get('convert', '-')):10}{part:11.2f}  {s.get('error', '')}")


def acquire_lock():
    if LOCK.exists():
        try:
            pid = int(LOCK.read_text())
            os.kill(pid, 0)
            sys.exit(f"another worker is running (pid {pid}). Wait for it or kill it.")
        except (ValueError, ProcessLookupError, PermissionError):
            pass                                    # stale lock
    RAW.mkdir(parents=True, exist_ok=True)
    LOCK.write_text(str(os.getpid()))


def release(*_):
    if _child and _child.poll() is None:
        _child.kill()
    LOCK.unlink(missing_ok=True)
    log("worker stopped (progress is saved; rerun to resume)")
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--convert", action="store_true", help="also convert to FASTQ.gz")
    ap.add_argument("--approve-gb", type=float, default=1.0, help="max total GB you approve to download now")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--max-attempts", type=int, default=20)
    ap.add_argument("--stall-sec", type=int, default=300)
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()

    if a.status:
        return print_status()
    allowed = SRR_LIST.read_text().split()
    runs = allowed if a.all else a.runs
    if not runs:
        sys.exit("give run accessions or --all")
    bad = [r for r in runs if r not in allowed]
    if bad:
        sys.exit(f"not in the selected 27-run list: {bad}")
    for tool in ("prefetch", "vdb-validate") + (("fasterq-dump",) if a.convert else ()):
        if not shutil.which(tool):
            sys.exit(f"{tool} not found (brew install sratoolkit)")

    sz = sizes_mb()
    todo_gb = sum(sz[r] for r in runs if not sra_path(r).exists()) / 1024
    if todo_gb > a.approve_gb:
        sys.exit(f"remaining download is {todo_gb:.1f} GB but only --approve-gb {a.approve_gb} was given. Refusing.")
    log(f"plan: {len(runs)} runs, ~{todo_gb:.1f} GB still to download, convert={a.convert}")

    acquire_lock()
    signal.signal(signal.SIGTERM, release)
    signal.signal(signal.SIGINT, release)
    st = load_state()
    fails = 0
    for srr in runs:
        if download(srr, st, a.max_attempts, a.stall_sec):
            if a.convert and not convert(srr, st, sz[srr], a.threads, 3):
                fails += 1
        else:
            fails += 1
    LOCK.unlink(missing_ok=True)
    log(f"finished: {len(runs) - fails} ok, {fails} not finished")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
