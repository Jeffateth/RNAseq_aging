"""Builds download_raw_data.ipynb (run: python scripts/make_download_notebook.py)."""
import nbformat as nbf

C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: C.append(nbf.v4.new_code_cell(s.strip("\n")))

md("""
# Resumable raw-data download (27 SRA runs) - run locally, hand heavy work to Euler later

**What this notebook does.** It downloads the 27 selected mouse-cerebellum runs (PRJNA875066 + PRJNA946653) with the SRA Toolkit,
optionally converts them to FASTQ, and tells you at every step how much disk and time is involved.

**Why it is safe on flaky Wi-Fi / a closing laptop lid**
| Risk | Protection |
|---|---|
| Wi-Fi drops mid-download | The worker waits until NCBI is reachable again, then retries with backoff (up to 20 attempts per run). |
| Connection silently hangs | Watchdog kills `prefetch` if the partial file stops growing for 5 min, then retries. |
| Corrupt / truncated file | A run counts as done only if `vdb-validate` passes; otherwise it is deleted and retried. |
| Notebook/kernel closed | The worker runs **detached in the background**; reopen the notebook and rerun the *Monitor* cell. |
| Mac goes to sleep | The worker is wrapped in `caffeinate -i` (keeps the Mac awake while it works; keep the charger plugged in). |
| Restart from scratch | Never needed: rerunning any cell skips finished runs (state lives on disk in `raw_data/`). |
| Accidental 50 GB download | Nothing starts without an explicit `APPROVE_GB` that covers the remaining size. |

**Recommended order:** 1 Setup -> 2 Preflight -> 3 Plan -> 4 Test run (1 run) -> 5 Monitor -> *(your decision)* 6 All runs -> 7 optional FASTQ -> 8 Verify -> 9 Euler hand-off.
Raw `.sra` files are never modified.
""")

md("## 1. Setup\nEdit only the settings in this cell. Everything is off by default.")
code("""
import os, sys, json, shutil, subprocess, time
from pathlib import Path
import pandas as pd

# --- find the project folder (the one that contains metadata/srr_list.txt) ---
ROOT = next((p for p in [Path.cwd(), Path.cwd().parent, Path.home()/"RNAseq_aging", Path.home()/"Desktop/RNAseq_aging",
                         Path.home()/"Desktop/RNAseq_tutorial/RNAseq_aging"] if (p/"metadata/srr_list.txt").exists()), None)
assert ROOT, "Could not find the project folder. Set ROOT = Path('/path/to/RNAseq_aging') by hand."
WORKER = ROOT/"scripts/sra_resumable.py"
RAW, STATE = ROOT/"raw_data", ROOT/"raw_data/state.json"
os.environ["PATH"] = "/opt/homebrew/bin:/usr/local/bin:" + os.environ["PATH"]

# --- settings ---------------------------------------------------------------
TEST_RUN   = "SRR21355398"   # smallest run (~0.7 GB) for the workflow test
APPROVE_GB = 1.0             # maximum GB you allow the worker to download in the next launch (raise it deliberately)
CONVERT    = False           # also convert to FASTQ.gz (needs ~7x the .sra size free per run)
THREADS    = 6               # conversion threads (M2 Pro: 6 is a good, cool setting)
print("project:", ROOT)
""")

md("""
## 2. Preflight: tools, disk, network
`prefetch` / `vdb-validate` / `fasterq-dump` come from SRA Toolkit (`brew install sratoolkit`). The disk must hold the `.sra` files,
and, if you convert, an additional ~65 GB of FASTQ plus temporary space.
""")
code("""
for t in ["prefetch", "vdb-validate", "fasterq-dump", "caffeinate"]:
    print(f"{t:14}", shutil.which(t) or "MISSING")
free = shutil.disk_usage(RAW if RAW.exists() else ROOT).free / 1024**3
print(f"free disk: {free:.0f} GB")
import urllib.request
try: urllib.request.urlopen("https://ftp.ncbi.nlm.nih.gov/", timeout=8).close(); print("NCBI reachable: yes")
except Exception as e: print("NCBI reachable: NO ->", e)
""")

md("""
## 3. Plan: what is downloaded and how big
Sizes come from SRA metadata (`metadata/sra_metadata.csv`). *Free after all downloads* tells you whether the whole set fits.
""")
code("""
meta = pd.read_csv(ROOT/"metadata/sra_metadata.csv")[["SRR","BioProject","sex","age_months","condition","sra_size_MB"]]
meta["size_GB"] = (meta.sra_size_MB/1024).round(2)
total = meta.size_GB.sum()
print(f"total .sra: {total:.1f} GB | +FASTQ.gz if converted: ~{total*1.3:.0f} GB | free disk now: {free:.0f} GB")
print("-> everything fits (.sra only)" if free > total + 10 else "-> NOT enough space for everything")
print("-> .sra + FASTQ would NOT fit together; convert run-by-run and delete .sra afterwards, or use Euler" if free < total*2.3+20 else "")
meta.sort_values("size_GB")
""")

md("""
## 4. Helpers: launch the worker in the background, and monitor it
`launch()` starts the worker as a detached process (it keeps running if you close this notebook).
`monitor()` shows live status; **interrupting the monitor (stop button) does not stop the download.**
""")
code("""
def worker_pid():
    lock = RAW/"worker.lock"
    try: pid = int(lock.read_text()); os.kill(pid, 0); return pid
    except Exception: return None

def launch(runs=(), all_runs=False, approve_gb=APPROVE_GB, convert=CONVERT):
    if worker_pid(): print("a worker is already running, pid", worker_pid()); return
    cmd = [sys.executable, "-u", str(WORKER), *runs, "--approve-gb", str(approve_gb), "--threads", str(THREADS)]
    if all_runs: cmd.append("--all")
    if convert:  cmd.append("--convert")
    if shutil.which("caffeinate"): cmd = ["caffeinate", "-i"] + cmd     # prevent idle sleep while working
    (RAW/"logs").mkdir(parents=True, exist_ok=True)
    out = open(RAW/"logs/worker_stdout.log", "a")
    subprocess.Popen(cmd, stdout=out, stderr=subprocess.STDOUT, start_new_session=True)
    time.sleep(3); print("launched:", " ".join(cmd[-6:]), "| pid", worker_pid() or "(exited - check log below)")
    print(subprocess.run(["tail", "-n", "5", str(RAW/"logs/worker_stdout.log")], capture_output=True, text=True).stdout)

def status():
    print(subprocess.run([sys.executable, str(WORKER), "--status"], capture_output=True, text=True).stdout)

def monitor(refresh=20, tail=6):
    from IPython.display import clear_output
    try:
        while True:
            clear_output(wait=True)
            print(time.strftime("%H:%M:%S"), "| worker pid:", worker_pid() or "not running"); status()
            print(subprocess.run(["tail", "-n", str(tail), str(RAW/"logs/sra_resumable.log")], capture_output=True, text=True).stdout)
            if not worker_pid(): print("worker finished/stopped."); break
            time.sleep(refresh)
    except KeyboardInterrupt:
        print("monitor stopped; the worker keeps running in the background.")

def stop():
    pid = worker_pid()
    if pid: os.kill(pid, 15); print("stop signal sent; progress is saved, rerun launch() to resume")
""")

md("""
## 5. Test the workflow with ONE run (~0.7 GB)
This stays under the 1 GB threshold. Start it, watch it, and (optionally) test resilience: turn Wi-Fi off for a minute and back on;
the log should say *waiting for connection* and then *network is back*.
""")
code("launch(runs=[TEST_RUN], approve_gb=1.0)")
code("monitor()   # safe to interrupt; rerun any time, also after restarting the notebook")

md("""
## 6. Download all 27 runs (~50 GB) - only after your explicit decision
Check the plan above first. To authorise it, set `APPROVE_GB` to at least the total (55 is enough), run the cell once, and then use *monitor()*.
It continues where it stopped if interrupted: simply run this cell again.
""")
code("""
RUN_ALL = False   # <-- change to True deliberately
if RUN_ALL:
    launch(all_runs=True, approve_gb=55, convert=False)
else:
    print("RUN_ALL is False - nothing started.")
""")

md("""
## 7. (Optional) Convert to FASTQ
Only needed if you want to map reads yourself. Do **not** convert all runs locally unless you have the space:
each run needs ~7x its `.sra` size temporarily, and the final `.fastq.gz` pair is ~1.3x. Convert a few runs, or let Euler do it (section 9).
Output goes to `raw_data/fastq/`, written atomically (a file only appears when complete).
""")
code("""
CONVERT_RUNS = []   # e.g. ["SRR21355398"]
if CONVERT_RUNS:
    launch(runs=CONVERT_RUNS, approve_gb=APPROVE_GB, convert=True)
else:
    print("no runs selected for conversion")
""")

md("""
## 8. Verify and write a manifest
Re-validates every downloaded `.sra` and records sizes in `metadata/raw_data_manifest.csv`, so the state of your raw data is documented.
""")
code("""
rows = []
for srr in (ROOT/"metadata/srr_list.txt").read_text().split():
    p = RAW/"sra"/srr/f"{srr}.sra"
    ok = False
    if p.exists():
        r = subprocess.run(["vdb-validate", str(p)], capture_output=True, text=True)
        ok = r.returncode == 0 and "inconsistent" not in (r.stdout+r.stderr).lower()
    rows.append(dict(SRR=srr, sra_present=p.exists(), sra_valid=ok, sra_GB=round(p.stat().st_size/1024**3, 2) if p.exists() else None,
                     fastq_pair=(RAW/"fastq"/f"{srr}_1.fastq.gz").exists() and (RAW/"fastq"/f"{srr}_2.fastq.gz").exists()))
man = pd.DataFrame(rows); man.to_csv(ROOT/"metadata/raw_data_manifest.csv", index=False)
print(f"{man.sra_valid.sum()}/27 valid .sra, {man.fastq_pair.sum()}/27 FASTQ pairs"); man
""")

md("""
## 9. Hand-off to Euler (heavy steps)
Local machine: download, QC, exploration. Euler: conversion, mapping, anything RAM-heavy (STAR mouse index needs ~30 GB RAM).
The cell below writes a SLURM array script (one task per run) that downloads and converts on Euler's scratch (`$SCRATCH`; files older than 15 days are deleted there,
so map and move results soon). Nothing is submitted from here.

On Euler (after `ssh <user>@euler.ethz.ch`, with SRA Toolkit installed in conda as in the tutorial, `sra-tools=3.1`):
```bash
mkdir -p $SCRATCH/rnaseq_aging/{logs,rawdata} && cd $SCRATCH/rnaseq_aging
# copy metadata/srr_list.txt and this script here (scp), then:
sbatch euler_fetch_array.sbatch
```
Next steps on Euler (not generated here): STAR index for mouse (GRCm39 + GENCODE vM-release) and `STAR --quantMode GeneCounts`, or kallisto for a light-weight alternative.
Note the GEO count tables already exist, so remapping is only worthwhile for QC, isoforms or a consistent pipeline.
""")
code("""
sbatch = '''#!/bin/bash
#SBATCH --job-name=sra_fetch
#SBATCH --array=1-27%6
#SBATCH --ntasks=1 --cpus-per-task=6
#SBATCH --mem-per-cpu=3G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out

set -euo pipefail
cd $SCRATCH/rnaseq_aging
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" srr_list.txt)
# resumable: skip if FASTQ pair exists; prefetch resumes partial files; retried on failure
if [[ -s rawdata/${SRR}_1.fastq.gz && -s rawdata/${SRR}_2.fastq.gz ]]; then echo "$SRR done"; exit 0; fi
for i in 1 2 3 4 5; do
  prefetch "$SRR" --output-directory sra --max-size 100G && vdb-validate "sra/$SRR/$SRR.sra" && break
  echo "attempt $i failed, retrying"; sleep $((30*i))
done
fasterq-dump "sra/$SRR/$SRR.sra" --split-3 --threads 6 --temp $TMPDIR --outdir rawdata
gzip -f rawdata/${SRR}_*.fastq
'''
(ROOT/"scripts/euler_fetch_array.sbatch").write_text(sbatch)
print("written:", ROOT/"scripts/euler_fetch_array.sbatch")
""")

md("""
## Troubleshooting
- **`prefetch not found`** -> `brew install sratoolkit`, restart the notebook.
- **Worker says "another worker is running"** -> it really is (use `monitor()`), or call `stop()`. Stale locks are detected automatically.
- **Repeated "no progress ... killing stalled process"** -> unstable network or NCBI throttling; the worker keeps retrying. Try again later or move to Euler.
- **Log files:** `raw_data/logs/sra_resumable.log` (what happened), `raw_data/logs/tool_output.log` (prefetch/fasterq-dump output).
- **Start a run over:** delete `raw_data/sra/<SRR>/` (and its entry in `raw_data/state.json` if you like).
""")

nb = nbf.v4.new_notebook(cells=C, metadata={"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"}})
nbf.validate(nb)
nbf.write(nb, "download_raw_data.ipynb")
print("wrote download_raw_data.ipynb")
