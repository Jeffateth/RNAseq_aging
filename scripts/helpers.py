"""Python helpers for notebooks 00, 01 and 02. Every pipeline step is idempotent: finished work is skipped.

Contents
  Setup and paths ..... ROOT, SRA, FASTQ, QC, KALLISTO, STAR ...; load_runs() reads metadata/samples.csv
  Download ............ download_all (prefetch with progress bar), validate_all, status
  QC and trimming ..... fastq (.sra -> FASTQ), fastqc, multiqc, fastqc_summary, trim (fastp)
  Quantification ...... build_index, quant (kallisto), process (one run end to end), check_consistent, qc_summary
  Euler results ....... unpack_euler (tarballs copied back from the cluster), load_star (STAR counts + mapping statistics)
Terms: docs/glossary.md."""
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

# kallisto/fastp live in the pixi env (not activated when VS Code starts the kernel); sra-tools come from Homebrew
os.environ["PATH"] = f"{Path(sys.executable).parent}:{os.environ['PATH']}:/opt/homebrew/bin"
ROOT = Path(__file__).resolve().parents[1]
SRA, FASTQ, TMP = ROOT / "raw_data/sra", ROOT / "raw_data/fastq", ROOT / "raw_data/tmp"
QC, KALLISTO, LOGS = ROOT / "results/qc", ROOT / "results/kallisto", ROOT / "results/logs"
REF = ROOT / "reference"
TRANSCRIPTS = REF / "gencode.vM25.transcripts.fa.gz"   # GENCODE vM25 = GRCm38/mm10, as in the paper
INDEX = REF / "gencode.vM25.kallisto.idx"
for d in (SRA, FASTQ, TMP, QC, KALLISTO, LOGS):
    d.mkdir(parents=True, exist_ok=True)


def load_runs() -> pd.DataFrame:
    """The 27 assigned runs with age/sex/group, indexed by SRR."""
    return pd.read_csv(ROOT / "metadata/samples.csv").set_index("SRR")


def sh(cmd: str, log: str) -> None:
    """Run a shell command; stdout/stderr go to results/logs/<log>.log (shown only on failure)."""
    with open(LOGS / f"{log}.log", "a") as f:
        f.write(f"\n$ {cmd}\n")
        f.flush()
        r = subprocess.run(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode:
        tail = "".join(open(LOGS / f"{log}.log").readlines()[-15:])
        raise RuntimeError(f"failed: {cmd}\n{tail}")


def timed(label):
    def deco(fn):
        def wrapper(srr, *a, **k):
            t = time.time()
            out = fn(srr, *a, **k)
            if out != "skip":
                print(f"  {srr} {label}: {time.time() - t:.0f}s")
            return out
        return wrapper
    return deco


# ---- stage checks -----------------------------------------------------------------------
def has_sra(srr):  return (SRA / srr / f"{srr}.sra").is_file()
def has_quant(srr): return (KALLISTO / srr / "abundance.h5").is_file()


def was_trimmed(srr):
    """True/False = quantified from fastp-trimmed/raw reads (read back from kallisto's own record); None if not quantified."""
    import json
    f = KALLISTO / srr / "run_info.json"
    return "_trim_" in json.load(open(f)).get("call", "") if f.is_file() else None


def status(runs: pd.DataFrame) -> pd.DataFrame:
    """Per-run progress table: downloaded -> quantified."""
    return pd.DataFrame({
        "group": runs["group"],
        "sra_GB": [round((SRA / s / f"{s}.sra").stat().st_size / 1e9, 2) if has_sra(s) else None for s in runs.index],
        "quantified": [has_quant(s) for s in runs.index],
        "trimmed": [was_trimmed(s) for s in runs.index],
    })


def check_consistent(runs: pd.DataFrame, trim_reads: bool):
    """Raise if any quantified run was processed differently from the intended setting (avoids a technical confound)."""
    bad = [s for s in runs.index if was_trimmed(s) not in (None, trim_reads)]
    if bad:
        raise RuntimeError(f"{len(bad)} runs were quantified with trimmed={not trim_reads}, expected {trim_reads}: {bad}. "
                           f"Delete results/kallisto/<SRR> for these and re-run.")
    print(f"consistent: all quantified runs used trimmed={trim_reads}")


# ---- stages -----------------------------------------------------------------------------
def _dir_bytes(d):
    return sum(f.stat().st_size for f in d.glob("*") if f.is_file()) if d.is_dir() else 0


def download_all(runs: pd.DataFrame):
    """prefetch every missing .sra (resumes partial files, validates) with a live progress bar in GB.
    Logs: results/logs/<SRR>.log. Safe to interrupt and re-run."""
    from tqdm.auto import tqdm
    exp = runs["sra_size_MB"] / 1e3 * 1.048576      # NCBI lists size_MB in MiB; x1.048576 converts to decimal GB (files on disk match)
    todo = [s for s in runs.index if not has_sra(s)]
    if not todo:
        return print(f"all {len(runs)} runs present")
    bar = tqdm(total=round(exp[todo].sum(), 1), unit="GB", bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} GB [{elapsed}<{remaining}]")
    finished = 0.0
    for i, srr in enumerate(todo, 1):
        bar.set_description(f"{srr} ({i}/{len(todo)})")
        with open(LOGS / f"{srr}.log", "a") as f:
            p = subprocess.Popen(f"prefetch {srr} --output-directory {SRA} --max-size 100G",
                                 shell=True, stdout=f, stderr=subprocess.STDOUT)
            try:
                while p.poll() is None:
                    bar.n = round(min(finished + _dir_bytes(SRA / srr) / 1e9, bar.total), 1)   # never exceed the total
                    bar.refresh()
                    time.sleep(1)
            except BaseException:
                p.kill()                                    # no orphaned prefetch if the cell is interrupted
                raise
        if p.returncode or not has_sra(srr):
            raise RuntimeError(f"prefetch failed for {srr}; see {LOGS / (srr + '.log')}")
        sh(f"vdb-validate {SRA / srr / (srr + '.sra')}", srr)
        finished += exp[srr]
        bar.n = round(min(finished, bar.total), 1)
        bar.refresh()
    bar.close()
    print(f"all {len(runs)} runs present" if all(has_sra(s) for s in runs.index) else "some runs missing")


def validate_all(runs: pd.DataFrame):
    """vdb-validate every .sra; returns the runs that fail (empty list = all OK)."""
    from tqdm.auto import tqdm
    bad = []
    for s in tqdm(runs.index, desc="validate"):
        r = subprocess.run(f"vdb-validate {SRA / s / (s + '.sra')}", shell=True, capture_output=True, text=True)
        if r.returncode or "consistent" not in (r.stdout + r.stderr):
            bad.append(s)
    return bad


@timed("fastq")
def fastq(srr):
    """.sra -> raw paired FASTQ (uncompressed, temporary)."""
    if (FASTQ / f"{srr}_1.fastq").is_file() and (FASTQ / f"{srr}_2.fastq").is_file():
        return "skip"
    sh(f"fasterq-dump {SRA / srr / (srr + '.sra')} --split-3 --threads 6 --temp {TMP} --outdir {FASTQ}", srr)


@timed("fastqc")
def fastqc(srr, n_pairs=None):
    """FastQC on ALL reads of a run (n_pairs=None) or only the first n_pairs pairs -> results/qc/raw/.
    Full runs are converted to a temporary FASTQ (~7 GB) that is deleted afterwards."""
    out = QC / "raw"
    out.mkdir(exist_ok=True)
    if (out / f"{srr}_2_fastqc.zip").is_file():
        return "skip"
    tmp = TMP / f"qc_{srr}"
    tmp.mkdir(exist_ok=True)
    sra = SRA / srr / f"{srr}.sra"
    if n_pairs:
        sh(f"fastq-dump -X {n_pairs} --split-files -O {tmp} {sra}", srr)
    else:
        sh(f"fasterq-dump {sra} --split-3 --threads 6 --temp {tmp} --outdir {tmp}", srr)
    sh(f"fastqc -q -t 4 --outdir {out} {tmp}/{srr}_1.fastq {tmp}/{srr}_2.fastq", srr)
    shutil.rmtree(tmp)


def multiqc():
    """One combined report for all FastQC results -> results/qc/multiqc_raw.html"""
    sh(f"multiqc -q -f -o {QC} -n multiqc_raw {QC / 'raw'}", "multiqc")
    return QC / "multiqc_raw.html"


KEY_MODULES = ["Per base sequence quality", "Per sequence quality scores", "Adapter Content",
               "Overrepresented sequences", "Sequence Duplication Levels"]


def fastqc_summary(srrs, modules=KEY_MODULES) -> pd.DataFrame:
    """FastQC PASS/WARN/FAIL for the key modules (rows: run_mate)."""
    import zipfile
    rows = {}
    for s in srrs:
        for m in (1, 2):
            z = QC / "raw" / f"{s}_{m}_fastqc.zip"
            if z.is_file():
                with zipfile.ZipFile(z) as zf:
                    txt = zf.read(f"{s}_{m}_fastqc/summary.txt").decode()
                rows[f"{s}_{m}"] = {l.split("\t")[1]: l.split("\t")[0] for l in txt.splitlines()}
    return pd.DataFrame(rows).T[modules] if rows else pd.DataFrame(columns=modules)


@timed("trim")
def trim(srr):
    """Optional: fastp adapter/quality trimming of the raw FASTQ -> <SRR>_trim_{1,2}.fastq.gz + fastp report."""
    out1, out2 = FASTQ / f"{srr}_trim_1.fastq.gz", FASTQ / f"{srr}_trim_2.fastq.gz"
    if out1.is_file() and out2.is_file():
        return "skip"
    fastq(srr)
    sh(f"fastp -i {FASTQ}/{srr}_1.fastq -I {FASTQ}/{srr}_2.fastq -o {out1} -O {out2} "
       f"--detect_adapter_for_pe --trim_poly_g --length_required 20 --thread 6 "
       f"--json {QC}/{srr}.fastp.json --html {QC}/{srr}.fastp.html", srr)


# adapters fastp detected (Nextera transposase end, the Tn5 library prep): mate 1 / mate 2 share the first 19 bases
NEXTERA = {1: "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", 2: "CTGTCTCTTATACACATCTGACGCTGCCGACGA"}


def adapter_check(srr, n_pairs=600_000, motif_len=12, polyg_len=12) -> pd.DataFrame:
    """Evidence that what fastp cuts is the library adapter, and what happens to poly-G tails (first n_pairs read pairs).
    Per mate: share of reads with the adapter's first 12 bases (and how many are expected by chance), whether the bases after
    them continue this mate's own adapter, and the fate of reads with >= 12 G in a row after adapter trimming ALONE (no --trim_poly_g, unlike trim()): this shows what
    the adapter cut leaves behind, the reason poly-G trimming is switched on in the pipeline."""
    tmp = TMP / f"adcheck_{srr}"
    tmp.mkdir(exist_ok=True)
    sh(f"fastq-dump -X {n_pairs} --split-files -O {tmp} {SRA / srr / (srr + '.sra')}", srr)
    sh(f"fastp -i {tmp}/{srr}_1.fastq -I {tmp}/{srr}_2.fastq -o {tmp}/t_1.fq -O {tmp}/t_2.fq "
       f"--detect_adapter_for_pe --length_required 20 --thread 8 --json /dev/null --html /dev/null", srr)

    def read(fn):
        out = {}
        with open(fn) as f:
            for name in f:
                seq = f.readline().rstrip(); f.readline(); f.readline()
                out[name.split()[0]] = seq
        return out

    mism = lambda a, b: sum(x != y for x, y in zip(a, b))
    G, rows = "G" * polyg_len, []
    for m in (1, 2):
        raw, kept = read(f"{tmp}/{srr}_{m}.fastq"), read(f"{tmp}/t_{m}.fq")
        own, other = NEXTERA[m], NEXTERA[3 - m]
        motif = own[:motif_len]
        hit = [n for n, q in raw.items() if motif in q]
        cls = {"own": 0, "other": 0, "neither": 0}
        for n in hit:
            q = raw[n]; i = q.find(motif); tail = q[i + 19:i + 33]     # the two mates' adapters are identical for 19 bases, then differ
            if len(tail) < 8:
                continue
            a, b = mism(tail, own[19:19 + len(tail)]), mism(tail, other[19:19 + len(tail)])
            cls["own" if a <= 1 < b else "other" if b <= 1 < a else "neither"] += 1
        pg = [n for n, q in raw.items() if G in q]
        after_ad = [n for n in pg if motif in raw[n] and raw[n].index(motif) < raw[n].index(G)]
        dropped = sum(n not in kept for n in pg)
        still = sum(n in kept and G in kept[n] for n in pg)
        pct = lambda x, d: round(100 * x / max(d, 1), 1)
        rows.append({"mate": m, "reads": len(raw), "adapter_start_%": pct(len(hit), len(raw)),
                     "expected_by_chance_reads": round(len(raw) * (150 - motif_len + 1) / 4 ** motif_len, 1),
                     "continues_own_adapter_%": pct(cls["own"], sum(cls.values())),
                     "matches_other_mate_adapter_%": pct(cls["other"], sum(cls.values())),
                     "polyG_reads_%": pct(len(pg), len(raw)), "polyG_after_adapter_%": pct(len(after_ad), len(pg)),
                     "polyG_dropped_%": pct(dropped, len(pg)), "polyG_removed_%": pct(len(pg) - dropped - still, len(pg)),
                     "polyG_remaining_%": pct(still, len(pg))})
    shutil.rmtree(tmp)
    return pd.DataFrame(rows).set_index("mate")


@timed("kallisto")
def quant(srr, trimmed=False):
    """Pseudoalign reads (raw or fastp-trimmed) to the transcriptome."""
    if has_quant(srr):
        return "skip"
    r1, r2 = (f"{FASTQ}/{srr}_trim_1.fastq.gz", f"{FASTQ}/{srr}_trim_2.fastq.gz") if trimmed \
        else (f"{FASTQ}/{srr}_1.fastq", f"{FASTQ}/{srr}_2.fastq")
    sh(f"kallisto quant -i {INDEX} -o {KALLISTO / srr} --threads 8 {r1} {r2}", srr)


def process(srr, trim_reads=False, keep_fastq=False):
    """sra -> FASTQ -> [fastp] -> kallisto -> delete FASTQ (unless keep_fastq). Skips finished runs."""
    if has_quant(srr):
        return
    fastq(srr)
    if trim_reads:
        trim(srr)
    quant(srr, trimmed=trim_reads)
    if not keep_fastq:
        for f in FASTQ.glob(f"{srr}_*fastq*"):
            f.unlink()


def build_index():
    if INDEX.is_file():
        return print("index exists")
    sh(f"kallisto index -i {INDEX} {TRANSCRIPTS}", "kallisto_index")


def clean_tmp():
    shutil.rmtree(TMP, ignore_errors=True)
    TMP.mkdir(parents=True, exist_ok=True)


def download_reference():
    """GENCODE vM25 mouse transcriptome (kallisto) and annotation (tx -> gene map for tximport)."""
    base = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25"
    REF.mkdir(exist_ok=True)
    for f in ("gencode.vM25.transcripts.fa.gz", "gencode.vM25.annotation.gtf.gz"):
        if not (REF / f).is_file():
            sh(f"curl -sL -o {REF / f} {base}/{f}", "reference")


def qc_summary(runs: pd.DataFrame) -> pd.DataFrame:
    """Per-run kallisto pseudoalignment rate (+ fastp read pairs/Q30/duplication if trimmed)."""
    import json
    rows = {}
    for s in runs.index:
        kj, fj = KALLISTO / s / "run_info.json", QC / f"{s}.fastp.json"
        if not kj.is_file():
            continue
        k = json.load(open(kj))
        rows[s] = {"pairs_M": k["n_processed"] / 1e6, "aligned_%": k["p_pseudoaligned"]}
        if fj.is_file():
            f = json.load(open(fj))
            rows[s] |= {"q30_%": 100 * f["summary"]["after_filtering"]["q30_rate"], "dup_%": 100 * f["duplication"]["rate"]}
    if not rows:
        return pd.DataFrame(columns=["pairs_M", "aligned_%", "group"])
    return pd.DataFrame(rows).T.round(1).join(runs["group"])


# ---- results copied back from Euler ----------------------------------------------------------
EULER = ROOT / "results/euler"      # put the tarballs from euler/04_collect.sh here (git-ignored)


def unpack_euler(src: Path = EULER):
    """Unpack the tarballs from euler/04_collect.sh: kallisto_results and qc_results into results/ (-> results/kallisto, results/qc),
    star_results into results/star/ (-> results/star/mapping). Missing tarballs are reported and skipped."""
    import tarfile
    targets = {"kallisto_results.tar.gz": ROOT / "results", "qc_results.tar.gz": ROOT / "results", "star_results.tar.gz": ROOT / "results/star"}
    for name, dest in targets.items():
        f = next((p for p in (src / name, ROOT / "results/star" / name) if p.is_file()), None)    # star tarball may sit in results/star/
        if f is None:
            print(f"not found: {name}")
            continue
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(f) as t:
            t.extractall(dest, filter="data")
        print(f"unpacked {name} -> {dest.relative_to(ROOT)}")


# ---- STAR counts and mapping statistics (from the unpacked Euler results) ----------------------------------------------------
STAR = ROOT / "results/star"


def load_star(runs: pd.DataFrame, root: Path = STAR):
    """Parse results/star/mapping/<SRR>/{ReadsPerGene.out.tab,Log.final.out} (from euler/04_collect.sh).
    Returns (counts genes x runs [unstranded column], qc table with mapping rates and a strandedness check)."""
    import re
    counts, qc = {}, {}
    for s in runs.index:
        d = root / "mapping" / s
        t = pd.read_csv(d / "ReadsPerGene.out.tab", sep="\t", header=None, index_col=0)
        g = t.iloc[4:]                                      # first 4 rows are N_unmapped / multimapping / noFeature / ambiguous
        counts[s] = g[1]                                    # col 1 = unstranded (Smart-seq2 libraries are unstranded)
        log = (d / "Log.final.out").read_text()
        pick = lambda key: float(re.search(rf"{re.escape(key)}\s*\|\s*([\d.]+)", log).group(1))
        qc[s] = {"input_pairs_M": pick("Number of input reads") / 1e6,
                 "unique_%": pick("Uniquely mapped reads %"),
                 "multi_%": pick("% of reads mapped to multiple loci"),
                 "assigned_to_gene_%": 100 * g[1].sum() / (pick("Number of input reads")),
                 "frac_read1_same_strand": g[2].sum() / (g[2].sum() + g[3].sum())}   # ~0.5 = unstranded
    counts = pd.DataFrame(counts).astype(int)
    counts.index.name = "gene_id"
    qc = pd.DataFrame(qc).T.round(2).join(runs["group"])
    qc.index.name = "SRR"
    return counts, qc
