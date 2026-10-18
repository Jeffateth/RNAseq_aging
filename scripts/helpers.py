"""Python helpers for notebooks 01 and 02.

Contents
  Setup and paths ..... ROOT, SRA, TMP, QC, KALLISTO, STAR, RSEM; load_runs() reads metadata/samples.csv
  Euler results ....... unpack_euler (tarballs copied back from the cluster), load_star (STAR counts + mapping statistics)
  Checks .............. status, check_consistent, qc_summary (kallisto + fastp), check_rsem
  QC reports .......... multiqc (combined FastQC report), fastqc_summary (FastQC verdicts)
  Adapter check ....... download_run (one .sra), adapter_check (evidence that fastp cuts the library adapter)
Terms: docs/glossary.md."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# fastp and multiqc live in the pixi env (not activated when VS Code starts the kernel); sra-tools come from Homebrew
os.environ["PATH"] = f"{Path(sys.executable).parent}:{os.environ['PATH']}:/opt/homebrew/bin"
ROOT = Path(__file__).resolve().parents[1]
SRA, TMP = ROOT / "raw_data/sra", ROOT / "raw_data/tmp"
QC, KALLISTO, LOGS = ROOT / "results/qc", ROOT / "results/kallisto", ROOT / "results/logs"
for d in (SRA, TMP, QC, KALLISTO, LOGS):
    d.mkdir(parents=True, exist_ok=True)

# ---- small helpers ---------------------------------------------------------------------


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


# ---- stage checks -----------------------------------------------------------------------


def has_sra(srr):  return (SRA / srr / f"{srr}.sra").is_file()


def has_quant(srr): return (KALLISTO / srr / "abundance.h5").is_file()


def was_trimmed(srr):
    """True/False = quantified from fastp-trimmed/raw reads (read back from kallisto's own record); None if not quantified."""
    import json
    f = KALLISTO / srr / "run_info.json"
    return "_trim_" in json.load(open(f)).get("call", "") if f.is_file() else None


def status(runs: pd.DataFrame) -> pd.DataFrame:
    """Per-run table: is kallisto's result present (quantified) and were the reads fastp-trimmed (trimmed)?"""
    return pd.DataFrame({
        "group": runs["group"],
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


# ---- QC reports and the adapter check ---------------------------------------------------


def download_run(srr: str):
    """prefetch one run from NCBI (resumes a partial file) and check it with vdb-validate; needs sra-tools (brew install sratoolkit)."""
    if not has_sra(srr):
        sh(f"prefetch {srr} --output-directory {SRA} --max-size 100G", srr)
    sh(f"vdb-validate {SRA / srr / (srr + '.sra')}", srr)
    print(f"{srr}.sra present and valid")


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
    """Unpack the tarballs from euler/04_collect.sh: kallisto_results, rsem_results and qc_results into results/ (-> results/kallisto,
    results/rsem, results/qc), star_results into results/star/ (-> results/star/mapping). Missing tarballs are reported and skipped."""
    import tarfile
    targets = {"kallisto_results.tar.gz": ROOT / "results", "rsem_results.tar.gz": ROOT / "results", "qc_results.tar.gz": ROOT / "results",
               "star_results.tar.gz": ROOT / "results/star"}
    for name, dest in targets.items():
        f = next((p for p in (src / name, ROOT / "results/star" / name) if p.is_file()), None)    # star tarball may sit in results/star/
        if f is None:
            print(f"not found: {name}")
            continue
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(f) as t:
            t.extractall(dest, filter="data")
        print(f"unpacked {name} -> {dest.relative_to(ROOT)}")


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
                 "STAR_assigned_to_gene_%": 100 * g[1].sum() / (pick("Number of input reads")),
                 "frac_read1_same_strand": g[2].sum() / (g[2].sum() + g[3].sum())}   # ~0.5 = unstranded
    counts = pd.DataFrame(counts).astype(int)
    counts.index.name = "gene_id"
    qc = pd.DataFrame(qc).T.round(2).join(runs["group"])
    qc.index.name = "SRR"
    return counts, qc


RSEM = ROOT / "results/rsem"


def check_rsem(runs: pd.DataFrame) -> pd.DataFrame:
    """One row per run from results/rsem/<SRR>/<SRR>.genes.results (euler/04_collect.sh): genes in the reference, genes with
    any reads, and total expected counts in millions. RSEM shares multi-mapping reads between genes, so counts are not integers."""
    rows = {}
    for s in runs.index:
        t = pd.read_csv(RSEM / s / f"{s}.genes.results", sep="\t", usecols=["gene_id", "expected_count"])
        rows[s] = {"genes": len(t), "genes_with_reads": int((t.expected_count > 0).sum()), "counts_M": round(t.expected_count.sum() / 1e6, 2)}
    return pd.DataFrame.from_dict(rows, orient="index").join(runs["group"])
