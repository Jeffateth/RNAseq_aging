"""Python helpers for notebooks 01 and 02.

Contents
  Setup and paths ..... ROOT, SRA, TMP, QC, KALLISTO, STAR, RSEM; load_runs() reads metadata/samples.csv
  Euler results ....... unpack_euler (tarballs copied back from the cluster), load_star (STAR counts + mapping statistics)
  Checks .............. status, check_consistent, qc_summary (kallisto + fastp), check_rsem
  QC reports .......... multiqc (combined FastQC report), fastqc_summary (FastQC verdicts), fastqc_metrics (the numbers behind them); raw or trimmed reads
  Adapter check ....... download_run (one .sra), adapter_check (evidence that fastp cuts the library adapter), tail_check (poly-A / poly-G runs before and after trimming), tail_anatomy (where the parts of a poly-G read lie)
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


def _fastqc_zip(s, m, which):
    """Path of the FastQC report of run s, mate m: `which` = "raw" (reads as sequenced) or "trimmed" (after fastp)."""
    return QC / which / (f"{s}_{m}_fastqc.zip" if which == "raw" else f"{s}_trim_{m}_fastqc.zip")


def multiqc(which="raw"):
    """One combined report for all FastQC results of `which` ("raw" or "trimmed") -> results/qc/multiqc_<which>.html"""
    sh(f"multiqc -q -f -o {QC} -n multiqc_{which} {QC / which}", "multiqc")
    return QC / f"multiqc_{which}.html"


KEY_MODULES = ["Per base sequence quality", "Per sequence quality scores", "Adapter Content",
               "Overrepresented sequences", "Sequence Duplication Levels"]


def fastqc_summary(srrs, modules=KEY_MODULES, which="raw") -> pd.DataFrame:
    """FastQC PASS/WARN/FAIL for the key modules (rows: run_mate)."""
    import zipfile
    rows = {}
    for s in srrs:
        for m in (1, 2):
            z = _fastqc_zip(s, m, which)
            if z.is_file():
                with zipfile.ZipFile(z) as zf:
                    txt = zf.read(next(n for n in zf.namelist() if n.endswith("summary.txt"))).decode()
                rows[f"{s}_{m}"] = {l.split("\t")[1]: l.split("\t")[0] for l in txt.splitlines()}
    return pd.DataFrame(rows).T[modules] if rows else pd.DataFrame(columns=modules)


def fastqc_metrics(srrs, which="raw") -> pd.DataFrame:
    """The numbers behind FastQC's verdicts, one row per read file (run_mate), from the reports of `which` ("raw" or "trimmed").
    base_median_min / base_lowerq_min / base_mean_min: lowest median / lower quartile / mean base quality over read positions (Phred);
    seq_mean_q: most common mean quality per read (Phred); nextera_max, polyA_max, polyG_max: highest share of reads (%) with that adapter-table
    entry at any position; overrep_top: share of reads (%) of the most common overrepresented sequence; dedup_pct: % of reads remaining after
    deduplication (FastQC's estimate from the first 100,000 reads; depends on file size)."""
    import zipfile
    rows = {}
    for s in srrs:
        for m in (1, 2):
            z = _fastqc_zip(s, m, which)
            if not z.is_file():
                continue
            with zipfile.ZipFile(z) as zf:
                t = zf.read(next(n for n in zf.namelist() if n.endswith("fastqc_data.txt"))).decode().split("\n")

            def section(title):
                out, hdr, on = [], None, False
                for line in t:
                    if line.startswith(">>" + title):
                        on = True
                    elif on and line.startswith(">>END_MODULE"):
                        break
                    elif on and line.startswith("#"):
                        hdr = line.lstrip("#").split("\t")
                    elif on:
                        out.append(line.split("\t"))
                return hdr, out
            _, pb = section("Per base sequence quality")
            _, ps = section("Per sequence quality scores")
            hdr, ad = section("Adapter Content")
            _, ov = section("Overrepresented sequences")
            adapt = {c: max(float(r[j]) for r in ad) for j, c in enumerate(hdr) if j}
            rows[f"{s}_{m}"] = {
                "base_median_min": min(float(r[2]) for r in pb), "base_lowerq_min": min(float(r[3]) for r in pb),
                "base_mean_min": round(min(float(r[1]) for r in pb), 2),
                "seq_mean_q": float(max(ps, key=lambda r: float(r[1]))[0]),
                "nextera_max": adapt["Nextera Transposase Sequence"], "polyA_max": adapt["PolyA"], "polyG_max": adapt["PolyG"],
                "overrep_top": max((float(r[2]) for r in ov), default=0.0),
                "dedup_pct": round(float(next(l for l in t if l.startswith("#Total Deduplicated Percentage")).split("\t")[1]), 1)}
    return pd.DataFrame(rows).T


NEXTERA = {1: "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC", 2: "CTGTCTCTTATACACATCTGACGCTGCCGACGA"}


def adapter_check(srr, n_pairs=600_000, motif_len=12, polyg_len=10) -> pd.DataFrame:
    """Evidence that what fastp cuts is the library adapter, and what happens to poly-G tails (first n_pairs read pairs).
    Per mate: share of reads with the adapter's first 12 bases (and how many are expected by chance), whether the bases after
    them continue this mate's own adapter, and the fate of poly-G reads (a tail of >= `polyg_len` G at the read end, fastp's own rule) after adapter trimming
    ALONE (no --trim_poly_g, unlike the pipeline): this shows what the adapter cut leaves behind, the reason poly-G trimming is switched on."""
    import re
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
    g_tail, rows = re.compile(f"G{{{polyg_len},}}$"), []
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
        pg = [n for n, q in raw.items() if g_tail.search(q)]
        after_ad = [n for n in pg if motif in raw[n] and raw[n].index(motif) < g_tail.search(raw[n]).start()]
        dropped = sum(n not in kept for n in pg)
        still = sum(n in kept and bool(g_tail.search(kept[n])) for n in pg)
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


def tail_anatomy(srr, n_pairs=600_000, motif_len=12, polyg_len=10) -> pd.DataFrame:
    """Where the parts of a read with a poly-G tail lie (first n_pairs read pairs, raw reads, only reads whose adapter comes before the G tail).
    Per mate, median [middle half] in bases: `insert` = bases before the adapter, `adapter_to_tail` = adapter start to G-tail start (adapter +
    barcode + rest of the library construct), `G_tail` = length of the G tail; `barcode_share_%` = reads with the commonest 8 bases right after the adapter."""
    import re
    import statistics as st
    from collections import Counter
    tmp = TMP / f"anat_{srr}"
    tmp.mkdir(exist_ok=True)
    sh(f"fastq-dump -X {n_pairs} --split-files -O {tmp} {SRA / srr / (srr + '.sra')}", srr)
    g_tail, rows = re.compile(f"G{{{polyg_len},}}$"), {}
    for m in (1, 2):
        motif, alen = NEXTERA[m][:motif_len], len(NEXTERA[m])
        ins, span, tail, bar = [], [], [], Counter()
        with open(tmp / f"{srr}_{m}.fastq") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    q = line.strip(); g = g_tail.search(q)
                    a = q.find(motif)
                    if g and 0 <= a < g.start():
                        ins.append(a); span.append(g.start() - a); tail.append(len(q) - g.start())
                        if a + alen + 8 <= g.start():
                            bar[q[a + alen:a + alen + 8]] += 1
        iqr = lambda v: f"{st.median(v):.0f} [{st.quantiles(v, n=4)[0]:.0f}-{st.quantiles(v, n=4)[2]:.0f}]"
        rows[f"mate {m}"] = {"reads": len(ins), "adapter_length": alen, "insert": iqr(ins), "adapter_to_tail": iqr(span), "G_tail": iqr(tail),
                             "barcode_share_%": round(100 * bar.most_common(1)[0][1] / sum(bar.values()))}
    shutil.rmtree(tmp)
    return pd.DataFrame(rows).T


def tail_check(srr, n_pairs=600_000, run=12) -> pd.DataFrame:
    """Homopolymer runs in the first n_pairs read pairs, before and after fastp (same settings as the pipeline).
    Columns (% of reads): `A_run` / `T_run` = a run of >= `run` A / T anywhere (poly-T is poly-A on the other strand; the library is unstranded);
    `ends_in_A` = read ends in >= 10 A; `starts_with_T` = begins with >= 10 T; `ends_in_G` = read ends in >= 10 G (the poly-G tail fastp trims)."""
    import re
    tmp = TMP / f"polya_{srr}"
    tmp.mkdir(exist_ok=True)
    sh(f"fastq-dump -X {n_pairs} --split-files -O {tmp} {SRA / srr / (srr + '.sra')}", srr)
    sh(f"fastp -i {tmp}/{srr}_1.fastq -I {tmp}/{srr}_2.fastq -o {tmp}/t_1.fq -O {tmp}/t_2.fq "
       f"--detect_adapter_for_pe --trim_poly_g --length_required 20 --thread 4 --json /dev/null --html /dev/null", srr)
    pats = {"A_run": re.compile(f"A{{{run},}}"), "T_run": re.compile(f"T{{{run},}}"), "ends_in_A": re.compile("A{10,}$"), "starts_with_T": re.compile("^T{10,}"), "ends_in_G": re.compile("G{10,}$")}
    rows = {}
    for label, name in (("raw", f"{srr}_%d.fastq"), ("trimmed", "t_%d.fq")):
        for mate in (1, 2):
            n, hits = 0, dict.fromkeys(pats, 0)
            with open(tmp / (name % mate)) as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        n += 1
                        for k, rx in pats.items():
                            hits[k] += bool(rx.search(line.strip()))
            rows[f"{label} mate {mate}"] = {"reads": n} | {k: round(100 * v / n, 2) for k, v in hits.items()}
    shutil.rmtree(tmp)
    return pd.DataFrame(rows).T


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
    """Unpack the tarballs from euler/04_collect.sh: kallisto_results, rsem_results, qc_results and qc_trimmed_results into results/ (-> results/kallisto,
    results/rsem, results/qc, results/qc/trimmed), star_results into results/star/ (-> results/star/mapping). Missing tarballs are reported and skipped."""
    import tarfile
    targets = {"kallisto_results.tar.gz": ROOT / "results", "rsem_results.tar.gz": ROOT / "results", "qc_results.tar.gz": ROOT / "results",
               "qc_trimmed_results.tar.gz": ROOT / "results", "star_results.tar.gz": ROOT / "results/star"}
    for name, dest in targets.items():
        f = src / name
        if not f.is_file():
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
