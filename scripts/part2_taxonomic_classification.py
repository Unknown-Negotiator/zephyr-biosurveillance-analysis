#!/usr/bin/env python3
"""Part 2 Zephyr respiratory-virus classification and coverage analysis.

This script downloads a recent set of Zephyr respiratory-read FASTA files,
maps reads to a curated panel of common respiratory-virus references with
minimap2, and summarizes virus calls plus approximate genome coverage.
"""

from __future__ import annotations

import gzip
import html
import io
import re
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
META_DIR = DATA_DIR / "metadata"
READS_DIR = DATA_DIR / "reads"
REF_DIR = DATA_DIR / "references"
MAP_DIR = ROOT / "results" / "mappings"
RESULTS_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures" / "coverage"

ZEPHYR_URL = "https://data.securebio.org/zephyr/"
READ_BASE_URL = "https://data.securebio.org/zephyr/"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

N_SAMPLES = 12
MIN_TABLE_READS = 200
MIN_ALIGNMENT_IDENTITY = 0.70
MIN_ALIGNMENT_LENGTH = 100

PYTHON = shutil.which("python") or "python"
MINIMAP2 = MINIMAP2 = shutil.which("minimap2") or "minimap2"


@dataclass(frozen=True)
class Reference:
    accession: str
    virus: str
    segment: str


REFERENCES = [
    Reference("NC_045512", "SARS-CoV-2", "genome"),
    Reference("NC_001617", "Rhinovirus A", "genome"),
    Reference("NC_001490", "Rhinovirus B", "genome"),
    Reference("NC_009996", "Rhinovirus C", "genome"),
    Reference("NC_002016", "Influenza A", "segment_1_PB2"),
    Reference("NC_002017", "Influenza A", "segment_2_PB1"),
    Reference("NC_002018", "Influenza A", "segment_3_PA"),
    Reference("NC_002019", "Influenza A", "segment_4_HA"),
    Reference("NC_002020", "Influenza A", "segment_5_NP"),
    Reference("NC_002021", "Influenza A", "segment_6_NA"),
    Reference("NC_002022", "Influenza A", "segment_7_MP"),
    Reference("NC_002023", "Influenza A", "segment_8_NS"),
    Reference("NC_002204", "Influenza B", "segment_1_PB1"),
    Reference("NC_002205", "Influenza B", "segment_2_PB2"),
    Reference("NC_002206", "Influenza B", "segment_3_PA"),
    Reference("NC_002207", "Influenza B", "segment_4_HA"),
    Reference("NC_002208", "Influenza B", "segment_5_NP"),
    Reference("NC_002209", "Influenza B", "segment_6_NA"),
    Reference("NC_002210", "Influenza B", "segment_7_MP"),
    Reference("NC_002211", "Influenza B", "segment_8_NS"),
    Reference("NC_006307", "Influenza C", "segment_1_PB2"),
    Reference("NC_006308", "Influenza C", "segment_2_PB1"),
    Reference("NC_006309", "Influenza C", "segment_3_P3"),
    Reference("NC_006310", "Influenza C", "segment_4_HEF"),
    Reference("NC_006311", "Influenza C", "segment_5_NP"),
    Reference("NC_006312", "Influenza C", "segment_6_MP"),
    Reference("NC_001803", "Respiratory syncytial virus A", "genome"),
    Reference("NC_001781", "Respiratory syncytial virus B", "genome"),
    Reference("NC_004148", "Human metapneumovirus", "genome"),
    Reference("NC_003461", "Human parainfluenza virus 1", "genome"),
    Reference("NC_003443", "Human parainfluenza virus 2", "genome"),
    Reference("NC_001796", "Human parainfluenza virus 3", "genome"),
    Reference("NC_002645", "Human coronavirus 229E", "genome"),
    Reference("NC_006213", "Human coronavirus OC43", "genome"),
    Reference("NC_005831", "Human coronavirus NL63", "genome"),
    Reference("NC_006577", "Human coronavirus HKU1", "genome"),
    Reference("NC_001405", "Human adenovirus C", "genome"),
    Reference("NC_007455", "Human bocavirus 1", "genome"),
]


def ensure_dirs() -> None:
    for path in [META_DIR, READS_DIR, REF_DIR, MAP_DIR, RESULTS_DIR, FIG_DIR, ROOT / ".cache" / "matplotlib"]:
        path.mkdir(parents=True, exist_ok=True)


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    tmp = path.with_suffix(path.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 biodef-part2"})
    with urllib.request.urlopen(request, timeout=120) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(path)


def strip_tags(value: str) -> str:
    value = re.sub(r"<[^>]+>", "", value)
    return html.unescape(value).strip()


def parse_read_count(text: str) -> int:
    match = re.search(r"([0-9,]+)\s+read", text)
    if not match:
        return 0
    return int(match.group(1).replace(",", ""))


def parse_zephyr_table(index_html: Path) -> pd.DataFrame:
    text = index_html.read_text()
    rows = []
    for row_html in re.findall(r"<tr[^>]*>(.*?)</tr>", text, flags=re.S):
        cells = re.findall(r"<td[^>]*>(.*?)</td>", row_html, flags=re.S)
        if len(cells) < 4:
            continue
        link_match = re.search(r'href="([^"]+\.respiratory\.fasta\.gz)"', cells[3])
        read_count = parse_read_count(cells[3])
        if not link_match or read_count == 0:
            continue
        file_path = html.unescape(link_match.group(1))
        sample_id = Path(file_path).name.replace(".respiratory.fasta.gz", "")
        rows.append(
            {
                "date": strip_tags(cells[0]),
                "location": strip_tags(cells[1]),
                "participant_samples": int(strip_tags(cells[2]).replace(",", "")),
                "published_respiratory_viral_reads": read_count,
                "relative_path": file_path,
                "sample_id": sample_id,
                "url": urllib.parse.urljoin(READ_BASE_URL, file_path),
            }
        )
    if not rows:
        raise RuntimeError("No respiratory FASTA links found in Zephyr table.")
    return pd.DataFrame(rows)


def count_fasta_reads(path: Path) -> tuple[int, int, float]:
    read_count = 0
    total_bases = 0
    lengths = []
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            length = len(record.seq)
            read_count += 1
            total_bases += length
            lengths.append(length)
    return read_count, total_bases, float(median(lengths)) if lengths else 0.0


def ncbi_fetch(accession: str) -> str:
    query = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": accession,
            "rettype": "fasta",
            "retmode": "text",
        }
    )
    url = f"{NCBI_EFETCH}?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "biodef-part2/0.1"})
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                data = response.read().decode()
            break
        except Exception as exc:
            last_error = exc
            time.sleep(1 + attempt)
    else:
        raise RuntimeError(f"NCBI fetch failed for {accession}: {last_error}")
    if not data.startswith(">"):
        raise RuntimeError(f"NCBI did not return FASTA for {accession}: {data[:120]}")
    return data


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def description_matches_virus(description: str, virus: str) -> bool:
    desc = description.lower()
    virus_lower = virus.lower()
    synonym_groups = {
        "sars-cov-2": ["sars-cov-2", "severe acute respiratory syndrome coronavirus 2"],
        "respiratory syncytial virus a": ["respiratory syncytial virus", "orthopneumovirus"],
        "respiratory syncytial virus b": ["respiratory syncytial virus", "orthopneumovirus"],
        "human metapneumovirus": ["human metapneumovirus"],
        "human parainfluenza virus 1": ["parainfluenza virus 1"],
        "human parainfluenza virus 2": ["rubulavirus 2", "parainfluenza virus 2"],
        "human parainfluenza virus 3": ["parainfluenza virus 3"],
        "human coronavirus 229e": ["coronavirus 229e"],
        "human coronavirus oc43": ["coronavirus oc43"],
        "human coronavirus nl63": ["coronavirus nl63"],
        "human coronavirus hku1": ["coronavirus hku1"],
        "human adenovirus c": ["adenovirus c"],
        "human bocavirus 1": ["bocavirus", "bocaparvovirus"],
    }
    if virus_lower in synonym_groups:
        return any(term in desc for term in synonym_groups[virus_lower])
    terms = [term for term in virus_lower.replace("human ", "").split() if term]
    return all(term in desc for term in terms[:2])


def build_references() -> pd.DataFrame:
    fasta_path = REF_DIR / "respiratory_reference_panel.fasta"
    meta_path = REF_DIR / "respiratory_reference_panel.csv"
    if fasta_path.exists() and meta_path.exists():
        meta = pd.read_csv(meta_path)
        expected_accessions = {ref.accession for ref in REFERENCES}
        found_accessions = set(meta["accession"])
        if expected_accessions.issubset(found_accessions):
            return meta
        missing = ", ".join(sorted(expected_accessions - found_accessions))
        print(f"Rebuilding reference panel; existing metadata is missing: {missing}", flush=True)

    tmp_fasta_path = fasta_path.with_suffix(".fasta.tmp")
    tmp_meta_path = meta_path.with_suffix(".csv.tmp")
    metadata = []
    skipped = []
    with tmp_fasta_path.open("w") as out:
        for ref in REFERENCES:
            print(f"Fetching reference {ref.accession} ({ref.virus}, {ref.segment})", flush=True)
            try:
                fasta = ncbi_fetch(ref.accession)
            except RuntimeError as exc:
                skipped.append({"accession": ref.accession, "virus": ref.virus, "segment": ref.segment, "reason": str(exc)})
                print(f"Skipping {ref.accession}: {exc}", flush=True)
                continue
            record = next(SeqIO.parse(io.StringIO(fasta), "fasta"))
            if not description_matches_virus(record.description, ref.virus):
                skipped.append(
                    {
                        "accession": ref.accession,
                        "virus": ref.virus,
                        "segment": ref.segment,
                        "reason": f"description mismatch: {record.description}",
                    }
                )
                print(f"Skipping {ref.accession}: description mismatch ({record.description})", flush=True)
                continue
            target_id = f"{safe_id(ref.virus)}|{safe_id(ref.segment)}|{ref.accession}"
            sequence = str(record.seq).upper()
            out.write(f">{target_id} {record.description}\n")
            for i in range(0, len(sequence), 80):
                out.write(sequence[i : i + 80] + "\n")
            metadata.append(
                {
                    "target_id": target_id,
                    "accession": ref.accession,
                    "virus": ref.virus,
                    "segment": ref.segment,
                    "length": len(sequence),
                    "ncbi_description": record.description,
                }
            )
            time.sleep(0.34)
    meta = pd.DataFrame(metadata)
    if meta.empty:
        raise RuntimeError("No references were successfully fetched.")
    meta.to_csv(tmp_meta_path, index=False)
    tmp_fasta_path.replace(fasta_path)
    tmp_meta_path.replace(meta_path)
    pd.DataFrame(skipped).to_csv(REF_DIR / "skipped_references.csv", index=False)
    return meta


def run_minimap(sample_id: str, read_path: Path, ref_fasta: Path) -> Path:
    paf_path = MAP_DIR / f"{sample_id}.paf"
    if paf_path.exists() and paf_path.stat().st_size > 0 and paf_path.stat().st_mtime >= ref_fasta.stat().st_mtime:
        return paf_path
    cmd = [
        str(MINIMAP2),
        "-t",
        "4",
        "-x",
        "map-ont",
        "--secondary=no",
        "-c",
        str(ref_fasta),
        str(read_path),
    ]
    with paf_path.open("w") as out:
        subprocess.run(cmd, check=True, stdout=out)
    return paf_path


def parse_paf(paf_path: Path, ref_meta: pd.DataFrame) -> pd.DataFrame:
    target_to_virus = ref_meta.set_index("target_id")["virus"].to_dict()
    rows = []
    with paf_path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qname = parts[0]
            qlen = int(parts[1])
            tname = parts[5]
            tlen = int(parts[6])
            tstart = int(parts[7])
            tend = int(parts[8])
            matches = int(parts[9])
            aln_len = int(parts[10])
            mapq = int(parts[11])
            identity = matches / aln_len if aln_len else 0.0
            if identity < MIN_ALIGNMENT_IDENTITY or aln_len < MIN_ALIGNMENT_LENGTH:
                continue
            rows.append(
                {
                    "read_id": qname,
                    "read_length": qlen,
                    "target_id": tname,
                    "target_length": tlen,
                    "target_start": tstart,
                    "target_end": tend,
                    "aligned_bases": aln_len,
                    "matches": matches,
                    "identity": identity,
                    "mapq": mapq,
                    "virus": target_to_virus.get(tname, "unknown"),
                }
            )
    if not rows:
        return pd.DataFrame(
            columns=[
                "read_id",
                "read_length",
                "target_id",
                "target_length",
                "target_start",
                "target_end",
                "aligned_bases",
                "matches",
                "identity",
                "mapq",
                "virus",
            ]
        )
    df = pd.DataFrame(rows)
    df = df.sort_values(["read_id", "matches", "identity", "mapq", "aligned_bases"], ascending=[True, False, False, False, False])
    return df.drop_duplicates("read_id", keep="first")


def interval_coverage(intervals: list[tuple[int, int]], length: int) -> tuple[int, float]:
    if not intervals or length <= 0:
        return 0, 0.0
    diff = np.zeros(length + 1, dtype=np.int32)
    for start, end in intervals:
        start = max(0, min(length, start))
        end = max(0, min(length, end))
        if end <= start:
            continue
        diff[start] += 1
        diff[end] -= 1
    depth = np.cumsum(diff[:-1])
    covered = int((depth > 0).sum())
    mean_depth = float(depth.mean())
    return covered, mean_depth


def confidence_label(reads: int, breadth: float, median_identity_value: float) -> str:
    if reads >= 50 and breadth >= 0.20 and median_identity_value >= 0.85:
        return "high"
    if reads >= 10 and breadth >= 0.02 and median_identity_value >= 0.75:
        return "medium"
    if reads >= 3 and median_identity_value >= 0.70:
        return "low"
    return "very low"


def summarize_sample(sample: pd.Series, alignments: pd.DataFrame, ref_meta: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    call_rows = []
    coverage_rows = []
    total_reads = int(sample["actual_read_count"])
    ref_lengths_by_virus = ref_meta.groupby("virus")["length"].sum().to_dict()

    for virus, virus_df in alignments.groupby("virus"):
        reads = int(virus_df["read_id"].nunique())
        ref_bases = int(ref_lengths_by_virus[virus])
        covered_bases = 0
        depth_weighted_sum = 0.0
        for target_id, target_df in virus_df.groupby("target_id"):
            length = int(ref_meta.loc[ref_meta["target_id"] == target_id, "length"].iloc[0])
            intervals = list(zip(target_df["target_start"].astype(int), target_df["target_end"].astype(int)))
            covered, mean_depth = interval_coverage(intervals, length)
            covered_bases += covered
            depth_weighted_sum += mean_depth * length
            coverage_rows.append(
                {
                    "sample_id": sample["sample_id"],
                    "virus": virus,
                    "target_id": target_id,
                    "covered_bases": covered,
                    "reference_bases": length,
                    "coverage_breadth": covered / length if length else 0,
                    "mean_depth": mean_depth,
                    "reads_assigned": int(target_df["read_id"].nunique()),
                }
            )
        breadth = covered_bases / ref_bases if ref_bases else 0.0
        mean_depth = depth_weighted_sum / ref_bases if ref_bases else 0.0
        med_identity = float(virus_df["identity"].median())
        call_rows.append(
            {
                "sample_id": sample["sample_id"],
                "date": sample["date"],
                "location": sample["location"],
                "virus": virus,
                "reads_assigned": reads,
                "percent_of_sample_reads": 100 * reads / total_reads if total_reads else 0,
                "median_identity": med_identity,
                "mean_identity": float(virus_df["identity"].mean()),
                "genome_coverage_breadth": breadth,
                "mean_depth": mean_depth,
                "confidence": confidence_label(reads, breadth, med_identity),
            }
        )
    return call_rows, coverage_rows


def plot_coverage(sample_id: str, virus: str, alignments: pd.DataFrame, ref_meta: pd.DataFrame) -> Path:
    virus_refs = ref_meta[ref_meta["virus"] == virus].copy()
    n = len(virus_refs)
    fig, axes = plt.subplots(n, 1, figsize=(10, max(2.2, 1.35 * n)), sharex=False)
    if n == 1:
        axes = [axes]
    for ax, (_, ref) in zip(axes, virus_refs.iterrows()):
        length = int(ref["length"])
        diff = np.zeros(length + 1, dtype=np.int32)
        subset = alignments[alignments["target_id"] == ref["target_id"]]
        for _, row in subset.iterrows():
            start = max(0, min(length, int(row["target_start"])))
            end = max(0, min(length, int(row["target_end"])))
            if end > start:
                diff[start] += 1
                diff[end] -= 1
        depth = np.cumsum(diff[:-1])
        x = np.arange(1, length + 1)
        ax.fill_between(x, depth, step="pre", color="#2f6f9f", alpha=0.85)
        ax.set_xlim(1, length)
        ax.set_ylabel(ref["segment"].replace("_", " "), fontsize=8)
        ax.grid(axis="y", alpha=0.25)
    axes[-1].set_xlabel("Reference position (nt)")
    fig.suptitle(f"{sample_id}: {virus} read coverage", fontsize=12)
    fig.tight_layout()
    out = FIG_DIR / f"{sample_id}__{safe_id(virus)}.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def main() -> None:
    ensure_dirs()
    index_path = META_DIR / "zephyr_index.html"
    download(ZEPHYR_URL, index_path)

    all_samples = parse_zephyr_table(index_path)
    all_samples.to_csv(META_DIR / "zephyr_respiratory_read_files.csv", index=False)
    selected = all_samples[all_samples["published_respiratory_viral_reads"] >= MIN_TABLE_READS].head(N_SAMPLES).copy()
    selected.to_csv(META_DIR / "selected_samples.csv", index=False)

    read_stats = []
    for _, sample in selected.iterrows():
        read_path = READS_DIR / Path(sample["relative_path"]).name
        download(sample["url"], read_path)
        count, bases, med_len = count_fasta_reads(read_path)
        read_stats.append(
            {
                "sample_id": sample["sample_id"],
                "read_file": str(read_path.relative_to(ROOT)),
                "actual_read_count": count,
                "total_bases": bases,
                "median_read_length": med_len,
            }
        )
    selected = selected.merge(pd.DataFrame(read_stats), on="sample_id")
    selected.to_csv(META_DIR / "selected_samples_with_read_stats.csv", index=False)

    ref_meta = build_references()
    ref_fasta = REF_DIR / "respiratory_reference_panel.fasta"

    all_call_rows = []
    all_coverage_rows = []
    plot_rows = []
    for _, sample in selected.iterrows():
        read_path = READS_DIR / Path(sample["relative_path"]).name
        paf = run_minimap(sample["sample_id"], read_path, ref_fasta)
        alignments = parse_paf(paf, ref_meta)
        alignments.to_csv(MAP_DIR / f"{sample['sample_id']}.best_alignments.csv", index=False)
        call_rows, coverage_rows = summarize_sample(sample, alignments, ref_meta)
        all_call_rows.extend(call_rows)
        all_coverage_rows.extend(coverage_rows)
        for call in call_rows:
            if call["confidence"] in {"high", "medium"} or call["reads_assigned"] >= 25:
                virus_alignments = alignments[alignments["virus"] == call["virus"]]
                plot_path = plot_coverage(sample["sample_id"], call["virus"], virus_alignments, ref_meta)
                plot_rows.append(
                    {
                        "sample_id": sample["sample_id"],
                        "virus": call["virus"],
                        "figure": str(plot_path.relative_to(ROOT)),
                    }
                )

    calls = pd.DataFrame(all_call_rows)
    if not calls.empty:
        calls = calls.sort_values(["sample_id", "reads_assigned"], ascending=[True, False])
    calls.to_csv(RESULTS_DIR / "virus_calls_by_sample.csv", index=False)

    coverage = pd.DataFrame(all_coverage_rows)
    coverage.to_csv(RESULTS_DIR / "coverage_by_reference_segment.csv", index=False)

    if calls.empty:
        summary = pd.DataFrame()
    else:
        summary = (
            calls.groupby("virus")
            .agg(
                positive_pools=("sample_id", "nunique"),
                total_reads_assigned=("reads_assigned", "sum"),
                max_reads_in_pool=("reads_assigned", "max"),
                median_identity=("median_identity", "median"),
                max_genome_coverage_breadth=("genome_coverage_breadth", "max"),
                median_genome_coverage_breadth=("genome_coverage_breadth", "median"),
            )
            .reset_index()
            .sort_values(["positive_pools", "total_reads_assigned"], ascending=[False, False])
        )
    summary.to_csv(RESULTS_DIR / "virus_summary_across_samples.csv", index=False)
    pd.DataFrame(plot_rows).to_csv(RESULTS_DIR / "coverage_figures.csv", index=False)

    print(f"Selected samples: {len(selected)}")
    print(f"Reference records: {len(ref_meta)}")
    print(f"Virus calls: {len(calls)}")
    print(f"Wrote: {RESULTS_DIR / 'virus_calls_by_sample.csv'}")
    print(f"Wrote: {RESULTS_DIR / 'virus_summary_across_samples.csv'}")


if __name__ == "__main__":
    main()
