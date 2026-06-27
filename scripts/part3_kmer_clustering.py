#!/usr/bin/env python3
"""Part 3B lightweight k-mer embedding and clustering.

This embeds both respiratory-virus reference sequences and sample-derived read
sets from the Part 2 mappings. Each sample-virus read set is represented by
normalized 5-mer composition and compared with the reference panel.
"""

from __future__ import annotations

import gzip
from itertools import product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
REF_FASTA = ROOT / "data" / "references" / "respiratory_reference_panel.fasta"
REF_META = ROOT / "data" / "references" / "respiratory_reference_panel.csv"
READS_DIR = ROOT / "data" / "reads"
MAPPING_DIR = ROOT / "results" / "mappings"
CALLS_PATH = ROOT / "results" / "virus_calls_by_sample.csv"
RESULTS_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures"
K = 5
N_CLUSTERS = 7
RANDOM_SEED = 7
MIN_SAMPLE_VIRUS_READS = 25


FAMILY_BY_VIRUS = {
    "SARS-CoV-2": "Coronaviridae",
    "Human coronavirus 229E": "Coronaviridae",
    "Human coronavirus OC43": "Coronaviridae",
    "Human coronavirus NL63": "Coronaviridae",
    "Human coronavirus HKU1": "Coronaviridae",
    "Rhinovirus A": "Picornaviridae",
    "Rhinovirus B": "Picornaviridae",
    "Rhinovirus C": "Picornaviridae",
    "Influenza A": "Orthomyxoviridae",
    "Influenza B": "Orthomyxoviridae",
    "Influenza C": "Orthomyxoviridae",
    "Respiratory syncytial virus A": "Pneumoviridae",
    "Respiratory syncytial virus B": "Pneumoviridae",
    "Human metapneumovirus": "Pneumoviridae",
    "Human parainfluenza virus 1": "Paramyxoviridae",
    "Human parainfluenza virus 2": "Paramyxoviridae",
    "Human parainfluenza virus 3": "Paramyxoviridae",
    "Human adenovirus C": "Adenoviridae",
    "Human bocavirus 1": "Parvoviridae",
}


def safe_id(value: str) -> str:
    return "".join(char if char.isalnum() or char in "_.-" else "_" for char in value).strip("_")


def ensure_dirs() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)


def kmer_index() -> tuple[list[str], dict[str, int]]:
    alphabet = "ACGT"
    kmers = ["".join(chars) for chars in product(alphabet, repeat=K)]
    return kmers, {kmer: i for i, kmer in enumerate(kmers)}


def add_kmers(seq: str, counts: np.ndarray, index: dict[str, int]) -> int:
    observed = 0
    seq = seq.upper()
    for i in range(0, len(seq) - K + 1):
        kmer = seq[i : i + K]
        col = index.get(kmer)
        if col is None:
            continue
        counts[col] += 1
        observed += 1
    return observed


def normalize_counts(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    return counts / total if total else counts


def kmer_matrix(records: list[SeqIO.SeqRecord], index: dict[str, int]) -> np.ndarray:
    matrix = np.zeros((len(records), len(index)), dtype=np.float64)
    for row, record in enumerate(records):
        counts = np.zeros(len(index), dtype=np.float64)
        add_kmers(str(record.seq), counts, index)
        matrix[row] = normalize_counts(counts)
    return matrix


def sample_read_profiles(index: dict[str, int]) -> tuple[pd.DataFrame, np.ndarray]:
    calls = pd.read_csv(CALLS_PATH)
    selected_calls = calls[calls["reads_assigned"] >= MIN_SAMPLE_VIRUS_READS].copy()
    rows = []
    profiles = []

    for sample_id, sample_calls in selected_calls.groupby("sample_id", sort=True):
        mapping_path = MAPPING_DIR / f"{sample_id}.best_alignments.csv"
        read_path = READS_DIR / f"{sample_id}.respiratory.fasta.gz"
        if not mapping_path.exists() or not read_path.exists():
            continue

        alignments = pd.read_csv(mapping_path)
        viruses = set(sample_calls["virus"])
        subset = alignments[alignments["virus"].isin(viruses)]
        read_to_virus = subset.set_index("read_id")["virus"].to_dict()
        counts_by_virus = {virus: np.zeros(len(index), dtype=np.float64) for virus in viruses}
        reads_found = {virus: 0 for virus in viruses}
        bases_found = {virus: 0 for virus in viruses}

        with gzip.open(read_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                virus = read_to_virus.get(record.id)
                if virus is None:
                    continue
                seq = str(record.seq)
                observed = add_kmers(seq, counts_by_virus[virus], index)
                if observed:
                    reads_found[virus] += 1
                    bases_found[virus] += len(seq)

        for _, call in sample_calls.iterrows():
            virus = call["virus"]
            counts = counts_by_virus[virus]
            if counts.sum() == 0:
                continue
            rows.append(
                {
                    "record_id": f"{sample_id}|{safe_id(virus)}|reads",
                    "record_type": "sample_reads",
                    "sample_id": sample_id,
                    "virus": virus,
                    "family": FAMILY_BY_VIRUS.get(virus, "Unknown"),
                    "segment": "assigned_reads",
                    "accession": "",
                    "length": bases_found[virus],
                    "reads_assigned": int(call["reads_assigned"]),
                    "reads_found": reads_found[virus],
                    "median_identity": float(call["median_identity"]),
                    "genome_coverage_breadth": float(call["genome_coverage_breadth"]),
                    "confidence": call["confidence"],
                }
            )
            profiles.append(normalize_counts(counts))

    if not rows:
        return pd.DataFrame(), np.zeros((0, len(index)), dtype=np.float64)
    return pd.DataFrame(rows), np.vstack(profiles)


def svd_embedding(matrix: np.ndarray, n_components: int = 6) -> np.ndarray:
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    return centered @ vh[:n_components].T


def kmeans(points: np.ndarray, n_clusters: int, seed: int = RANDOM_SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    best_labels = None
    best_inertia = np.inf
    for _ in range(50):
        centers = points[rng.choice(len(points), size=n_clusters, replace=False)].copy()
        labels = np.zeros(len(points), dtype=int)
        for _ in range(100):
            distances = ((points[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
            new_labels = distances.argmin(axis=1)
            if np.array_equal(new_labels, labels):
                break
            labels = new_labels
            for cluster in range(n_clusters):
                if np.any(labels == cluster):
                    centers[cluster] = points[labels == cluster].mean(axis=0)
        inertia = ((points - centers[labels]) ** 2).sum()
        if inertia < best_inertia:
            best_inertia = inertia
            best_labels = labels.copy()
    if best_labels is None:
        raise RuntimeError("k-means did not produce labels")
    return best_labels


def cluster_summary(embedded: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cluster, group in embedded.groupby("cluster"):
        family_counts = group["family"].value_counts()
        majority_family = family_counts.index[0]
        reference_counts = group.loc[group["record_type"] == "reference", "family"].value_counts()
        reference_majority_family = reference_counts.index[0] if not reference_counts.empty else ""
        rows.append(
            {
                "cluster": cluster,
                "n_records": len(group),
                "n_reference_records": int((group["record_type"] == "reference").sum()),
                "n_sample_read_sets": int((group["record_type"] == "sample_reads").sum()),
                "majority_family": majority_family,
                "family_purity": family_counts.iloc[0] / len(group),
                "reference_majority_family": reference_majority_family,
                "families": "; ".join(f"{name}:{count}" for name, count in family_counts.items()),
                "viruses": "; ".join(sorted(group["virus"].unique())),
            }
        )
    return pd.DataFrame(rows).sort_values(["cluster"])


def sample_reference_comparison(embedded: pd.DataFrame) -> pd.DataFrame:
    refs = embedded[embedded["record_type"] == "reference"].copy()
    samples = embedded[embedded["record_type"] == "sample_reads"].copy()
    rows = []
    dims = [col for col in embedded.columns if col.startswith("pc")]
    ref_points = refs[dims].to_numpy()
    for _, sample in samples.iterrows():
        point = sample[dims].to_numpy(dtype=float)
        distances = ((ref_points - point) ** 2).sum(axis=1) ** 0.5
        nearest = refs.iloc[int(distances.argmin())]
        cluster_refs = refs[refs["cluster"] == sample["cluster"]]
        if cluster_refs.empty:
            cluster_reference_majority_family = ""
            cluster_matches_assigned_family = False
        else:
            cluster_reference_majority_family = cluster_refs["family"].value_counts().index[0]
            cluster_matches_assigned_family = cluster_reference_majority_family == sample["family"]
        rows.append(
            {
                "sample_id": sample["sample_id"],
                "assigned_virus": sample["virus"],
                "assigned_family": sample["family"],
                "reads_assigned": int(sample["reads_assigned"]),
                "reads_found": int(sample["reads_found"]),
                "median_identity": sample["median_identity"],
                "genome_coverage_breadth": sample["genome_coverage_breadth"],
                "confidence": sample["confidence"],
                "cluster": int(sample["cluster"]),
                "nearest_reference_virus": nearest["virus"],
                "nearest_reference_family": nearest["family"],
                "nearest_reference_segment": nearest["segment"],
                "nearest_reference_distance": float(distances.min()),
                "nearest_reference_matches_virus": nearest["virus"] == sample["virus"],
                "nearest_reference_matches_family": nearest["family"] == sample["family"],
                "cluster_reference_majority_family": cluster_reference_majority_family,
                "cluster_matches_assigned_family": cluster_matches_assigned_family,
            }
        )
    return pd.DataFrame(rows).sort_values(["sample_id", "assigned_virus"])


def plot_embedding(embedded: pd.DataFrame) -> Path:
    colors = {
        "Adenoviridae": "#4e79a7",
        "Coronaviridae": "#f28e2b",
        "Orthomyxoviridae": "#e15759",
        "Paramyxoviridae": "#76b7b2",
        "Parvoviridae": "#59a14f",
        "Picornaviridae": "#edc948",
        "Pneumoviridae": "#b07aa1",
    }
    fig, ax = plt.subplots(figsize=(8, 5.6))
    refs = embedded[embedded["record_type"] == "reference"]
    samples = embedded[embedded["record_type"] == "sample_reads"]
    for family, group in refs.groupby("family"):
        ax.scatter(
            group["pc1"],
            group["pc2"],
            s=48,
            label=f"{family} ref",
            color=colors.get(family, "#777777"),
            alpha=0.65,
        )
    for family, group in samples.groupby("family"):
        ax.scatter(
            group["pc1"],
            group["pc2"],
            s=88,
            marker="x",
            linewidths=1.8,
            label=f"{family} sample",
            color=colors.get(family, "#777777"),
            alpha=0.95,
        )
    for _, row in refs.iterrows():
        if row["segment"] == "genome":
            label = row["virus"].replace("Human ", "").replace(" coronavirus", " CoV")
            ax.annotate(label, (row["pc1"], row["pc2"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
    for _, row in samples.iterrows():
        if row["confidence"] in {"high", "medium"}:
            label = f"{row['sample_id'].split('-')[0]} {row['virus'].replace('Human coronavirus ', 'CoV ')}"
            ax.annotate(label, (row["pc1"], row["pc2"]), fontsize=6, xytext=(3, -7), textcoords="offset points")
    ax.set_xlabel("k-mer PC1")
    ax.set_ylabel("k-mer PC2")
    ax.set_title("Reference and sample-read 5-mer embedding")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=6, frameon=False, ncols=2)
    fig.tight_layout()
    out = FIG_DIR / "part3_sample_reference_kmer_clusters.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def main() -> None:
    ensure_dirs()
    meta = pd.read_csv(REF_META)
    records = list(SeqIO.parse(REF_FASTA, "fasta"))
    _, index = kmer_index()
    ref_matrix = kmer_matrix(records, index)

    ref_rows = meta.copy()
    ref_rows["record_id"] = ref_rows["target_id"]
    ref_rows["record_type"] = "reference"
    ref_rows["sample_id"] = ""
    ref_rows["family"] = ref_rows["virus"].map(FAMILY_BY_VIRUS).fillna("Unknown")
    ref_rows["reads_assigned"] = 0
    ref_rows["reads_found"] = 0
    ref_rows["median_identity"] = np.nan
    ref_rows["genome_coverage_breadth"] = np.nan
    ref_rows["confidence"] = ""

    sample_rows, sample_matrix = sample_read_profiles(index)
    matrix = np.vstack([ref_matrix, sample_matrix]) if len(sample_rows) else ref_matrix
    rows = pd.concat([ref_rows, sample_rows], ignore_index=True, sort=False)

    coords = svd_embedding(matrix, n_components=6)
    labels = kmeans(coords[:, :6], N_CLUSTERS)

    embedded = rows.copy()
    for i in range(coords.shape[1]):
        embedded[f"pc{i + 1}"] = coords[:, i]
    embedded["cluster"] = labels
    embedded.to_csv(RESULTS_DIR / "part3_combined_kmer_embeddings.csv", index=False)
    embedded[embedded["record_type"] == "reference"].to_csv(RESULTS_DIR / "part3_kmer_reference_embeddings.csv", index=False)

    summary = cluster_summary(embedded)
    summary.to_csv(RESULTS_DIR / "part3_cluster_summary.csv", index=False)
    comparison = sample_reference_comparison(embedded)
    comparison.to_csv(RESULTS_DIR / "part3_sample_reference_comparison.csv", index=False)
    figure = plot_embedding(embedded)

    weighted_purity = (summary["family_purity"] * summary["n_records"]).sum() / summary["n_records"].sum()
    family_match_rate = comparison["nearest_reference_matches_family"].mean() if not comparison.empty else float("nan")
    virus_match_rate = comparison["nearest_reference_matches_virus"].mean() if not comparison.empty else float("nan")
    print(f"Reference records embedded: {(embedded['record_type'] == 'reference').sum()}")
    print(f"Sample-virus read sets embedded: {(embedded['record_type'] == 'sample_reads').sum()}")
    print(f"Clusters: {N_CLUSTERS}")
    print(f"Weighted family purity: {weighted_purity:.3f}")
    print(f"Nearest-reference family match rate for sample read sets: {family_match_rate:.3f}")
    print(f"Nearest-reference virus match rate for sample read sets: {virus_match_rate:.3f}")
    print(f"Wrote: {RESULTS_DIR / 'part3_combined_kmer_embeddings.csv'}")
    print(f"Wrote: {RESULTS_DIR / 'part3_kmer_reference_embeddings.csv'}")
    print(f"Wrote: {RESULTS_DIR / 'part3_cluster_summary.csv'}")
    print(f"Wrote: {RESULTS_DIR / 'part3_sample_reference_comparison.csv'}")
    print(f"Wrote: {figure}")


if __name__ == "__main__":
    main()
