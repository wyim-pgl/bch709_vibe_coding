#!/usr/bin/env python3
"""Analyze genomic feature distribution across yeast chromosomes."""

from __future__ import annotations

import argparse
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


FEATURE_TYPE_TO_GROUP = {
    "gene": "genes",
    "exon": "exons",
    "tRNA_gene": "tRNAs",
    "snoRNA_gene": "snoRNAs",
}

FEATURE_GROUPS = ["genes", "exons", "tRNAs", "snoRNAs"]


def read_chrom_sizes(path: Path) -> pd.DataFrame:
    chrom_df = pd.read_csv(path, sep="\t", header=None, names=["chromosome", "size_bp"])
    chrom_df = chrom_df.drop_duplicates(subset=["chromosome"], keep="first")
    return chrom_df


def count_features_from_gff(gff_path: Path) -> pd.DataFrame:
    counts = defaultdict(lambda: defaultdict(int))
    with gzip.open(gff_path, "rt") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            feature_type = parts[2]
            group = FEATURE_TYPE_TO_GROUP.get(feature_type)
            if group is not None:
                counts[chrom][group] += 1

    count_df = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(int)
    for group in FEATURE_GROUPS:
        if group not in count_df.columns:
            count_df[group] = 0
    count_df = count_df[FEATURE_GROUPS]
    count_df.index.name = "chromosome"
    return count_df


def make_density_table(chrom_df: pd.DataFrame, count_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = chrom_df.merge(count_df.reset_index(), on="chromosome", how="outer")
    for group in FEATURE_GROUPS:
        merged[group] = merged[group].fillna(0)
        merged[group] = merged[group].astype(int)
        merged[f"{group}_per_mb"] = np.where(
            merged["size_bp"].notna(),
            merged[group] / (merged["size_bp"] / 1_000_000.0),
            np.nan,
        )
    return merged, merged[["chromosome", "size_bp"] + [f"{g}_per_mb" for g in FEATURE_GROUPS]]


def correlation_table(density_df: pd.DataFrame, include_chrM: bool) -> pd.DataFrame:
    work = density_df.copy()
    subset = "all_chromosomes"
    if not include_chrM:
        work = work[~work["chromosome"].isin(["chrM", "chrmt"])]
        subset = "nuclear_only"
    work = work.dropna(subset=["size_bp"] + [f"{g}_per_mb" for g in FEATURE_GROUPS])

    rows = []
    for group in FEATURE_GROUPS:
        y = work[f"{group}_per_mb"].to_numpy(dtype=float)
        x = work["size_bp"].to_numpy(dtype=float)
        if np.isclose(np.nanstd(y), 0.0):
            pearson_r, pearson_p = np.nan, np.nan
            spearman_rho, spearman_p = np.nan, np.nan
        else:
            pearson_r, pearson_p = pearsonr(x, y)
            spearman_rho, spearman_p = spearmanr(x, y)
        rows.append(
            {
                "subset": subset,
                "feature": group,
                "pearson_r": pearson_r,
                "pearson_pvalue": pearson_p,
                "spearman_rho": spearman_rho,
                "spearman_pvalue": spearman_p,
            }
        )
    return pd.DataFrame(rows)


def plot_feature_counts(merged_df: pd.DataFrame, output_path: Path) -> None:
    melted = merged_df.melt(
        id_vars=["chromosome"],
        value_vars=FEATURE_GROUPS,
        var_name="feature",
        value_name="count",
    )
    plt.figure(figsize=(14, 6))
    sns.barplot(data=melted, x="chromosome", y="count", hue="feature")
    plt.title("Genomic Feature Counts by Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("Count")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_density_vs_size(
    density_df: pd.DataFrame,
    corr_nuclear: pd.DataFrame,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
    nuclear = density_df[~density_df["chromosome"].isin(["chrM", "chrmt"])].dropna(subset=["size_bp"])
    mitochondrial = density_df[density_df["chromosome"].isin(["chrM", "chrmt"])].dropna(subset=["size_bp"])

    for ax, group in zip(axes.flat, FEATURE_GROUPS):
        y_col = f"{group}_per_mb"
        sns.regplot(
            data=nuclear,
            x="size_bp",
            y=y_col,
            ax=ax,
            scatter_kws={"s": 45, "alpha": 0.85},
            line_kws={"color": "black", "linewidth": 1.5},
            ci=None,
        )
        if not mitochondrial.empty:
            ax.scatter(
                mitochondrial["size_bp"],
                mitochondrial[y_col],
                s=55,
                color="tomato",
                marker="D",
                label="mito (size-known)",
                zorder=5,
            )
            ax.legend(loc="best", frameon=False)
        r = corr_nuclear.loc[corr_nuclear["feature"] == group, "pearson_r"].iloc[0]
        p = corr_nuclear.loc[corr_nuclear["feature"] == group, "pearson_pvalue"].iloc[0]
        ax.set_title(f"{group}: r={r:.2f}, p={p:.3g}")
        ax.set_xlabel("Chromosome size (bp)")
        ax.set_ylabel("Density (count per Mb)")
        ax.grid(alpha=0.2)

    fig.suptitle("Feature Density vs Chromosome Size (Regression on Nuclear Chromosomes)", y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def print_summary(merged_df: pd.DataFrame, corr_df: pd.DataFrame) -> None:
    missing_size = merged_df[merged_df["size_bp"].isna()]["chromosome"].tolist()
    if missing_size:
        print("Chromosomes with feature counts but missing size (density not computed):")
        print("  " + ", ".join(missing_size))

    print("\nTop density chromosome for each feature (count per Mb):")
    for group in FEATURE_GROUPS:
        col = f"{group}_per_mb"
        valid = merged_df.dropna(subset=[col])
        if valid.empty or np.isclose(valid[col].max(), 0.0):
            print(f"  {group:7s}: no non-zero density found")
            continue
        row = valid.loc[valid[col].idxmax(), ["chromosome", col]]
        print(f"  {group:7s}: {row['chromosome']} ({row[col]:.2f}/Mb)")

    print("\nCorrelation summary (nuclear chromosomes only):")
    nuclear = corr_df[corr_df["subset"] == "nuclear_only"]
    for _, row in nuclear.iterrows():
        print(
            f"  {row['feature']:7s}: Pearson r={row['pearson_r']:.3f} (p={row['pearson_pvalue']:.3g}), "
            f"Spearman rho={row['spearman_rho']:.3f} (p={row['spearman_pvalue']:.3g})"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze genomic feature distributions and density-size correlations."
    )
    parser.add_argument("--gff", type=Path, required=True, help="Path to .gff.gz file")
    parser.add_argument("--chrom-sizes", type=Path, required=True, help="Path to chrom.sizes file")
    parser.add_argument("--outdir", type=Path, default=Path("results"), help="Output directory")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    chrom_df = read_chrom_sizes(args.chrom_sizes)
    count_df = count_features_from_gff(args.gff)
    merged_df, density_df = make_density_table(chrom_df, count_df)

    corr_all = correlation_table(density_df, include_chrM=True)
    corr_nuclear = correlation_table(density_df, include_chrM=False)
    corr_df = pd.concat([corr_all, corr_nuclear], ignore_index=True)

    merged_df.to_csv(args.outdir / "feature_counts_and_density_by_chromosome.csv", index=False)
    density_df.to_csv(args.outdir / "feature_density_by_chromosome.csv", index=False)
    corr_df.to_csv(args.outdir / "feature_density_correlations.csv", index=False)

    sns.set_theme(style="whitegrid")
    plot_feature_counts(merged_df, args.outdir / "feature_counts_by_chromosome.png")
    plot_density_vs_size(density_df, corr_nuclear, args.outdir / "feature_density_vs_size.png")

    print_summary(merged_df, corr_df)
    print(f"\nWrote results to: {args.outdir.resolve()}")


if __name__ == "__main__":
    main()
