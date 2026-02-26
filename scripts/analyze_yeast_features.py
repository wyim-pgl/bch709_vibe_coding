#!/usr/bin/env python3
"""Analyze genomic feature distribution across yeast chromosomes.

This script reads a GFF3 annotation file and a chromosome-sizes file for
Saccharomyces cerevisiae, counts genomic features (genes, CDS, tRNAs,
snoRNAs) per chromosome, computes feature densities (per megabase), and
tests for correlations between feature density and chromosome size.

Outputs:
    - CSV tables of counts, densities, and correlation statistics
    - A grouped bar chart of feature counts by chromosome
    - A 2x2 scatter-plot panel of density vs. chromosome size

Usage:
    python scripts/analyze_yeast_features.py \
        --gff data/saccharomyces_cerevisiae.gff.gz \
        --chrom-sizes data/chrom.sizes \
        --outdir results
"""

# Enable postponed evaluation of annotations (PEP 604 union syntax, etc.)
from __future__ import annotations

# Standard-library imports
import argparse  # command-line argument parsing
import gzip  # reading gzip-compressed GFF files
import logging  # structured log output instead of bare print()
import re  # regular expressions for Roman-numeral chromosome sorting
import sys  # for sys.exit on validation failures
from collections import defaultdict  # auto-initializing nested dicts
from pathlib import Path  # object-oriented filesystem paths

# Third-party imports
import matplotlib.pyplot as plt  # low-level plotting API
import numpy as np  # numerical array operations
import pandas as pd  # tabular data manipulation
import seaborn as sns  # high-level statistical plotting
from scipy.stats import pearsonr, spearmanr  # correlation coefficients

# ---------------------------------------------------------------------------
# Configure module-level logger so all messages go to stderr
# ---------------------------------------------------------------------------
logger = logging.getLogger(__name__)  # create logger scoped to this module

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Map raw GFF3 feature-type strings to human-friendly group names.
# NOTE: SGD GFF3 does not emit bare "exon" features; it uses "CDS" for
# protein-coding segments and "noncoding_exon" for ncRNA exons.  We count
# "CDS" instead of "exon" so the column is not all-zeros.
FEATURE_TYPE_TO_GROUP: dict[str, str] = {
    "gene": "genes",  # protein-coding and non-coding gene loci
    "CDS": "CDS",  # coding sequences (replaces the absent "exon" type)
    "tRNA_gene": "tRNAs",  # transfer-RNA gene loci
    "snoRNA_gene": "snoRNAs",  # small nucleolar RNA gene loci
}

# Ordered list of the group names we report on (matches dict values above)
FEATURE_GROUPS: list[str] = ["genes", "CDS", "tRNAs", "snoRNAs"]

# Canonical name we normalise the mitochondrial chromosome to.
# The GFF uses "chrmt" while chrom.sizes uses "chrM"; we unify to "chrM".
MITO_CANONICAL: str = "chrM"

# Set of alternative mitochondrial chromosome names found across data sources
MITO_ALIASES: set[str] = {"chrM", "chrmt", "chrMito", "MT"}

# Minimum number of observations required to compute a correlation
MIN_CORR_SAMPLES: int = 4


# ---------------------------------------------------------------------------
# Roman-numeral helper for biologically meaningful chromosome ordering
# ---------------------------------------------------------------------------

# Mapping of individual Roman-numeral characters to integer values
_ROMAN_VALUES: dict[str, int] = {
    "I": 1, "V": 5, "X": 10, "L": 50, "C": 100, "D": 500, "M": 1000,
}


def _roman_to_int(roman: str) -> int:
    """Convert a Roman-numeral string (e.g. 'XVI') to an integer (16).

    Parameters
    ----------
    roman : str
        Upper-case Roman numeral.

    Returns
    -------
    int
        Decimal equivalent.
    """
    total: int = 0  # accumulator for the converted value
    prev: int = 0  # value of the previous character (for subtractive notation)
    # Iterate from right to left so subtractive pairs (IV, IX, …) are handled
    for char in reversed(roman):
        value: int = _ROMAN_VALUES.get(char, 0)  # look up numeric value
        if value < prev:
            # Subtractive case: e.g. 'I' before 'V' means subtract 1
            total -= value
        else:
            # Normal additive case
            total += value
        prev = value  # remember this value for the next iteration
    return total


def chrom_sort_key(name: str) -> tuple[int, int, str]:
    """Return a sort key that orders yeast chromosomes biologically.

    Ordering logic:
        1. Nuclear chromosomes by Roman-numeral number (chrI … chrXVI)
        2. Mitochondrial chromosome (chrM) last among standard chroms
        3. Any unexpected names sorted alphabetically at the very end

    Parameters
    ----------
    name : str
        Chromosome name, e.g. ``"chrIV"`` or ``"chrM"``.

    Returns
    -------
    tuple[int, int, str]
        A three-element sort key: (priority-group, numeric-order, name).
    """
    # Check whether this is the mitochondrial chromosome
    if name in MITO_ALIASES or name == MITO_CANONICAL:
        # Place mito after the 16 nuclear chromosomes (group 1)
        return (1, 0, name)

    # Try to extract a Roman-numeral suffix after "chr"
    match = re.match(r"^chr([IVXLCDM]+)$", name, re.IGNORECASE)
    if match:
        # Nuclear chromosome — group 0, sorted by numeric value
        return (0, _roman_to_int(match.group(1).upper()), name)

    # Fallback for any non-standard chromosome names
    return (2, 0, name)


# ---------------------------------------------------------------------------
# Normalise mitochondrial chromosome names
# ---------------------------------------------------------------------------

def _normalise_mito(name: str) -> str:
    """Replace any mitochondrial alias with the canonical name.

    Parameters
    ----------
    name : str
        Raw chromosome name from an input file.

    Returns
    -------
    str
        ``MITO_CANONICAL`` if *name* is a known alias, otherwise *name* unchanged.
    """
    if name in MITO_ALIASES:
        return MITO_CANONICAL  # unify all aliases to one canonical name
    return name


# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------

def read_chrom_sizes(path: Path) -> pd.DataFrame:
    """Read a two-column TAB-delimited chromosome-sizes file.

    Expected format (no header):
        chrI\\t230218
        chrII\\t813184
        ...

    Duplicate chromosome names are dropped (first occurrence kept).

    Parameters
    ----------
    path : Path
        Filesystem path to the chrom.sizes file.

    Returns
    -------
    pd.DataFrame
        Columns: ``["chromosome", "size_bp"]``.
    """
    # Read the raw two-column TSV into a DataFrame
    chrom_df: pd.DataFrame = pd.read_csv(
        path,  # path to the input file
        sep="\t",  # columns separated by TAB
        header=None,  # file has no header row
        names=["chromosome", "size_bp"],  # assign our own column names
    )

    # Normalise mitochondrial chromosome name so it matches the GFF data
    chrom_df["chromosome"] = chrom_df["chromosome"].apply(_normalise_mito)

    # Drop any duplicate chromosome entries, keeping the first occurrence
    chrom_df = chrom_df.drop_duplicates(subset=["chromosome"], keep="first")

    # Log how many chromosomes were loaded for transparency
    logger.info("Loaded %d chromosome sizes from %s", len(chrom_df), path)

    return chrom_df


def count_features_from_gff(gff_path: Path) -> pd.DataFrame:
    """Parse a gzip-compressed GFF3 file and count features per chromosome.

    Only feature types listed in ``FEATURE_TYPE_TO_GROUP`` are tallied.
    Parsing stops at the ``##FASTA`` directive (sequence section).

    Parameters
    ----------
    gff_path : Path
        Filesystem path to the ``.gff.gz`` file.

    Returns
    -------
    pd.DataFrame
        Index: ``chromosome``; columns: one per group in ``FEATURE_GROUPS``,
        values are integer counts.
    """
    # Nested defaultdict: outer key = chromosome, inner key = group name
    counts: defaultdict[str, defaultdict[str, int]] = defaultdict(
        lambda: defaultdict(int)
    )

    # Open the gzip-compressed file in text mode
    with gzip.open(gff_path, "rt") as handle:
        for line in handle:
            # The ##FASTA directive marks the start of embedded sequences; stop
            if line.startswith("##FASTA"):
                break

            # Skip blank lines and comment lines (start with '#')
            if not line or line.startswith("#"):
                continue

            # GFF3 lines are TAB-separated with at least 9 columns
            parts: list[str] = line.rstrip("\n").split("\t")

            # Guard against malformed lines with fewer than 3 fields
            if len(parts) < 3:
                continue

            # Column 0 = chromosome, column 2 = feature type
            chrom: str = _normalise_mito(parts[0])  # unify mito names
            feature_type: str = parts[2]

            # Look up whether this feature type is one we care about
            group: str | None = FEATURE_TYPE_TO_GROUP.get(feature_type)

            if group is not None:
                # Increment the count for this chromosome + group pair
                counts[chrom][group] += 1

    # Convert the nested dict into a DataFrame (chromosomes as rows)
    count_df: pd.DataFrame = (
        pd.DataFrame.from_dict(counts, orient="index")  # rows = chroms
        .fillna(0)  # missing groups become zero
        .astype(int)  # ensure integer counts
    )

    # Guarantee every expected group column exists (even if zero everywhere)
    for group in FEATURE_GROUPS:
        if group not in count_df.columns:
            count_df[group] = 0  # add a zero-filled column

    # Reorder columns to the canonical group order
    count_df = count_df[FEATURE_GROUPS]

    # Name the index so it merges cleanly later
    count_df.index.name = "chromosome"

    # Log summary counts for verification
    logger.info(
        "Counted features across %d chromosomes from %s", len(count_df), gff_path,
    )

    return count_df


# ---------------------------------------------------------------------------
# Density calculation
# ---------------------------------------------------------------------------

def make_density_table(
    chrom_df: pd.DataFrame,
    count_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge chromosome sizes with feature counts and compute densities.

    Density is defined as ``count / (size_bp / 1,000,000)``, i.e. features
    per megabase.  Chromosomes present in only one of the two inputs will
    appear with ``NaN`` for the missing values and a warning is logged.

    Parameters
    ----------
    chrom_df : pd.DataFrame
        Output of :func:`read_chrom_sizes`.
    count_df : pd.DataFrame
        Output of :func:`count_features_from_gff`.

    Returns
    -------
    merged : pd.DataFrame
        Full table with raw counts, sizes, and density columns.
    density_only : pd.DataFrame
        Subset with ``chromosome``, ``size_bp``, and ``*_per_mb`` columns.
    """
    # Outer-merge keeps chromosomes from both sides; mismatches become NaN
    merged: pd.DataFrame = chrom_df.merge(
        count_df.reset_index(),  # move chromosome from index to column
        on="chromosome",  # join key
        how="outer",  # keep all chromosomes from both sides
    )

    # --- Warn the user about any mismatched chromosomes ---
    # Chromosomes in the GFF but not in chrom.sizes (will have NaN size)
    missing_size: list[str] = merged.loc[
        merged["size_bp"].isna(), "chromosome"
    ].tolist()
    if missing_size:
        logger.warning(
            "Chromosomes in GFF but NOT in chrom.sizes (density unavailable): %s",
            ", ".join(missing_size),
        )

    # Chromosomes in chrom.sizes but not in the GFF (will have NaN counts)
    for group in FEATURE_GROUPS:
        merged[group] = merged[group].fillna(0).astype(int)  # NaN → 0

    missing_feats: list[str] = merged.loc[
        merged[FEATURE_GROUPS].sum(axis=1) == 0, "chromosome"
    ].tolist()
    if missing_feats:
        logger.info(
            "Chromosomes with zero features across all groups: %s",
            ", ".join(missing_feats),
        )

    # Compute density columns: features per megabase
    for group in FEATURE_GROUPS:
        merged[f"{group}_per_mb"] = np.where(
            merged["size_bp"].notna(),  # only compute where size is known
            merged[group] / (merged["size_bp"] / 1_000_000.0),  # count / Mb
            np.nan,  # no density without a known chromosome size
        )

    # Sort rows in biological chromosome order for readable output
    merged = merged.sort_values(
        "chromosome",
        key=lambda col: col.map(chrom_sort_key),  # apply our custom sort
    ).reset_index(drop=True)

    # Build the density-only subset (chromosome + size + per-Mb columns)
    density_cols: list[str] = ["chromosome", "size_bp"] + [
        f"{g}_per_mb" for g in FEATURE_GROUPS
    ]
    density_only: pd.DataFrame = merged[density_cols]

    return merged, density_only


# ---------------------------------------------------------------------------
# Correlation analysis
# ---------------------------------------------------------------------------

def correlation_table(
    density_df: pd.DataFrame,
    include_chrM: bool,
) -> pd.DataFrame:
    """Compute Pearson and Spearman correlations between density and size.

    Parameters
    ----------
    density_df : pd.DataFrame
        Output density table from :func:`make_density_table`.
    include_chrM : bool
        If ``False``, exclude mitochondrial chromosomes before computing.

    Returns
    -------
    pd.DataFrame
        One row per feature group with r, rho, and p-values.
    """
    # Work on a copy to avoid mutating the caller's data
    work: pd.DataFrame = density_df.copy()

    # Label the subset for output identification
    subset: str = "all_chromosomes"

    if not include_chrM:
        # Exclude mitochondrial chromosome(s) from the analysis
        work = work[work["chromosome"] != MITO_CANONICAL]
        subset = "nuclear_only"

    # Drop rows where size or any density is NaN (incomplete data)
    required_cols: list[str] = ["size_bp"] + [
        f"{g}_per_mb" for g in FEATURE_GROUPS
    ]
    work = work.dropna(subset=required_cols)

    # Collect one result row per feature group
    rows: list[dict] = []

    for group in FEATURE_GROUPS:
        # Extract the density and size arrays as float64
        y: np.ndarray = work[f"{group}_per_mb"].to_numpy(dtype=float)
        x: np.ndarray = work["size_bp"].to_numpy(dtype=float)

        # Guard: need at least MIN_CORR_SAMPLES data points
        if len(x) < MIN_CORR_SAMPLES:
            logger.warning(
                "Only %d samples for %s (%s); skipping correlation",
                len(x), group, subset,
            )
            pearson_r = pearson_p = spearman_rho = spearman_p = np.nan

        # Guard: if either variable has zero variance, correlation is undefined
        elif np.isclose(np.nanstd(x), 0.0) or np.isclose(np.nanstd(y), 0.0):
            logger.warning(
                "Zero variance in %s (%s); correlation undefined", group, subset,
            )
            pearson_r = pearson_p = spearman_rho = spearman_p = np.nan

        else:
            # Compute Pearson (linear) and Spearman (rank) correlations
            pearson_r, pearson_p = pearsonr(x, y)
            spearman_rho, spearman_p = spearmanr(x, y)

        # Append the result row for this feature group
        rows.append(
            {
                "subset": subset,
                "feature": group,
                "n_chromosomes": len(x),
                "pearson_r": pearson_r,
                "pearson_pvalue": pearson_p,
                "spearman_rho": spearman_rho,
                "spearman_pvalue": spearman_p,
            }
        )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_feature_counts(
    merged_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Create a grouped bar chart of feature counts by chromosome.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Full merged table with raw counts.
    output_path : Path
        Where to save the PNG image.
    """
    # Melt wide-format counts into long format for seaborn
    melted: pd.DataFrame = merged_df.melt(
        id_vars=["chromosome"],  # keep chromosome as identifier
        value_vars=FEATURE_GROUPS,  # one row per feature group
        var_name="feature",  # new column holding the group name
        value_name="count",  # new column holding the count value
    )

    # Create a figure with generous width for many chromosome labels
    plt.figure(figsize=(14, 6))

    # Draw grouped bars; chromosome order comes from the pre-sorted DataFrame
    sns.barplot(
        data=melted,
        x="chromosome",  # x-axis: chromosome names
        y="count",  # y-axis: feature counts
        hue="feature",  # colour-code by feature group
        order=merged_df["chromosome"].tolist(),  # preserve biological sort
    )

    # Add informative title and axis labels
    plt.title("Genomic Feature Counts by Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("Count")

    # Rotate x-tick labels so long names don't overlap
    plt.xticks(rotation=45, ha="right")

    # Compact the layout to prevent label clipping
    plt.tight_layout()

    # Save figure at 200 DPI for reasonable file size + readability
    plt.savefig(output_path, dpi=200)

    # Close the figure to free memory
    plt.close()

    # Log confirmation
    logger.info("Saved feature-count bar chart to %s", output_path)


def plot_density_vs_size(
    density_df: pd.DataFrame,
    corr_nuclear: pd.DataFrame,
    output_path: Path,
) -> None:
    """Create a 2x2 panel of density-vs-size scatter plots with regression.

    Nuclear chromosomes are shown as dots with a regression line; the
    mitochondrial chromosome (if present and sized) is plotted as a red
    diamond for visual distinction.

    Parameters
    ----------
    density_df : pd.DataFrame
        Density table from :func:`make_density_table`.
    corr_nuclear : pd.DataFrame
        Correlation stats (nuclear-only) from :func:`correlation_table`.
    output_path : Path
        Where to save the PNG image.
    """
    # Create a 2-row × 2-column grid of axes sharing the x-axis
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)

    # Separate nuclear and mitochondrial chromosomes
    nuclear: pd.DataFrame = density_df[
        density_df["chromosome"] != MITO_CANONICAL
    ].dropna(subset=["size_bp"])  # need a valid size to plot

    mitochondrial: pd.DataFrame = density_df[
        density_df["chromosome"] == MITO_CANONICAL
    ].dropna(subset=["size_bp"])  # only plot if size is known

    # Iterate over the four feature groups, one per subplot
    for ax, group in zip(axes.flat, FEATURE_GROUPS):
        # Column name for this group's density
        y_col: str = f"{group}_per_mb"

        # Draw the regression on nuclear chromosomes only
        sns.regplot(
            data=nuclear,
            x="size_bp",  # x-axis: chromosome size in base pairs
            y=y_col,  # y-axis: density (features per Mb)
            ax=ax,  # target this specific subplot
            scatter_kws={"s": 45, "alpha": 0.85},  # dot size and opacity
            line_kws={"color": "black", "linewidth": 1.5},  # regression line
            ci=None,  # suppress confidence interval band for clarity
        )

        # Overlay mitochondrial chromosome as a distinct marker
        if not mitochondrial.empty:
            ax.scatter(
                mitochondrial["size_bp"],  # mito chromosome size
                mitochondrial[y_col],  # mito density value
                s=55,  # slightly larger dot
                color="tomato",  # red to stand out
                marker="D",  # diamond shape
                label="mitochondrial",  # legend text
                zorder=5,  # draw on top of regression elements
            )
            ax.legend(loc="best", frameon=False)  # show legend without frame

        # Annotate subplot with Pearson r and p-value from the correlation table
        corr_row = corr_nuclear.loc[corr_nuclear["feature"] == group]
        if not corr_row.empty:
            r_val: float = corr_row["pearson_r"].iloc[0]
            p_val: float = corr_row["pearson_pvalue"].iloc[0]
            ax.set_title(f"{group}: r={r_val:.2f}, p={p_val:.3g}")
        else:
            # Fallback if correlation data is missing for this group
            ax.set_title(f"{group}: correlation N/A")

        # Label axes
        ax.set_xlabel("Chromosome size (bp)")
        ax.set_ylabel("Density (count per Mb)")

        # Light grid for readability
        ax.grid(alpha=0.2)

    # Overall figure title (slightly above the top row)
    fig.suptitle(
        "Feature Density vs Chromosome Size (Regression on Nuclear Chromosomes)",
        y=1.02,
    )

    # Adjust subplot spacing so titles and labels don't overlap
    fig.tight_layout()

    # Save with tight bounding box to include the suptitle
    fig.savefig(output_path, dpi=220, bbox_inches="tight")

    # Free memory
    plt.close(fig)

    # Log confirmation
    logger.info("Saved density-vs-size scatter plots to %s", output_path)


# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------

def print_summary(
    merged_df: pd.DataFrame,
    corr_df: pd.DataFrame,
) -> None:
    """Print a human-readable summary of the analysis to stdout.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Full merged table (counts + densities).
    corr_df : pd.DataFrame
        Combined correlation results (all + nuclear subsets).
    """
    # --- Report chromosomes with missing size data ---
    missing_size: list[str] = merged_df[
        merged_df["size_bp"].isna()
    ]["chromosome"].tolist()

    if missing_size:
        print("Chromosomes with feature counts but missing size (density not computed):")
        print("  " + ", ".join(missing_size))

    # --- Top-density chromosome for each feature group ---
    print("\nTop density chromosome for each feature (count per Mb):")

    for group in FEATURE_GROUPS:
        col: str = f"{group}_per_mb"  # density column name

        # Only consider rows where density is valid and non-zero
        valid: pd.DataFrame = merged_df.dropna(subset=[col])

        if valid.empty or np.isclose(valid[col].max(), 0.0):
            # No meaningful density data for this feature group
            print(f"  {group:7s}: no non-zero density found")
            continue

        # Find the row with the highest density
        best = valid.loc[valid[col].idxmax(), ["chromosome", col]]
        print(f"  {group:7s}: {best['chromosome']} ({best[col]:.2f}/Mb)")

    # --- Correlation summary (nuclear chromosomes only) ---
    print("\nCorrelation summary (nuclear chromosomes only):")

    # Filter to the nuclear-only subset
    nuclear: pd.DataFrame = corr_df[corr_df["subset"] == "nuclear_only"]

    for _, row in nuclear.iterrows():
        # Print Pearson and Spearman statistics side by side
        print(
            f"  {row['feature']:7s}: "
            f"Pearson r={row['pearson_r']:.3f} (p={row['pearson_pvalue']:.3g}), "
            f"Spearman rho={row['spearman_rho']:.3f} (p={row['spearman_pvalue']:.3g})"
        )


# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """Parse and validate command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with validated file paths.
    """
    # Create the top-level parser with a description
    parser = argparse.ArgumentParser(
        description=(
            "Analyze genomic feature distributions and density-size "
            "correlations across yeast chromosomes."
        ),
    )

    # Required: path to the gzip-compressed GFF3 annotation file
    parser.add_argument(
        "--gff",
        type=Path,
        required=True,
        help="Path to gzip-compressed GFF3 file (.gff.gz)",
    )

    # Required: path to the two-column chromosome-sizes file
    parser.add_argument(
        "--chrom-sizes",
        type=Path,
        required=True,
        help="Path to TAB-separated chrom.sizes file (chrom<TAB>size)",
    )

    # Optional: output directory (created automatically if absent)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results"),
        help="Output directory for CSVs and plots (default: results/)",
    )

    # Parse the arguments from sys.argv
    args: argparse.Namespace = parser.parse_args()

    # --- Validate that input files exist before proceeding ---
    if not args.gff.exists():
        parser.error(f"GFF file not found: {args.gff}")

    if not args.chrom_sizes.exists():
        parser.error(f"Chromosome-sizes file not found: {args.chrom_sizes}")

    # Warn if the GFF file does not look gzip-compressed
    if not args.gff.name.endswith(".gz"):
        logger.warning(
            "GFF path %s does not end with '.gz'; expecting gzip-compressed input",
            args.gff,
        )

    return args


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Run the full analysis pipeline: read, count, correlate, plot, report."""

    # --- Set up logging to stderr at INFO level ---
    logging.basicConfig(
        level=logging.INFO,  # show INFO and above
        format="%(levelname)s: %(message)s",  # simple format without timestamp
        stream=sys.stderr,  # keep stdout clean for data output
    )

    # --- Parse and validate CLI arguments ---
    args: argparse.Namespace = parse_args()

    # Create the output directory (and any missing parents) if needed
    args.outdir.mkdir(parents=True, exist_ok=True)

    # --- Set the seaborn theme early, before any plotting ---
    sns.set_theme(style="whitegrid")

    # --- Step 1: Read input data ---
    logger.info("Reading chromosome sizes …")
    chrom_df: pd.DataFrame = read_chrom_sizes(args.chrom_sizes)

    logger.info("Counting features from GFF …")
    count_df: pd.DataFrame = count_features_from_gff(args.gff)

    # --- Step 2: Compute densities ---
    logger.info("Computing feature densities …")
    merged_df, density_df = make_density_table(chrom_df, count_df)

    # --- Step 3: Correlation analysis (with and without mitochondria) ---
    logger.info("Computing correlations …")
    corr_all: pd.DataFrame = correlation_table(density_df, include_chrM=True)
    corr_nuclear: pd.DataFrame = correlation_table(density_df, include_chrM=False)

    # Concatenate both subsets into one table
    corr_df: pd.DataFrame = pd.concat(
        [corr_all, corr_nuclear], ignore_index=True,
    )

    # --- Step 4: Write CSV outputs ---
    merged_csv: Path = args.outdir / "feature_counts_and_density_by_chromosome.csv"
    density_csv: Path = args.outdir / "feature_density_by_chromosome.csv"
    corr_csv: Path = args.outdir / "feature_density_correlations.csv"

    merged_df.to_csv(merged_csv, index=False)  # full counts + density table
    density_df.to_csv(density_csv, index=False)  # density-only table
    corr_df.to_csv(corr_csv, index=False)  # correlation statistics

    logger.info("Wrote %s", merged_csv)
    logger.info("Wrote %s", density_csv)
    logger.info("Wrote %s", corr_csv)

    # --- Step 5: Generate plots ---
    counts_png: Path = args.outdir / "feature_counts_by_chromosome.png"
    density_png: Path = args.outdir / "feature_density_vs_size.png"

    plot_feature_counts(merged_df, counts_png)
    plot_density_vs_size(density_df, corr_nuclear, density_png)

    # --- Step 6: Print human-readable summary to stdout ---
    print_summary(merged_df, corr_df)
    print(f"\nWrote results to: {args.outdir.resolve()}")


# Standard Python idiom: only run main() when executed as a script
if __name__ == "__main__":
    main()
