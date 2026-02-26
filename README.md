# bch709_vibe_coding

Manual for analyzing genomic feature distribution across yeast chromosomes using:
- `data/saccharomyces_cerevisiae.gff.gz`
- `data/chrom.sizes`

## What this project does

The script computes per-chromosome:
- Feature counts: `genes`, `exons`, `tRNAs`, `snoRNAs`
- Feature densities: count per Mb
- Correlation between chromosome size and feature density (Pearson + Spearman)

It also generates summary plots.

## Important feature definitions used

These are strict rules in `scripts/analyze_yeast_features.py`:
- `genes` = GFF type `gene`
- `exons` = GFF type `exon` only
- `tRNAs` = GFF type `tRNA_gene`
- `snoRNAs` = GFF type `snoRNA_gene`

Notes:
- `noncoding_exon` is **not** counted as `exon`.
- `chrmt` and `chrM` are treated as independent chromosome names.
- If a chromosome appears in GFF but has no size in `chrom.sizes`, densities are set to `NaN` and excluded from correlation calculations.

## File structure

- `environment.yml`: Conda environment
- `scripts/analyze_yeast_features.py`: Main analysis script
- `data/`: Input files
- `results/`: Output tables and figures

## Setup (Conda)

From the project root:

```bash
cd /home/wyim/bin/bch709_vibe_coding
conda env create -f environment.yml
conda activate bch709_vibe_coding
```

If the environment already exists and you changed dependencies:

```bash
conda env update -f environment.yml --prune
conda activate bch709_vibe_coding
```

## Run the analysis

```bash
python scripts/analyze_yeast_features.py \
  --gff data/saccharomyces_cerevisiae.gff.gz \
  --chrom-sizes data/chrom.sizes \
  --outdir results
```

## Outputs

The script writes:

- `results/feature_counts_and_density_by_chromosome.csv`
- `results/feature_density_by_chromosome.csv`
- `results/feature_density_correlations.csv`
- `results/feature_counts_by_chromosome.png`
- `results/feature_density_vs_size.png`

## How to read results

- `feature_counts_and_density_by_chromosome.csv`
  - One row per chromosome
  - Contains `size_bp`, raw counts, and density columns (`*_per_mb`)
- `feature_density_correlations.csv`
  - Includes `all_chromosomes` and `nuclear_only` subsets
  - For `nuclear_only`, both `chrM` and `chrmt` are excluded
  - If a feature has constant density (for example all zero), correlation values are `NaN`

## Troubleshooting

- `ModuleNotFoundError` (e.g., matplotlib):
  - Activate env first: `conda activate bch709_vibe_coding`
- Conda cannot download packages:
  - Check network/DNS access to conda channels (`conda-forge`)
- Empty/zero `exons`:
  - Your GFF may not have true `exon` entries (only `noncoding_exon`)
