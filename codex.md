# bch709_vibe_coding Quick Run Guide

## 1. Create and activate Conda environment

```bash
cd /home/wyim/bin/bch709_vibe_coding
conda env create -f environment.yml
conda activate bch709_vibe_coding
```

If environment already exists:

```bash
conda env update -f environment.yml --prune
conda activate bch709_vibe_coding
```

## 2. Run the analysis code

```bash
python scripts/analyze_yeast_features.py \
  --gff data/saccharomyces_cerevisiae.gff.gz \
  --chrom-sizes data/chrom.sizes \
  --outdir results
```

## 3. Output files

After running, check:

- `results/feature_counts_and_density_by_chromosome.csv`
- `results/feature_density_by_chromosome.csv`
- `results/feature_density_correlations.csv`
- `results/feature_counts_by_chromosome.png`
- `results/feature_density_vs_size.png`

## Notes

- `exons` counts only true `exon` features (not `noncoding_exon`).
- `chrmt` and `chrM` are treated as independent chromosome names.
