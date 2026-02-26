# bch709_vibe_coding

This project analyzes genomic feature distribution across yeast chromosomes using:
- `data/saccharomyces_cerevisiae.gff.gz`
- `data/chrom.sizes`

## 1) Installation and Environment Setup (First)

### Prerequisite
- Conda is installed and available in your shell (`conda --version`).

### Recommended: create environment from `environment.yml`

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

### Verify installation (recommended)

```bash
python -c "import pandas, numpy, scipy, matplotlib, seaborn; print('Python OK')"
```

### Optional reference command style from BCH709 Step 0C
If you prefer explicit package install commands (instead of YAML), BCH709 Step 0C uses this style:

```bash
conda create -n bch709_vibe_coding -y -c conda-forge \
  python=3.11 pandas numpy scipy matplotlib seaborn
conda activate bch709_vibe_coding
python -c "import pandas, numpy, scipy, matplotlib, seaborn; print('Python OK')"
```

### Keep environment reproducible (recommended)
```bash
conda env export -n bch709_vibe_coding > environment.yml
```

### If installation fails
- Confirm active env: `conda activate bch709_vibe_coding`
- Retry install/update command
- Check network access to `conda-forge`

## 2) How to Run

```bash
python scripts/analyze_yeast_features.py \
  --gff data/saccharomyces_cerevisiae.gff.gz \
  --chrom-sizes data/chrom.sizes \
  --outdir results
```

## 3) Output Files

- `results/feature_counts_and_density_by_chromosome.csv`
- `results/feature_density_by_chromosome.csv`
- `results/feature_density_correlations.csv`
- `results/feature_counts_by_chromosome.png`
- `results/feature_density_vs_size.png`

## 4) What the Script Does

Per chromosome:
- Counts `genes`, `exons`, `tRNAs`, `snoRNAs`
- Density (count per Mb)
- Correlation of feature density vs chromosome size (Pearson and Spearman)

## 5) Feature Definition Rules in This Project

- `genes` = GFF type `gene`
- `exons` = GFF type `exon` only
- `tRNAs` = GFF type `tRNA_gene`
- `snoRNAs` = GFF type `snoRNA_gene`
- `noncoding_exon` is not counted as `exons`
- `chrmt` and `chrM` are treated as independent chromosome names
- If a chromosome has counts but no size in `chrom.sizes`, density is `NaN` and excluded from correlation

## 6) Data Check (if needed)
If your data files are missing, download with:

```bash
mkdir -p data
curl -L -o data/saccharomyces_cerevisiae.gff.gz \
  http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz
curl -L -o data/chrom.sizes \
  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.chrom.sizes
```

## 7) Reference
Environment setup details above were aligned with BCH709 Vibe Coding Step 0C:
- https://bch709.plantgenomicslab.org/vibe_coding/index.html#step-0c-environment-creation-commands
