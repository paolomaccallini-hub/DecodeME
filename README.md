# DecodeME Pipeline

This repository contains R scripts for processing and fine-mapping DecodeME GWAS-1 summary statistics. The pipeline requires both R and Python and has been tested with **R 4.4.1** and **Python 3.12**.
## Overview
DecodeME is an R-based workflow for processing and fine-mapping DecodeME GWAS-1 summary statistics. It harmonises inputs, converts genome builds, and retrieves linkage disequilibrium (LD) information for downstream analyses. The pipeline has been tested with **R 4.4.1** and **Python 3.12**.

## Repository structure

- `DecodeME_main.R` – orchestrates filtering, format conversion, and fine mapping.
- `DecodeME_func.R` – helper functions and package setup.
- `DecodeME_config.yml` – configuration file specifying filters, API keys, and phenotype settings.
- `My_genes_DEcodeME.csv` – example output file created by `DecodeME_main.R`.
- `DecodeME_main.R` – orchestrates filtering, liftover, and fine-mapping.
- `DecodeME_func.R` – helper functions and package installation.
- `DecodeME_config.yml` – configuration for filters, API keys, and phenotypes.

## Data sources

- **Summary statistics** – bundled in `Data/DecodeME/DecodeME_summary.zip` and
  derived from DecodeME GWAS-1. The pipeline unpacks the relevant phenotype file
  before processing.
- **LD matrices** – fetched on demand from the UK Biobank LD reference panel
  hosted by the Broad Institute at
  `https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/` and converted to
  RDS files for fine-mapping.
- **Summary statistics** – packaged in `Data/DecodeME/DecodeME_summary.zip`, derived from DecodeME GWAS‑1.
- **LD matrices** – downloaded on demand from the UK Biobank LD reference hosted by the Broad Institute (<https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/>).
- **RegulomeDB scores** – obtain from <https://www.regulomedb.org/regulome-search/> and place the score file in the repository root.

## Requirements
- R (tested with 4.4.1)
- Python (tested with 3.12)
- ≥10 GB free disk space
- Internet access to retrieve LD matrices

- R (tested with version 4.4.1)
- Python (tested with version 3.12)
- At least 10 GB of free disk space
- RegulomeDB score file downloaded from [https://www.regulomedb.org/regulome-search/](https://www.regulomedb.org/regulome-search/) and placed in the repository root
- The internal folder structure is generated automatically by the scripts

## Setup
The internal directory tree is created automatically by `DecodeME_main.R`.

1. Install [R](https://www.r-project.org/) and [Python](https://www.python.org/) along with the required R packages. The scripts make use of Bioconductor resources such as `SNPlocs.Hsapiens.dbSNP155.GRCh38` and `BSgenome.Hsapiens.NCBI.GRCh38`.
## Quick start
1. Install [R](https://www.r-project.org/) and [Python](https://www.python.org/) plus required R packages (e.g., `SNPlocs.Hsapiens.dbSNP155.GRCh38`, `BSgenome.Hsapiens.NCBI.GRCh38`).
2. Download the RegulomeDB score file and place it in the repository root.
3. Adjust `DecodeME_config.yml` to provide paths, API keys, and filtering thresholds.
4. Run the main script (it will automatically build the required folder hierarchy):
3. Configure paths and filters in `DecodeME_config.yml`.
4. Run:
   ```bash
   Rscript DecodeME_main.R
   ```

## Liftover and harmonization

`DecodeME_main.R` standardizes the summary statistics with
[`format_sumstats`](https://github.com/neurogenomics/MungeSumstats), which
harmonizes allele columns and **lifts coordinates from GRCh38 to GRCh37**. The
conversion ensures that downstream steps operate on GRCh37 positions, matching
the genome build used by the UK Biobank LD matrices.
`DecodeME_main.R` invokes [`format_sumstats`](https://github.com/neurogenomics/MungeSumstats) to standardize column names, align alleles, and **convert coordinates from GRCh38 to GRCh37** to match the UK Biobank LD reference.

## INFO column approximation
Because the DecodeME summary statistics lack an INFO field, the workflow estimates a proxy:
\[
\mathrm{INFO_{proxy}} = \frac{1}{\mathrm{SE}^2 \times N_{\mathrm{eff}} \times 2p(1-p)}
\]
where \(p\) is the minor allele frequency and \(N_{\mathrm{eff}} = N \times \pi \times (1-\pi)\) with \(\pi = N_{\text{cases}}/N\). Values are truncated to [0,1].

The original summary statistics do not include an INFO score. The pipeline estimates a proxy INFO value using:

\[ \text{INFO\_proxy} = \frac{1}{\text{SE}^2 \times N_{\text{eff}} \times 2p(1-p)} \]

where \( p \) is the minor allele frequency and \( N_{\text{eff}} = N \times \pi \times (1-\pi) \) with \( \pi = N_{\text{cases}} / N \). Values are truncated to lie within [0, 1]. This approximation is employed because INFO values are not provided in the input data.
## Output
An example output `My_genes_DEcodeME.csv` lists genes prioritized after fine‑mapping and annotation steps.

## License
This project is distributed without an explicit license; contact the authors for reuse permissions.

This project is distributed without an explicit license. Please contact the authors for reuse permissions.
