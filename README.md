# DecodeME Pipeline

This repository contains R scripts for processing and fine-mapping DecodeME GWAS-1 summary statistics.

## Repository structure

- `DecodeME_main.R` – orchestrates filtering, format conversion, and fine mapping.
- `DecodeME_func.R` – helper functions and package setup.
- `DecodeME_config.yml` – configuration file specifying filters, API keys, and phenotype settings.
- `My_genes_DEcodeME.csv` – example output file created by `DecodeME_main.R`.

## Setup

1. Install [R](https://www.r-project.org/) and the required R packages. The scripts make use of Bioconductor resources such as `SNPlocs.Hsapiens.dbSNP155.GRCh38` and `BSgenome.Hsapiens.NCBI.GRCh38`.
2. Adjust `DecodeME_config.yml` to provide paths, API keys, and filtering thresholds.
3. Run the main script:
   ```bash
   Rscript DecodeME_main.R
   ```

## INFO column approximation

The original summary statistics do not include an INFO score. The pipeline estimates a proxy INFO value using the following approximation:

\[
\text{INFO\_proxy} = \frac{1}{\text{SE}^2 \times N_{\text{eff}} \times 2p(1-p)}
\]

where \( p \) is the minor allele frequency and \( N_{\text{eff}} = N \times \pi \times (1-\pi) \) with \( \pi = N_{\text{cases}} / N \). Values are truncated to lie within the range [0, 1]. This approximation is employed because INFO values are not provided in the input data.

## License

This project is distributed without an explicit license. Please contact the authors for reuse permissions.
+   ```
+
+## INFO column approximation
+
+The original summary statistics do not include an INFO score. The pipeline estimates a proxy INFO value using the following approximation:
+
+\[
+\text{INFO\_proxy} = \frac{1}{\text{SE}^2 \times N_{\text{eff}} \times 2p(1-p)}
+\]
+
+where \( p \) is the minor allele frequency and \( N_{\text{eff}} = N \times \pi \times (1-\pi) \) with \( \pi = N_{\text{cases}} / N \). Values are truncated to lie within the range [0, 1]. This approximation is employed because INFO values are not provided in the input data.
+
+## License
+
+This project is distributed without an explicit license. Please contact the authors for reuse permissions.
 
EOF
)
