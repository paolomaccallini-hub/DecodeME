# DecodeME Pipeline

This repository contains R scripts for processing and fine-mapping DecodeME summary statistics. The pipeline requires both R and Python and has been tested with **R 4.4.1** and **Python 3.12**.
## Overview
DecodeME is an R-based workflow for processing and fine-mapping DecodeME GWAS-1, GWAS-2, GWAS-infectious, GWAS-non-infectious, GWAS-female, GWAS-male summary statistics. It munges inputs, converts genome builds, and retrieves linkage disequilibrium (LD) information for downstream analyses. The pipeline has been tested with **R 4.4.1** and **Python 3.12**.

This pipeline represents a replication of the results discussed in the preprint released by [DecodeME](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-). 

## Repository structure

- `DecodeME_main.R` – orchestrates filtering, format conversion, and fine mapping.
- `DecodeME_func.R` – helper functions and package setup.
- `DecodeME_config.yml` – configuration file specifying filters, API keys, and phenotype settings.
- `My_genes_DEcodeME.csv` – example output file created by `DecodeME_main.R`.

## Data sources

- **Summary statistics** – bundled in `Data/DecodeME/DecodeME_summary.zip` and
  derived from DecodeME (https://osf.io/rgqs3/files/osfstorage). The pipeline unpacks the relevant phenotype file
  before processing.
- **LD matrices** – fetched on demand from the UK Biobank LD reference panel
  hosted by the Broad Institute at
  `https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/` and converted to
  RDS files for fine-mapping.
- **RegulomeDB scores** – obtain from <https://www.regulomedb.org/regulome-search/> and place the score file in the repository root.

## Requirements
- R (tested with version 4.4.1)
- Python (tested with version 3.12)
- At least 10 GB of free disk space
- RegulomeDB score file downloaded from [https://www.regulomedb.org/regulome-search/](https://www.regulomedb.org/regulome-search/) and placed in the repository root
- The internal folder structure is generated automatically by the scripts

## Setup
The internal directory tree is created automatically by `DecodeME_main.R`.

1. Install [R](https://www.r-project.org/) and [Python](https://www.python.org/) along with the required R packages. The scripts make use of Bioconductor resources such as `SNPlocs.Hsapiens.dbSNP155.GRCh38` and `BSgenome.Hsapiens.NCBI.GRCh38`.
2. Download the RegulomeDB score file and place it in the repository root.
3. Adjust `DecodeME_config.yml` to provide paths, API keys, and filtering thresholds.
4. Run the main script (it will automatically build the required folder hierarchy):
3. Configure paths and filters in `DecodeME_config.yml`.
4. Run:
   ```bash
   Rscript DecodeME_main.R
   ```

## Munging and Liftover

`DecodeME_main.R` standardises the summary statistics with
[`format_sumstats`](https://github.com/neurogenomics/MungeSumstats), which
munges allele columns and **lifts coordinates from GRCh38 to GRCh37**. The
conversion ensures that downstream steps operate on GRCh37 positions, matching
the genome build used by the UK Biobank LD matrices.
`DecodeME_main.R` invokes [`format_sumstats`](https://github.com/neurogenomics/MungeSumstats) to standardize column names, align alleles, and **convert coordinates from GRCh38 to GRCh37** to match the UK Biobank LD reference.

## Output 

This analysis generates two supplementary loci (on chr 10 and chr 15) that are not present in the results by [DecodeME](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-). This may be due to the lift-over from GRCh38 to GRCh37 that was performed on the summary statistics to allow for the use of available UKB LD matrices. I manually removed these two loci from the results reported below.

### Output 1: fine-mapping
The following images represent the output of fine-mapping on all the DecodeME cohorts (GWAS-1, GWAS-2, GWAS-infectious, GWAS-non-infectious, GWAS-female, GWAS-male). Posterior Inclusion Probability (PIP) is reported on the y-axis, and it represents the probability of being a causal variant. Credible sets are defined as sets of variants whose PIPs sum up to 95%. The legend indicates the number of credible sets, the size of each credible set, and the mean linkage disequilibrium (measured as |R|) between each possible pair of variants from the same credible set. Coordinates on the x-axis are with respect to GRCh37!

<img width="1000" height="500" alt="gwas_1_infectious_onset_1" src="https://github.com/user-attachments/assets/76a270d8-2f16-49f0-aa3d-94c0fbc198a1" />
<img width="1000" height="500" alt="gwas_1_infectious_onset_2" src="https://github.com/user-attachments/assets/0c061aa5-3bcc-457c-9150-92bf787a32e2" />
<img width="1000" height="500" alt="gwas_1_infectious_onset_3" src="https://github.com/user-attachments/assets/ad0ea391-c2f9-48ed-b4c8-a64a8d632fd0" />
<img width="1000" height="500" alt="gwas_1_female_1" src="https://github.com/user-attachments/assets/45a2ca8a-f6e0-43c9-9ff2-4d00e0e5ed89" />
<img width="1000" height="500" alt="gwas_1_female_2" src="https://github.com/user-attachments/assets/5bd98ad8-cc28-48a5-832f-516eda00dadd" />
<img width="1000" height="500" alt="gwas_1_female_3" src="https://github.com/user-attachments/assets/15a4404e-e195-44a8-9797-7213395d080c" />
<img width="1000" height="500" alt="gwas_2_1" src="https://github.com/user-attachments/assets/cd2101c7-776b-44a4-bf19-0c0ed5a11720" />
<img width="1000" height="500" alt="gwas_2_2" src="https://github.com/user-attachments/assets/93ccefb8-1b62-4da3-8e4d-b8c36a243c4a" />
<img width="1000" height="500" alt="gwas_2_3" src="https://github.com/user-attachments/assets/a8cfc424-63da-444e-8742-ae5bb3d01165" />
<img width="1000" height="500" alt="gwas_2_4" src="https://github.com/user-attachments/assets/efd34825-230d-42da-ac72-77bb287e057c" />
<img width="1000" height="500" alt="gwas_1_1" src="https://github.com/user-attachments/assets/711bf8fc-3292-44b9-be76-2f7170de3189" />
<img width="1000" height="500" alt="gwas_1_2" src="https://github.com/user-attachments/assets/77d7bc6d-2786-4484-ac71-3dd16a55a09a" />
<img width="1000" height="500" alt="gwas_1_3" src="https://github.com/user-attachments/assets/d3047c90-d1ca-4d6c-bab9-13969dd4900e" />
<img width="1000" height="500" alt="gwas_1_4" src="https://github.com/user-attachments/assets/b41de0ba-dc2f-46fd-b097-759c6b974a71" />
<img width="1000" height="500" alt="gwas_1_5" src="https://github.com/user-attachments/assets/ee92b447-2f2f-464e-b91e-d0a12446ed4d" />

### Output 2: gene-mapping
The first rows of the output file `My_genes_DecodeME.csv` are reported below. We find the following 18 genes: ABT1, ANKRD45, ARFGEF2, BTN2A2, CSE1L, DARS2, KLHL20, PRDX6, RABGAP1L, RC3H1, SERPINC1, SLC9C2, STAU1, TNFSF4, TRIM38, ZBTB37, ZNFX1, OLFM4. This is a subset of the 32 genes proposed as candidate genes in the [DecodeME preprint](https://www.research.ed.ac.uk/en/publications/initial-findings-from-the-decodeme-genome-wide-association-study-).  

### GWAS-associated eQTL gene list (excerpt)

| name    | NCBI.id | weight | Phenotype     | Cases | Description | Variant    | Tissues                                                                                                                                                              | 
|---------|---------|--------|---------------|-------|-------------|------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ABT1    | 29777   | 1      | gwas_1        | 15579 | eQTL        | rs10946808 | Lung / Cells_Cultured_fibroblasts                                                                                                                                    |
| ABT1    | 29777   | 1      | gwas_1_female | 12833 | eQTL        | rs61534839 | Muscle_Skeletal                                                                                                                                                      |
| ABT1    | 29777   | 1      | gwas_1_female | 12833 | eQTL        | rs9379845  | Muscle_Skeletal                                                                                                                                                      | 
| ANKRD45 | 339416  | 1      | gwas_1        | 15579 | eQTL        | rs1322775  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta |
| ANKRD45 | 339416  | 1      | gwas_1        | 15579 | eQTL        | rs1322777  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta | 
| ANKRD45 | 339416  | 1      | gwas_1        | 15579 | eQTL        | rs2065171  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta | 
| ANKRD45 | 339416  | 1      | gwas_1        | 15579 | eQTL        | rs9425435  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta | 
| ANKRD45 | 339416  | 1      | gwas_1        | 15579 | eQTL        | rs9425757  | Skin_Not_Sun_Exposed_Suprapubic / Skin_Sun_Exposed_Lower_leg / Artery_Tibial / Adipose_Subcutaneous / Cells_Cultured_fibroblasts / Thyroid / Breast_Mammary_Tissue / Artery_Aorta / Nerve_Tibial |
| ANKRD45 | 339416  | 1      | gwas_2        | 15579 | eQTL        | rs1322775  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta |
| ANKRD45 | 339416  | 1      | gwas_2        | 15579 | eQTL        | rs1322777  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta |
| ANKRD45 | 339416  | 1      | gwas_2        | 15579 | eQTL        | rs2065171  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta |
| ANKRD45 | 339416  | 1      | gwas_2        | 15579 | eQTL        | rs9425435  | Artery_Tibial / Cells_Cultured_fibroblasts / Nerve_Tibial / Thyroid / Skin_Not_Sun_Exposed_Suprapubic / Breast_Mammary_Tissue / Adipose_Subcutaneous / Skin_Sun_Exposed_Lower_leg / Artery_Aorta |
| ANKRD45 | 339416  | 1      | gwas_2        | 15579 | eQTL        | rs9425757  | Skin_Not_Sun_Exposed_Suprapubic / Skin_Sun_Exposed_Lower_leg / Artery_Tibial / Adipose_Subcutaneous / Cells_Cultured_fibroblasts / Thyroid / Breast_Mammary_Tissue / Artery_Aorta / Nerve_Tibial |

## FUMA/MAGMA

The summary statistics generated after filtering with respect to INFO and after lift-over and munging can be used as input for [FUMA](https://fuma.ctglab.nl/), setting UKB release2b White British as reference population, positional mapping and eQTL mapping (GTEx v8) as gene mapping algorithms. This analysis is publicly available [here](https://fuma.ctglab.nl/browse/663962). FUMA will select risk loci and map them to genes. In this case, it selects a total of 45 genes. Note that FUMA does not apply any fine-mapping method, and this explains why it retrieves a high number of genes.

### Tissie enrichment 

One of FUMA's features is tissue-enrichment analysis against GTEx v.8 (54 tissues), as can be seen in the histogram below.

<img width="1145" height="598" alt="image" src="https://github.com/user-attachments/assets/b51c6be6-71b1-4997-8bd8-e176f66f9c68" />

### Cell enrichment 

Given the prevalence of brain tissues in the previous analysis, I proceed with cell-enrichment (Window Cell Type on FUMA), using human and mouse brains as datasets. The main output is reported below, with a description of cell types and samples.

<img width="697" height="432" alt="Screenshot 2025-10-02 172713" src="https://github.com/user-attachments/assets/72481e74-c8a3-4fde-9f60-1f84a9d36fe8" />

| Cell Type | Species | Brain Region | Cortical Layer | Biological Role |
|-----------|---------|---------------|----------------|------------------|
| **Neuron (broad)** | Human | Developing cortex (Carnegie Stage 18) | All | General neurons; early neurodevelopmental stages |
| **PC.Neuron_Sc17a7_Fermt1.2_12** |  Mouse | Whole Brain | Layers II/III–IV  | Excitatory glutamatergic pyramidal neurons |
| **L4_Plcxd2** | Mouse | Cortex / hippocampus | Layer IV | Layer IV excitatory projection neurons |
| **S1PyrL4** | Mouse | Primary somatosensory cortex (S1) | Layer IV | Layer IV pyramidal neurons |
