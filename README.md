# TCGA Survival Analysis Script

This TCGA Survival script utilizes the Kaplan-Meier estimate, a crucial methodology in disease research, to estimate the probability of patients undergoing defined clinical events. Specifically, it focuses on Overall Survival for lung cancer and Disease-Free Survival for prostate cancer.

## Gene-Specific Analysis

The script categorises cancer patients based on the copy number alteration (CNA) of your gene of interest.The expression levels of these genes, measured through CNA, serve as key factors in patient grouping.

In the example below the genes CCT3 and PTEN are analysed in the context of Lung and Prostate cancers respectively.

## Table of Contents
1. [Requirements](#requirements)
2. [Project Structure](#project-structure)
3. [Pipeline](#pipeline)
4. [How to Run](#how-to-run)
5. [Input Data](#input-data)
6. [Results](#results)
7. [Citation](#citation)

## Requirements

Before running the script, ensure you download the necessary data files for your chosen cancer type from [cBioPortal](https://www.cbioportal.org/) and place them in your working directory:

- **data_clinical_patient.txt:** Contains survival information for patients.
- **data_clinical_sample.txt:** Provides patient-level information, crucial for CNA data processing.

## Project Structure

.
├── OS_analysis.R
├── DFS_analysis.R
├── example_input/
│   ├── data_clinical_patient.txt
│   ├── data_clinical_sample.txt
│   └── data_dna.txt
├── results/
│   ├── coxph_summary.tsv
│   └── km_plot.png
└── README.md

## How to Run

1. Download the required data files mentioned above.
2. Set your working directory.
3. Run the function within the script to perform survival analysis on the specified cancer type, considering gene-specific copy number alterations of your gene of interest.

## Results

# PTEN loss in Prostate Cancer

The tumour supressor gene PTEN is one of the most commonly lost genes in prostate cancer. Here the "DFS" script determines if loss of PTEN copy number is able to estimate if a patient with prostate cancer is likely to progress with their disease. This progression is measured by DFS as the event.

The example uses TCGA-PRAD (PanCancer Atlas: Prostate Cancer data from c-bioportal)

![PTEN PC](https://user-images.githubusercontent.com/18528125/173332063-d6286cc4-9c33-4bf9-8e85-5e4e40c49435.png)

## Citation

If you use the data used in this example, please cite the cbioportal publications:

Cerami et al. The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. Cancer Discovery. May 2012 2; 401. PubMed.

Gao et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. Sci. Signal. 6, pl1 (2013). PubMed.

de Bruijn et al. Analysis and Visualization of Longitudinal Genomic and Clinical Data from the AACR Project GENIE Biopharma Collaborative in cBioPortal. Cancer Res (2023). PubMed.

Remember also to cite the source of the data if you are using a publicly available dataset.

https://www.cbioportal.org/