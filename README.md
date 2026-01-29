# Group-wise-Hybrid-KNN-MinProb-Imputation-Pipeline-for-LFQ-Proteomics-Data

This repository contains an analysis pipeline developed for the processing and statistical analysis of label-free quantitative (LFQ) proteomics data.

The pipeline is centered on a **sequential, group-aware missing value imputation strategy (KNN followed by MinProb)**, designed to robustly address both stochastic and low-abundance missingness commonly observed in LFQ proteomics datasets.

The workflow was developed to identify differentially abundant proteins in brown adipose tissue (BAT) from mouse models, comparing wild-type (WT) and *Trpa1*⁻/⁻ animals exposed to prolonged mild cold conditions.

## Overview

The analysis workflow includes:

* Preprocessing and filtering of LFQ proteomics data
* Log2 transformation and handling of zero-intensity values
* **Sequential missing value imputation using KNN and MinProb approaches**
* Group-wise statistical testing for differential protein abundance
* Log2 fold-change calculation and significance classification
* Integration of quantitative results with protein annotations

All steps were implemented with a strong emphasis on reproducibility and adherence to best practices in LFQ proteomics analysis.

## Experimental Design

* Tissue: Brown adipose tissue (BAT)
* Organism: Mouse (*Mus musculus*)
* Genotypes:

  * Wild-type (WT)
  * *Trpa1*⁻/⁻
* Temperature conditions:

  * Thermoneutrality (30ºC)
  * Mild Cold (22ºC)
* Biological replicates: 3 per group
* Quantification method: Label-free quantification (LFQ)

## Missing Value Imputation Strategy

Missing values were addressed using a **two-step, sequential imputation strategy applied independently within each experimental group**, allowing distinct mechanisms of missingness to be modeled explicitly.

### Step 1: KNN Imputation (Partially Observed Proteins)

Proteins quantified in **2 out of 3 biological replicates within a group** were imputed using a K-nearest neighbors (KNN) approach (`impute.knn`).

* Applied independently for each experimental group
* Restricted to proteins with sufficient observed values (≥2 replicates)
* The number of neighbors (k) was optimized per group by minimizing the root mean square error (RMSE) through repeated missing value simulations

This step primarily targets missing values arising from **stochastic sampling effects (missing at random, MAR)** rather than true biological absence.

### Step 2: MinProb Imputation (Low-Abundance Signals)

Remaining missing values were imputed using a **minimal probability (MinProb)** approach (`impute.MinProb`).

* Applied only to proteins detected in at least one replicate within a group
* Proteins completely absent across all replicates of a group were left as missing
* Imputed values were drawn from the lower tail of the intensity distribution (q = 0.01)

This step models missing values consistent with **low-abundance proteins below the detection limit (missing not at random, MNAR)**.

### Reproducibility Considerations

To ensure reproducibility of the imputation procedures:

* Protein ordering was fixed prior to imputation
* Random seeds were explicitly defined
* All imputation steps were performed deterministically within experimental groups

## Statistical Analysis

Differential protein abundance between temperature conditions within each genotype was assessed using Student’s *t*-test applied to log2-transformed, imputed intensities.

Log2 fold changes were calculated as differences between group means. Proteins were classified as significantly regulated based on combined thresholds for statistical significance and effect size.

## Repository Structure

* `data/` – Input LFQ data files generated with PEAKS Studio v11 (PEAKS Q mode)
* `scripts/` – Analysis scripts
* `results/` – Output tables and figures
* `docs/` – Additional documentation

## Notes

This repository is provided to support transparency and reproducibility of the analyses performed in the associated study. The corresponding mass spectrometry data are publicly available in MassIVE under accession **MSV000100650** (doi:10.25345/C5J38KX4M).

