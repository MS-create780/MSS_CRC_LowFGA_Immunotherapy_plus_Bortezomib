# MSS_CRC_LowFGA_Immunotherapy_plus_Bortezomib
Code and analysis pipeline for the study: "Microsatellite-stable Colorectal Cancer with Low Chromosomal Instability Responds to Immunotherapy and Enhanced by Proteasome Inhibition"

This repository contains the R scripts used for data analysis in the study.

## ## 1. Overview
The scripts provided here cover the bioinformatic analysis of various experimental models used in our study, including:
* **Single-cell RNA sequencing (scRNA-seq)**
* **Patient-derived tissue samples**
* **Patient-derived organoid models**
* **Colorectal cancer cell lines**

## ## 2. Data Availability
The raw sequencing data have been deposited in the China National GeneBankSequence Archive (CNSA) database. 

Due to the large file size, the processed Seurat objects for single-cell data are not hosted in this GitHub repository. These data are **available from the first author** (Please contact: mengyuanshipku@gmail.com).

## ## 3. Code Structure & Modules
The analysis scripts are organized by experimental models and analysis types. There is no strict execution order; each script can be reviewed for its specific logic and parameters.

## ## 4. Environment & Dependencies
These scripts were developed and tested in
R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26200)
Matrix products: default

Key R packages used:
* `Seurat_5.2.1` 
* `dplyr_1.1.4`
* `ggplot2_3.5.1 `
* `reshape2_1.4.4`
* `DESeq2_1.42.1`
* 'clusterProfiler_4.10.1'

## ## 5. Contact
For any inquiries regarding the code or analysis, please contact:
Mengyuan Shi, Key laboratory of Carcinogenesis and Translational Research, Department of Gastrointestinal Surgery III, Peking University Cancer Hospital & Institute, mengyuanshipku@gmail.com.
