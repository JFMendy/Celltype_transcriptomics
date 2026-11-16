## Comparative transcriptomic analysis of marker expression in different cell types using bulk RNA-seq and single-cell profiling
> **Master of Science (MSc) Thesis Project| Life Science Informatics**  
> **Author:** Joseph F Mendy  
> 
> **Institution:** Deggendorf Institute of Technology

### Project Overview:
This project investigates the transcriptomic characteristics of nine peripheral blood mononuclear cell (PBMC) immune-cell populations using bulk RNA sequencing (bulk RNA-seq), followed by further characterisation using single-cell RNA sequencing (scRNA-seq) and fluorescence-activated cell sorting (FACS)

The primary objective was to investigate the potential of bulk RNA-seq as a tool for characterising distinct immune-cell populations based on their gene-expression profiles

Single-cell transcriptomic analysis provided additional resolution of cellular heterogeneity and enabled the examination of cell-type-specific gene-expression patterns. FACS provided an additional experimental approach to characterising immune-cell populations using surface markers

Single-cell transcriptomic analysis provided additional resolution of cellular heterogeneity and enabled the examination of cell-type-specific gene-expression patterns. FACS provided another approach to characterising immune-cell populations using surface markers.

The project integrates experimental immunology with computational transcriptomics to investigate immune-cell identity and gene-expression signatures.

### Research Objectives 
1. Characterise nine PBMC immune-cell populations using bulk RNA-seq.
2. Identify differentially expressed genes and population-specific molecular signatures.
3. Further investigate immune-cell heterogeneity and marker expression using scRNA-seq.
4. Examine the correspondence between transcriptomic findings and established immune-cell markers.
5. Do the identified markers encode surface proteins used in the FACS panel? 

### Immune cells of interest
The bulk RNA-seq analysis investigated nine immune-cell populations:
-- B cells
-- Conventional dendritic cells type 1 (cDC1)
-- Conventional dendritic cells type 2 (cDC2)
-- Dendritic cells type 3 (DC3)
-- Monocytes type 1 (mo1)
-- Monocytes type 2 (mo2)
-- Intermediate monocytes (moint)
-- Plasmacytoid dendritic cells (pDC)
-- Transitional dendritic cells (tDC)

The single-cell analysis enabled further examination of cellular heterogeneity and immune-cell subpopulations.
### Computational Workflow
1. Bulk RNA-seq Analysis
   - Expression data preparation and quality assessment
   - Gene filtering and data preprocessing
   - Variance-stabilising transformation
   - Principal component analysis (PCA)
   - Differential gene expression analysis using DESeq2
   - Identification of population-specific gene-expression signatures
   - Visualisation of differentially expressed genes

2. Single-Cell RNA-seq Analysis
   - Quality control and filtering
   - Doublet detection using scDblFinder
   - Data normalisation and highly variable gene selection
   - Dimensionality reduction using PCA and UMAP
   - Cell clustering and population characterisation
   - Identification of differentially expressed genes and cell-type markers
   - Visualisation of cellular heterogeneity

3. Experimental Characterisation
   - FACS was used as an additional experimental approach to examine immune-cell populations through established surface-   marker expression.

### Technologies and Tools






### Data Files description:
- The processed bulk RNA-seq data is stored as a SummarizedExperiment object
- The single-cell RNA-seq data is stored as a Seurat object
- Both objects contain expression data and associated metadata required for downstream analysis.
- The processed datasets used in the analyses are not found on GitHub but on Zenodo. Only the scripts and plots are on GitHub.
- Also, neither dataset is tracked by git.  
  
### Data Access and Placement
Due to the large size of the scRNA-seq dataset, both the bulk RNA-seq and scRNA-seq datasets are hosted on Zenodo rather than GitHub. 
The data can be downloaded from this DOI: "**[10.5281/zenodo.21398652](https://zenodo.org/uploads/21398652)**".

## Download & Setup Instructions
1. Download the following processed data objects from the Zenodo repository link above:
- se.rds: processed bulk RNA-seq data as a SummarizedExperiment object
- hsobj_20221209_annotated.rds: Processed scRNA-seq data stored as a Seurat object (approx. 7 GB)

2. Place the downloaded files into their respective subfolders within the data/ directory

- data/bulkRNASeq/se.rds
- data/scRNASeq/hsobj_20221209_annotated.rds

### R environment setup:
#### System & Software Prerequisites
- R version 4.4.3
- RStudio (recommended)
- Git & GitHub
- Memory: High RAM capacity (at least 16 GB+ recommended) to load and process the 7 GB single-cell Seurat object
  
#### Dependencies management
Package dependencies are managed via renv to ensure full environment reproducibility.
1. Clone this repository:
- git clone [https://github.com/JFMendy/Celltype_transcriptomics.git](https://github.com/JFMendy/Celltype_transcriptomics.git)
  
2. Open MSc_thesis_analysis.Rproj in RStudio
3. Restore the R package environment by executing the following in your R console:
  
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::restore()

#### Additional Rendering Dependencies
The gt package and its rendering workflow may require additional dependencies, including webshot2, chromote, and a compatible Chrome/Chromium installation. If these dependencies are not available in the project environment, they can be installed using:

renv::install(c("gt", "webshot2", "chromote", "websocket", "AsioHeaders"))
renv::snapshot()

## Execution Workflow
The analysis pipeline must be executed sequentially
1. Run Bulk RNA-Seq Analysis First:
   
source("codes/bulk_rscript.R")

**Note:** This script creates an intermediate object in the global environment and saves a copy inside the result/ folder. This intermediate object is strictly required for the cross-methodological validation in the single-cell pipeline.

2. Run Single-Cell RNA-Seq Analysis Second:

source("codes/scrna_rscript.R")

**Note:** Running scrna_rscript.R without first executing bulk_rscript.R will cause an error due to missing intermediate variables

## Output Files

- __Bulk RNA-seq plots & tables__: Saved to plots/bulk_plot/

- __Single-cell RNA-seq plots & tables:__ Saved to plots/scrna_plot/

- __Intermediate objects__: Saved to result/

## Core Packages Used
***Bulk RNA-seq Analysis***
- __SummarizedExperiment:__ stores and manages all information (metadata) and the expression matrix for bulk RNA-seq analysis
- __DESeq2:__ performs differential gene expression in bulk RNA-seq data
- __biomaRt:__ retrieves gene annotations from the Ensembl online database

**Single-cell Analysis**
- __Seurat:__ Stores and manages single-cell data, metadata and analysis results and visualisations
- __presto:__ used to calculate single-cell gene expression analysis

***Data Manipulation***
__tidyverse (dplyr, tidyr and tibble):__ used for data manipulation, transformation or wrangling

***Visualisation***
- __ggplot2:__ used to generate graphics for visualisation
- __patchwork:__ combines many plots into a single layout
- __complexHeatmap:__ used to customise the heatmap for gene expression analysis
- __circlize:__ provides colour mapping and annotation for ComplexHeatmap
- __factoextra:__ used to visualise principal component analysis (PCA) results and other multivariate analyses
- __ggrepel:__ used for adding text labels in plots
- __viridis:__ used to provides plot colours
- __scales:__ used for formatting axes and scaling of colours for plots
- __grid:__ provide graphic functions for plots
- __gridExtra:__ manages multiple plots

***Tables***
- __gt:__ used to generate summary tables
- __kableExtra:__ enhances the summary table in HTML and PDFs
- __knitr:__ used for generating reports

**Related Resources**
In a nutshell: the source code can be found on the GitHub repository; the software release and data are on Zenodo, as stated below. 

 __Source code:__
https://github.com/JFMendy/MSc_thesis_analysis

__Analysis software release:__
Zenodo DOI: 10.5281/zenodo.22211143

__Processed bulk and single-cell RNA-seq data:__
Zenodo DOI: 10.5281/zenodo.21398652
