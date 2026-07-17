# Comparative transcriptomic analysis of marker expression in different cell types using bulk RNA-seq and single-cell profiling
> **Master of Science (MSc) Thesis Analysis**  
> **Author:** Joseph F Mendy  
> **Contact:** joseph.mendy@stud.th-deg.de 
> **Institution:** Deggendorf Institute of Technology

### Project goal:
The aim of the project was to access the capability of bulk RNA-seq as a standalone tool to characterise immune cells. While single-cell analysis was used to ascertain the results obtained from bulk RNA-seq. 
This repository computational analysis R scripts used for these analyses.

### Data Files description:
- The processed bulk RNA-seq data is stored as a SummarizedExperiment object
- The single-cell RNA-seq data is stored as Seurat object
- Both objects contains expression data and associated metadata required for downstream analysis.
- The proceessed datasets used in the analyses are not found in GitHub but Zenodo. Only the scripts and plots are in GitHub. 
  
### How to access the used data for this analysis
Due to the large files, the data is not tracked in this GitHUb respository. 
Both files can be downloaded from Zenodo. Please download the data from this urls: "[10.5281/zenodo.21398652](https://zenodo.org/uploads/21398652)"
When the files are downloaded, "se.rds"  => bulk RNA-seq SummerisedExperiment object and "hsobj_20221209_annotated.rds" => scRNA-seq Seurat object. The subfolders of this files should be placed in the "data" folder. 

### File path: 
Following the downloads, file paths would be: 
"data/bulkRNASeq/se.rds"   
"data/scRNASeq/hsobj_20221209_annotated.rds" 

### R environment setup:
A renv was was set up for analyses reproduciblity and package management. 
To recreate these analyses: 
- open the "MSc_thesis_analysis.Rproj"
- run this code: "renv::restore()" on the R console. R would automatically install required packages and their dependencies.

### Running the analyses: 
The analyses are structured sequencially. To recreate these analyses, the scripts should run as indicated below. 
- The seurat object is large (7GB) hence, would need a large memeory (RAM) to load the data.
- Note, first run the "bulk_rscript.R" then the "scrna_script.R"
- bulk RNA-seq analysis was performed by: codes/bulk_rscript.R
- Single-cell RNA-seq analysis: codes/scrna_rscript.R

