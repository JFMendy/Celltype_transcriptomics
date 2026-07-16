# MSc_thesis_analysis
## Life Science MSc thesis project

### Data Files description:
- The processed bulk RNA-seq data is stored as a SummarizedExperiment object
- The single-cell RNA-seq data is stored as Seurat object
- Both objects contains expression data and associated metadata required for downstream analysis.
- The proceessed datasets used in the analyses are not found in GitHub but Zenodo. Only the scripts and plots are in GitHub. 
  
### How to access the used data for this analysis
Both files can be downloaded from Zenodo. Please download the data from this urls: "https://zenodo.org/uploads/21396480"
When the files are downloaded place them in the "data" directory of the GitHub repository

"se.rds"  => bulk RNA-seq SummerisedExperiment object
"hsobj_20221209_annotated.rds" => scRNA-seq Seurat object.

### File path: 
Following the downloads, file paths would be: 
"data/bulkRNASeq/se.rds"   
"data/scRNASeq/hsobj_20221209_annotated.rds"  

### Running the analyses: 
- The seurat object is large (9 GB) hence, would need a large memeory (RAM) to load the data.
- Note, first run the "bulk_rscript.R" then the "scrna_script.R"
- bulk RNA-seq analysis was performed by: scripts/bulk_rscript.R
- Single-cell RNA-seq analysis: scripts/scrna_rscript.R


