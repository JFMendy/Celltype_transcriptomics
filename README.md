# MSc_thesis_analysis
## Life Science MSc thesis project

## How to access the used data for this analysis
The bulk and single cell RNA-seq dataset are provided as GItHub Release assets.  Please download the data from latest release as: 
"se.rds"  => bulk RNA-seq SummerisedExperiment object
"seu.rds" => scRNA-seq Seurat object.

Following the downloads, file paths would be: 
data/bulkRNASeq/se.rds   for the bulk RNA-seq dataset  and 
data/scRNASeq/seu.rds    for the single cell dataset

## Running the analyses:  
bulk RNA-seq analysis was performed by : scripts/bulk_rscript.R
Single-cell RNA-seq analysis: scripts/scrna_rscript.R
