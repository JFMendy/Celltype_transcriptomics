## Comparative transcriptomic analysis of marker expression in different cell types using bulk RNA-seq and single-cell profiling
> **Master of Science (MSc) Thesis Analysis**  
> **Author:** Joseph F Mendy  
> **Contact:** joseph.mendy@stud.th-deg.de
> 
> **Institution:** Deggendorf Institute of Technology

### Project goal:
This repository contains my MSc thesis in Life Science Informatics analysis documentation. The project includes two sets of data: bulk RNA sequencing (bulk RNA-seq) and single-cell RNA sequencing (scRNA-seq).
The aim of the project was to access the capability of bulk RNA-seq as a standalone tool to characterise immune cells. While single-cell analysis was used to ascertain the results obtained from bulk RNA-seq. 
This repository also included computational analysis R scripts used for these analyses.

### Data Files description:
- The processed bulk RNA-seq data is stored as a SummarizedExperiment object
- The single-cell RNA-seq data is stored as Seurat object
- Both objects contains expression data and associated metadata required for downstream analysis.
- The proceessed datasets used in the analyses are not found in GitHub but Zenodo. Only the scripts and plots are in GitHub.
- Also both datasets are not tracked by git.  
  
### How to access the used data for this analysis
Due to the large scRNA-seq file, both files were stored in Zenodo. And could be downloaded from this url: "**[10.5281/zenodo.21398652](https://zenodo.org/uploads/21398652)**". When the files are downloaded, "se.rds"  => bulk RNA-seq SummerisedExperiment object and "hsobj_20221209_annotated.rds" => scRNA-seq Seurat object. The subfolders of this files should be placed in the "data" folder in repository. 

### File paths: 
Following the downloads, file paths would be: 
- se.rds => ***data/bulkRNASeq/se.rds***  
- hsobj_20221209_annotated.rds => ***data/scRNASeq/hsobj_20221209_annotated.rds***
- bulk RNA-seq analysis source code => ***codes/bulk_rscript.R***
- Single-cell RNA-seq analysis => ***codes/scrna_rscript.R***

### R environment setup:
A renv was was set up for analyses reproduciblity and package management. 
To recreate these analyses: 
- use this code and link: ***git clone https://github.com/USERNAME/MSc_thesis_analysis.git***
- open the ***MSc_thesis_analysis.Rproj***
- Run the following codes on the R console:
  - ***install.packages("renv")***
  - ***renv::restore()***.
    
R would automatically install required exact versions of packages and their dependencies used in this project stored in the renv.lock file.

### Running the analyses: 
The analyses are structured sequencially. To recreate these analyses, the scripts should run as indicated below. 
- The seurat object is large (7GB) hence, would need a large memeory (RAM) to load the data.
- Note, first run the "bulk_rscript.R" then the "scrna_script.R"
- Each script also contain custom helper functions to perform repetitive analyses.

### Ouput Files 
During the running of the source codes, each script automatically generate's figures, tables and diffferential expression results. The generated plots and tables are stored as indicalted. 
- The bulk script would store its output in  => **plots/bulk_plot/**
- While the single-cell script => **plots/scrna_plot/**

**Note:** There is a generated output intermediate object by bulk RNA-seq source code required to compare single-cell analysis results against bulk RNA-seq result. This intermediate object is stored in the R global environment, hence the need to run the bulk RNA-seq source code first otherwise, the single-cell would generate an error during the course of running. 

### Software requirement
- R Version: 4.4.3
- Rstudio version:
- CRAN packages
- Bioconductor packages
- renv for dependency management
- Git and Git Hub

### Packages description and dependencies 

The packages were obtained from CRAN and Bioconductor. The Bioconductor packages were also managed through **BiocManager** as well renv.lock file.
- ***Dependencies:*** 
The following packages require dependencies as indicated by the arrow: 
gt packages needs => webshot2 => chromote => chrome(google chrome) /chromium => websocket => AsioHeaders and zip  to work properly. Therefore, this code would be needed if the above mentioned packages are not working as expected. Then run these two codes to get it running as expected. 
- __renv::install(c("gt",
"webshot2",
"chromote",
"websocket",
"AsioHeaders"
))__
- __renv::snapshot()__


***reproducibility***
- __renv:__ Manages the project-specific R environment and package versions

***Bulk RNA-seq Analysis***
- __SummarizedExperiment:__ stores and manages all information (metadata) and expression matrix for bulk RNA-seq analysis
- __DESeq2:__ perform the differential gene expression in bulk RNA-seq data
- __biomaRt:__ retrieves gene annotations from Ensembl online database

**Single-cell Analysis**
- __Seurat:__ Stores and manages single-cell data, metadata and analysis results and visulisations
- __presto:__ used to calculate single-cell gene expression analysis

***Data Manipulation***
__tidyverse (dplyr, tidyr and tibble):__ used for data manipulation, transformation or wrangling

***Visualisation***
- __ggplot2:__ used to generate graphics for visualisation
- __patchwork:__ combines many plots into a single layout
- __complexHeatmap:__ used to customised the heatmap for gene expression analysis
- __circlize:__ provides colour mapping and annotation for ComplexHeatmap
- __factoextra:__ used to visulaise principal component analysis (PCA) results and other multivariate analyses
- __ggrepel:__ used for adding text labels in plots
- __viridis:__ used to provides plot colours
- __scales:__ used for formatting axesand scaling of colours for plots
- __grid:__ provide graphic functions for plots
- __gridExtra:__ manages multiple plots

***Tables***
- __gt:__ used to generate summary tables
- __kableExtra:__ enehance the summary table in HTML and pdfs
- __knitr:__ used for generating reports


