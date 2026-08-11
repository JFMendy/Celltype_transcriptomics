

#=======MSc Thesis Analysis=============================================================================================================================
# MSc Thesis Analysis: Bulk and single cell RNA-seq dataset

# Bulk RNA-seq Analysis
# 0.1 package load =====================================================================================================================================
{library(DESeq2)
  library(biomaRt)
  library(tidyverse)
  library(factoextra)
  library(gridExtra) 
  library(grid)
  library(ggrepel)
  library(gt)
  library(scales)
  library(ComplexHeatmap)
  library(patchwork)
  library(circlize)
  library(UpSetR)
  library(grid)

  
  }

# 0.2 read in the data ===================================================================================================================================
## read in the bulk RNA-seq data 

mydata_bulk_RNA <- readRDS("data/bulkRNASeq/se.rds")




# analysis workflow 


# 1. Data inspect the data and extract the count data set ================================================================================================



## 1.1 make a copy of original data---- 
mydata_bulk_RNA2 <- mydata_bulk_RNA

rm(mydata_bulk_RNA)

## 1.2 Inspect the SummerisedExperiemnt structure --------------------------------------------------------------------------------------------------------

# Check available assays
assayNames(mydata_bulk_RNA2)  

# colData column  and rowData names 
names(colData(mydata_bulk_RNA2)) 
names(rowData(mydata_bulk_RNA2))

# check the column names  
colData(mydata_bulk_RNA2)

# check the row names 
rowData(mydata_bulk_RNA2)   

# show top 6 rows 
head(colData(mydata_bulk_RNA2)) 

# dimensions should match original
dim(mydata_bulk_RNA2) 

# summary of data file 
glimpse(mydata_bulk_RNA2) 

# detail structure of the data file 
str(mydata_bulk_RNA2)

#assy names
assayNames(mydata_bulk_RNA2)



# 2. data wrangling and clean ===============================================================================================================================

# check the sample information 
colData(mydata_bulk_RNA2)
rowData(mydata_bulk_RNA2)
# remove  columns of no interest 
colData(mydata_bulk_RNA2)$...7 <- NULL
colData(mydata_bulk_RNA2)$We.have.4.donors..Donor.252..253..254..255 <- NULL

# check for duplicated ids 
any(duplicated(rowData(mydata_bulk_RNA2)$id))

# extract the duplicates and view the duplicates and view 
duplicates <- (rowData(mydata_bulk_RNA2)$id[duplicated(rowData(mydata_bulk_RNA2)$id)])  
duplicates 
rm(duplicates)
# remove the duplicate
mydata_bulk_RNA2 <- mydata_bulk_RNA2[!duplicated(rowData(mydata_bulk_RNA2)$id), ]

# confirm duplicate removal 
any(duplicated(rowData(mydata_bulk_RNA2)$id))

# check the col_data characteristics 
str(colData(mydata_bulk_RNA2))

# convert the col_data columns of interest to factors 
colData(mydata_bulk_RNA2)$Cell.population <- as.factor(colData(mydata_bulk_RNA2)$Cell.population)

colData(mydata_bulk_RNA2)$Patient <- as.factor(colData(mydata_bulk_RNA2)$Patient)

colData(mydata_bulk_RNA2)$Timepoint.after.vaccination <- as.factor(colData(mydata_bulk_RNA2)$Timepoint.after.vaccination)

# remove genes with no symbols 
mydata_bulk_RNA2 <- mydata_bulk_RNA2[
  !is.na(rowData(mydata_bulk_RNA2)$symbol), ]

# extract the genes counts 
counts <- assay(mydata_bulk_RNA2, "counts")


# remove genes with less than 10 row sum count and sample less than 3 
keep <- rowSums(counts >= 10) >= 3
mydata_bulk_RNA2 <- mydata_bulk_RNA2[keep, ]

# check if the row and column names match
all(colnames(assay(mydata_bulk_RNA2, "counts")) == 
      rownames(colData(mydata_bulk_RNA2)))


rm(keep)


# 3 exploratory analysis and visualisation =====================================================================================================================

# create a saving function 
# wrap  saving function around ggsave 
save_plot_with_ggsave <- function(
    plot_obj,
    filename,
    width,
    height,
    units = "in",
    folder = "plots/bulk_plot"
) {
  
  # Create folder if needed
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  ggsave(
    filename = file.path(folder, filename),
    plot = plot_obj,
    width = width,
    height = height,
    units = units
  )
  
  message("Saved: ", file.path(folder, filename))
}


## 3.1 Exploratory analysis: --------------------------------------------------------------------------------------------------------------------------------------------

# Sample demographic

# extract data of interest to generate demographic data and gene count profile per sample 
meta <- as.data.frame(colData(mydata_bulk_RNA2))
head(meta)

# tabulate study participants
table(meta$Patient, meta$Timepoint.after.vaccination)

# number of samples per cell 
table(meta$Cell.population)

# gene count per total sample 
dim(mydata_bulk_RNA2)

# extract metadata and summarize
demo_table <- meta %>%
  dplyr::count(Patient, Cell.population) %>%   
  tidyr::pivot_wider(
    names_from = Cell.population,
    values_from = n,
    values_fill = 0
  ) %>%
  dplyr::mutate(total_samples = rowSums(dplyr::across(where(is.numeric)))) %>%
  dplyr::arrange(Patient) %>%
  as.data.frame()

# print the table 
head(demo_table)

# generate a row vector (totals) to contain samples per cell  
totals <- demo_table %>%
  dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
  dplyr::mutate(Patient = "Total") %>%
  dplyr::select(Patient, dplyr::everything())
totals

# join the totals to the demo_table  
demo_table_total <- dplyr::bind_rows(demo_table, totals)
demo_table_total 

# generate a more polish demographics table 
# convert standard table to a grob for plotting
# load the library(gridExtra) library(grid)

sample_demo_table <- gridExtra::tableGrob(demo_table_total, 
                                          rows = NULL, 
                                          theme = gridExtra::ttheme_default(
                                            core = list(fg_params = list(cex = 1)),
                                            colhead = list(fg_params = list(cex = 1.2, 
                                                                 fontface = "bold"
                                                                 )) ) )


# display table on a blank plot
grid::grid.newpage()
grid::grid.draw(sample_demo_table)

# to save the plot 
# wrap the plot in grid_arrange for ggplot to open it easily. 
# Create a grob object for saving
sample_demo_table_plot2 <- gridExtra::arrangeGrob(sample_demo_table)



# gene count analysis: 

# extract the gene counts information 
# generate gene count plot
# order samples by total_gene_counts for nicer visualization
# add a column called total_gene_count 


# calculate the total read per sample and save it in a column total_gene_counts 
meta$total_read_counts <- colSums(counts)    

#create a table with column names: sampleName, cell type, total counts
read_per_sample <- meta %>%
  
  # needed columns 
  dplyr::select(sampleName, Cell.population, total_read_counts) %>% 
  
  # order base on  total_read_counts
  dplyr::arrange(Cell.population,total_read_counts) %>%  
  
  # create new column name for sample
  dplyr::mutate(sampleName = factor(sampleName, levels = sampleName))   

# preview the table 
head(read_per_sample)


# remove any missing value or NAs
# prepare the data for plotting 
read_clean <- read_per_sample %>%
  filter(!is.na(total_read_counts) & !is.na(Cell.population))





# generate a plot for the gene count per sample 
read_plot <- ggplot(read_clean, 
       aes(x = sampleName, 
           y = total_read_counts, 
           color = Cell.population)) +
  
  # lollipop plot sticks
  geom_segment(aes(xend = sampleName, 
                   y = 0,
                   yend = total_read_counts),
               linewidth = 1.2) +
  
  # lollipop plot heads
  geom_point(size = 4) +
  
  # titles and axis labels
  labs(
    title = "read counts per Cell Type",
    x = "Sample ID",
    y = "Total read counts per Sample",
    color = "Cell Population"
  ) +
  
  # y-axis with comma formatting and 6 pretty breaks
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::comma
  ) +
  
  # color palette for many cell types
  viridis::scale_color_viridis(discrete = TRUE,
                               option = "D"
                               ) +
  
  # select theme: subtle grids, boxed axes
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, 
                              size = 25, 
                              face = "bold"),
    
    # text size and direction 
    axis.title = element_text(size = 20, 
                              face = "bold"),
    
    axis.text.y = element_text(size = 15, 
                               face = "bold"),
    
    axis.text.x = element_text(angle = 90, 
                               hjust = 1, 
                               size = 10),
    
    # subtle grid lines
    panel.grid.major = element_line(color = "gray85",
                                    linewidth = 0.3),
    
    panel.grid.minor = element_line(color = "gray90", 
                                    linewidth = 0.15),
    
    # boxed axes
    axis.line = element_line(color = "black", 
                             linewidth = 1),
    
    # legend styling
    legend.position = "right",
    legend.title = element_text(size = 15, 
                                face = "bold"),
    
    legend.text = element_text(size = 12,
                               face = "bold"),
    
    legend.key = element_rect(fill = "white", 
                              color = "black"),
    
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.y = unit(0.3, "cm")
  ) +
  
  # customise the legend symbols
  guides(color = guide_legend(override.aes = list(size = 4),  # bigger legend points
      ncol = 1                        # number of columns
    )
  )

read_plot


# arrange lollipop plot in ascending order 
# extract the count matrix 
df_counts <- assay(mydata_bulk_RNA2, "counts")

head(df_counts)

# turn into a tibble in long format
df_long <- as_tibble(df_counts, rownames = "gene") %>% 
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "count"
  )

head(df_long)

# Compute QC metrics per sample
df_detected_qc_metrics <- df_long %>%
  group_by(sample) %>%
  
  summarize(library_size = sum(count), 
            
            # total counts per sample
            # number of genes detected
            n_genes = sum(count > 0), 
            
            # number of zero-count genes
            n_zero = sum(count == 0),
            
            # percent of genes with zero count
            pct_zero = mean(count == 0) * 100, 
            
            # mean count per gene
            mean_count = mean(count), 
            
            # median count per gene
            median_count = median(count),          
            .groups = "drop"
            )

# Library Size Bar plot: Samples with lower total counts may indicate low sequencing depth or a problem with that sample
df_detected_qc_metrics

# convert to dataframe 
df_detected_qc_metrics <- as.data.frame(df_detected_qc_metrics)

# convert the sample column to a factor levels 
df_detected_qc_metrics$sample<- as.factor(df_detected_qc_metrics$sample
                                          )
# check the df structure 
str(df_detected_qc_metrics)

# load the forcats library 
# order the plots based on gene counts
df_detected_qc_metrics$sample <- fct_reorder(df_detected_qc_metrics$sample, 
                                             df_detected_qc_metrics$n_genes)



# generate display as Lollipop plot 
df_detected_qc_metrics_plot <- ggplot(df_detected_qc_metrics, 
                                      aes(x = sample, 
                                          y = n_genes)) +
  
  # Lollipop sticks - thinner and lighter
  geom_segment(aes(xend = sample,  
                   y = 0, 
                   yend = n_genes),
               
               # line colours
               color = "black", 
               
               # line thickness
               linewidth = 0.2) +        
  
  # Points - gradient colors with viridis scale
  # Color by gene count value
  geom_point(aes(fill = n_genes),
             # size of the point 
             size = 2, 
             
             # filled circles with outline
             shape = 21,              
             color = "black",
             
             # outline thickness 
             stroke = 0.2)  +  
  
  # plot background and customised display 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,  
                                   hjust = 1,
                                   
                                   # vertical adjusment  
                                   vjust = 1, 
                                   size = 3, 
                                   face = 'bold'),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = 25,
                                  
                                  # horizontal centering 
                                  hjust = 0.5, 
                                  margin = margin(b = 10)),
        # text size and bolden 
        axis.title.x = element_text(face = "bold", 
                                    size = 20,
                                    margin = margin(t = 10)),
        axis.title.y = element_text(face = "bold",
                                    size = 20,
                                    margin = margin(r = 15)),
        legend.title = element_text(size = 13,
                                    face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.caption = element_text(color = "gray40", 
                                    face = "italic", 
                                    hjust = 0),
        # margin round the plot 
        plot.margin = margin(15, 15, 15, 15)) + 
  
  # y-axis scaling 
  scale_y_continuous(
    limits = c(0, 
               max(df_detected_qc_metrics$n_genes) * 1.1),
    expand = c(0, 0),
    # Format numbers with commas
    labels = scales::comma  
  )+
  
  # plot labels
  labs(
    title = "Gene Detection Rate per Sample",
    x = "Samples (Ordered by Gene Count)",
    y = "Number of Detected Genes",
    caption = "Yellow Fever: Bulk RNA-Seq Dataset"
  )

df_detected_qc_metrics_plot 



## 3.2 call the saving function:-------------------------------------------------------------------------------------------------------
 
# call the save_plot_with_ggsave function 

# save the demographic table 
save_plot_with_ggsave(
  plot_obj = sample_demo_table_plot2,
  filename = "sample_demo_table_plot2.pdf",
  width = 12,
  height = 8
)

# save the gene count table 
save_plot_with_ggsave(
  plot_obj = read_plot,
  filename = "Read_count_plot.pdf",
  width = 40,
  height = 20
  )


# save gene count arrange in ascending order of  detected genes per sample 
save_plot_with_ggsave(
  plot_obj = df_detected_qc_metrics_plot,
  filename = "df_detected_qc_metrics_plot.pdf",
  width = 12,
  height = 8
)



# 4 Principle component analysis ==============================================================================

# Normalisation was already done before I received the data 

# extract TPM matrix
tpms <- assay(mydata_bulk_RNA2, 'tpms')

# filter out zero values 
tpms <- tpms[rowSums(tpms != 0) > 0, ]
sum(is.nan(tpms))

# add constant to avoid  negative infinity  because log2(0) give negative infinity 
# log-transform + constant value 
mat <- log2(tpms + 0.55)  

# transpose so samples are rows
mat_t <- t(mat)

# perform PCA
# again the data was also already scaled

pca_result <- prcomp(mat_t)  
summary(pca_result)

# extract PC scores + join metadata
pca_df <- as.data.frame(pca_result$x)

# create name call sampeName 
pca_df$sampleName <- rownames(pca_df)

# meta data extraction  in a df 
metadata_df <- as.data.frame(colData(mydata_bulk_RNA2))

# join the two dfs by samp name column 
pca_df <- merge(pca_df, metadata_df, by = "sampleName")  

# view the outcome of the operation 
head(pca_df)


# Plot PC1 vs PC2
# round up the % values 
percentVar <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

# provide colors for the cells 
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "black", "#AD7700", "#666666")


# display the plot in ggplot 
pca_cell_population_plot <- ggplot(pca_df, 
                               aes(x = PC1, 
                                   y = PC2, 
                                   colour = Cell.population)) +
  geom_point(size = 3) +
  
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  geom_point(size = 3) +
  
  # Add ellipses with multivariate t-distribution
  stat_ellipse(type = "t",
               
               # 95% confidence level
               level = 0.95, 
               
               # line thickness  
               linewidth = 1, 
               
               # line type
               linetype = "solid",
               
               # appearance 
               alpha = 0.5)+
  
  scale_color_manual(values = colorblind_palette) +
  theme_classic()+ 
  
  labs(title = "PCA with 95% Confidence Ellipses") +
  
  theme(plot.title = element_text(hjust = 0.5, 
                                  face = "bold", 
                                  size = 20),       
        axis.title.x = element_text(face = "bold", 
                                    size = 16),                    
        axis.title.y = element_text(face = "bold",
                                    size = 16),                    
        axis.text.x = element_text(face = "bold", 
                                   size = 14),                     
        axis.text.y = element_text(face = "bold", 
                                   size = 14),                     
        legend.text = element_text(size = 14),                                   
        legend.title = element_text(size = 16, 
                                    face = "bold"),  
        legend.position = "right"
  )
pca_cell_population_plot


# call the save_plot_with_ggsave function  to save the plot 
save_plot_with_ggsave(
  plot_obj = pca_cell_population_plot,
  filename = "pca_cell_population_plot.pdf",
  width = 12,
  height = 8
)


# generate a Scree plot with  fviz_eig function  

elbow_point <- 6

scree_plot <- fviz_eig(pca_result,
                       choice  = "variance",
                       ncp = 15,
                       geom = c("bar", "line"),
                       barfill = "gray50",
                       barcolor = "black",
                       linecolor = "red",
                       addlabels = TRUE,
                       hjust = -0.3,
                       linewidth = 10
                       ) +
  ylim(0, 20) +
  
  # horizontal threshold line
  geom_hline(
    aes(yintercept = 6, linetype = "5% variance threshold"),
    color = "blue",
    linewidth = 1
  ) +
  
  # elbow line
  geom_vline(
    aes(xintercept = elbow_point, linetype = "Elbow point (PC6)"),
    color = "green",
    size = 1
  ) +
  
  scale_linetype_manual(
    name = "Reference lines",
    values = c(
      "5% variance threshold" = "dashed",
      "Elbow point (PC6)" = "dashed"
    )
  ) +
  
  labs(
    title = "Scree Plot — Variance Explained",
    x = "Principal Components (PC)",
    y = "Percent Variance Explained"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold",
      color = "black",
      size = 12
    ),
    axis.text.y = element_text(
      face = "bold",
      color = "black",
      size = 11
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

print(scree_plot)

# call the save_plot_with_ggsave function to save the plot 
save_plot_with_ggsave(
  plot_obj = scree_plot,
  filename = "scree_plot.pdf",
  width = 12,
  height = 8
)



# 5 Differential gene expression analysis ==================================================================================================================================================================
# pairwise analysis was performed (1 vs 1)

# check coldata column name 
names(colData(mydata_bulk_RNA2))
names(rowData(mydata_bulk_RNA2))


# 5.1 Data wrangling ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Rename the 'Cell.population' column to 'cell_types'
colnames(colData(mydata_bulk_RNA2))[colnames(colData(mydata_bulk_RNA2)) == "Cell.population"] <- "cell_types"

# convert cell type column to factor
mydata_bulk_RNA2$cell_types <- factor(mydata_bulk_RNA2$cell_types)

# check cell types levels
levels(colData(mydata_bulk_RNA2)$cell_types)

# design and build DESeqDataSet object 
dds <- DESeqDataSet(mydata_bulk_RNA2, 
                             design = ~ Patient + cell_types)


# check experimental design 
design(dds)


dim(dds)
# drop PBMCs from the downstream analysis 
dds <- dds[, colData(dds)$cell_types != "PBMC"]

# drop the unused levels if it still there 
colData(dds) <- droplevels(colData(dds))
levels(colData(dds)$cell_types)


# filter out total genes count less 10 per row
dds <- dds[rowSums(counts(dds)) >= 10, ]


# remove the NAs from the data 
# extract the rownames (genes) from dds_pairwise
ensembl_ids <- rownames(dds)

# clean IDs ( remove any addition to name, strip suffix off)
ensembl_ids_clean <- sub("\\.\\d+$", "", ensembl_ids)

# query biomaRt for annotation
# connect to the database (bioMart) via the R package (biomaRt)

mart_obj <- useEnsembl(biomart = "genes",
                       dataset = "hsapiens_gene_ensembl"#,
                       #version = 109  
                       )

# get biomart saved in mart_obj and pass it to Ensembl database
# retrive the matching symbols and genes 
# attributes refers to the columns or Ensembl identifiers to match
annot <- getBM(attributes = c("ensembl_gene_id", 
                              "hgnc_symbol", 
                              "gene_biotype"),
               filters = "ensembl_gene_id",
               values = ensembl_ids_clean,
               mart = mart_obj)


# process annotation
# look for only protein_coding genes  and no duplicates 
annot_clean <- annot %>%
  filter(gene_biotype == "protein_coding",
         hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)


# create named vectors for look up
symbol_lookup <- setNames(annot_clean$hgnc_symbol, 
                          annot_clean$ensembl_gene_id)

biotype_lookup <- setNames(annot_clean$gene_biotype, 
                           annot_clean$ensembl_gene_id)

# create keep vector that stored the match symbols 
matched_genes <- ensembl_ids_clean %in% names(symbol_lookup)

# subset the latest inform from dds object previously built 
dds_obj <- dds[matched_genes, ]

# add annotations
rowData(dds_obj)$hgnc_symbol <- symbol_lookup[ensembl_ids_clean[matched_genes]]

rowData(dds_obj)$gene_biotype <- biotype_lookup[ensembl_ids_clean[matched_genes]]
rowData(dds_obj)

nrow(dds_obj)

# sample information after filtering 
dds_obj

### 5.2 further quality control assessment before DESeq2 ============================================================================================================================== 
#library  size check 


library_df <- data.frame(Sample = colnames(dds_obj),
                         Library_Size = colSums(counts(dds_obj))/1e6,
                         Condition = dds_obj$sampleName
                         )

head(library_df)

# check for outliers  
summary(library_df$Library_Size)


# check a box plot to identify outliers 

ggplot(library_df, aes(y = Library_Size)) +
  geom_boxplot() +
  theme_classic()




# calculate the outliers 
Q1 <- quantile(library_df$Library_Size, 0.25)
Q3 <- quantile(library_df$Library_Size, 0.75)

IQR_value <- Q3 - Q1
IQR_value

lower_limit <- Q1 - 1.5 * IQR_value
upper_limit <- Q3 + 1.5 * IQR_value

lower_limit
upper_limit

# check the which samples exhibit small and large library sizes.  

library_df %>%
  filter(
    Library_Size < lower_limit |
      Library_Size > upper_limit
  )



## create a general template for theme for all the QC plots 
publication_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )


# make changes to publication theme for fit best for the library size plot 

publication_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      size = 10,
      angle = 90,
      hjust = 1,
      face = "plain"
    ),
    plot.title = element_text(
      hjust = 0.5,
      size = 35,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 25,
      face = "bold")
  )

# visualise the library size 

library_size_plot <- ggplot(
  library_df,
  aes(
    x = Sample,
    y = Library_Size
  )
) +
  geom_col(width = 0.7) +
  labs(
    title = "Library size distribution",
    y = "Library size (million reads)",
    x = ""
  ) +
  scale_y_continuous(
    limits = c(0, max(library_df$Library_Size) * 1.1),
    expand = c(0, 0)
  ) +
  publication_theme



library_size_plot




# save the library size plot by calling the save function 


save_plot_with_ggsave(
  plot_obj = library_size_plot,
  filename = "library_size_plot.pdf",
  width = 30,
  height = 12
)


# variance-stabilizing transform for QC

vsd <- vst(dds_obj, blind = TRUE)

# perform the PCA 

pca_data <- plotPCA(
  vsd,
  intgroup = "cell_types",
  returnData = TRUE
)

nrow(vsd)
head(pca_data)

#check number of sample per group again to draw an ellipse 
table(pca_data$cell_types)

# convert PC fractions to percentage 
percent_var <- round(
  100 * attr(pca_data, "percentVar")
)


# generate a vst PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = cell_types,
    label = name
  )
) +
  geom_point(size = 3) +
  geom_text(
    size = 3,
    vjust = -1
  ) +
  stat_ellipse(
    aes(group = cell_types),
    type = "norm",
    linewidth = 1
  ) +
  labs(
    title = "PCA of RNA-seq samples",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    colour = "Cell Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # PCA title
    plot.title = element_text(
      size = 30,
      face = "bold",
      hjust = 0.5
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 25,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 25,
      face = "bold"
    ),
    
    # Axis text
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    
    # Legend inside plot (right side)
    legend.position = "inside",
    legend.position.inside = c(0.95, 0.2),
    legend.justification = c("right", "center"),
    
    # Legend formatting
    legend.title = element_text(
      size = 15,
      face = "bold"
    ),
    legend.text = element_text(
      size = 13
    ),
    legend.background = element_rect(
      fill = "white",
      colour = "black"
    )
  )




pca_plot

# call the save the pca plot 
save_plot_with_ggsave(
  plot_obj = pca_plot,
  filename = "deseq2_pca_plot.pdf",
  width = 20,
  height = 20
)






# check gene expression distribution using density plot 

vst_df <- as.data.frame(assay(vsd)) %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression"
  )

head(vst_df)


## add the cell_type column 
vst_df <- vst_df |>
  left_join(
    data.frame(
      sample = colnames(vsd),
      cell_types = colData(vsd)$cell_types
    ),
    by = "sample"
  )

head(vst_df)


# create ggplot of gene expression distribution 

celltype_colors <- c(
  cDC1 = "#377EB8",
  cDC2 = "#4DAF4A",
  DC3 = "#984EA3",
  B  = "#FF7F00",
  mo1 = "#E41A1C",
  mo2 = "#A65628",
  moint = "#F781BF",
  pDC = "#999999",
  tDC = "#66C2A5"
)


gene_expression_dist <- ggplot(
  vst_df,
  aes(
    x = expression,
    colour = cell_types,
    group = sample
  )
) +
  geom_density(
    linewidth = 1.2,
    alpha = 0.9
  ) +
  scale_colour_manual(
    values = celltype_colors
  ) +
  labs(
    title = "Distribution of VST-transformed gene expression",
    x = "Expression level (VST)",
    y = "Density",
    colour = "Cell type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      size = 30,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title.x = element_text(
      size = 25,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 25,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 20
    ),
    
    axis.text.y = element_text(
      size = 20
    ),
    
    legend.title = element_text(
      size = 25,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 20
    ),
    
    legend.key.height = unit(0.8, "cm")
  )

gene_expression_dist



save_plot_with_ggsave(
  plot_obj = gene_expression_dist,
  filename = "gene_expression_dist.pdf",
  width = 30,
  height = 25
)

# add facet wrap
gene_expression_dist_facet <- ggplot(
  vst_df,
  aes(
    x = expression,
    fill = cell_types,
    group = sample
  )
) +
  geom_density(
    alpha = 0.5,
    linewidth = 2
  ) +
  facet_wrap(
    ~ cell_types,
    ncol = 3,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = celltype_colors
  ) +
  labs(
    title = "Expression distribution by cell type",
    x = "VST expression",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    # Facet labels
    strip.text = element_text(
      size = 25,
      face = "bold"
    ),
    
    strip.background = element_rect(
      linewidth = 2
    ),
    
    # Remove legend
    legend.position = "none",
    
    # Plot title
    plot.title = element_text(
      size = 30,
      face = "bold",
      hjust = 0.5
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 22,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 22,
      face = "bold"
    ),
    
    # Axis tick labels
    axis.text.x = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 16,
      face = "bold"
    ),
    
    # Axis tick marks
    axis.ticks = element_line(
      linewidth = 1.2
    ),
    
    # Length of ticks
    axis.ticks.length = unit(
      0.25,
      "cm"
    )
  )

gene_expression_dist_facet

save_plot_with_ggsave(
  plot_obj = gene_expression_dist_facet,
  filename = "gene_expression_dist_facet.pdf",
  width = 30,
  height = 25
)






# Calculate sample distances
sampleDist <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDist)

# Sample annotation
annotation_df <- data.frame(
  Cell_type = dds_obj$cell_types
)

rownames(annotation_df) <- colnames(dds_obj)

head (annotation_df)

# Define annotation colours
celltype_colors <- c(
  B = "#E41A1C",
  cDC1 = "#377EB8",
  cDC2 = "#984EA3",
  DC3 = "#FF7F00",
  mo1 = "#4DAF4A",
  mo2 = "#A65628",
  moint = "#F781BF",
  pDC = "#999999",
  tDC = "#66C2A5"
)


# Annotation with customised legend
ha <- HeatmapAnnotation(df = annotation_df,
                        col = list(Cell_type = celltype_colors
                                   ),
                        annotation_name_gp = gpar(fontsize = 15,
                                                  fontface = "bold"
                                                  ),
                        annotation_legend_param = list(Cell_type = list(title = "Cell Type",
                                                                        title_gp = gpar(fontsize = 20,
                                                                                        fontface = "bold"
                                                                                        ),
                                                                        labels_gp = gpar(fontsize = 15
                                                                                         ),
                                                                        grid_height = unit(0.6,
                                                                                           "cm"
                                                                                           ),
                                                                        grid_width = unit(
                                                                          0.6,
                                                                          "cm"
                                                                          )
                                                                        )
                                                       )
                        )


# Create heatmap
sample_distance_heatmap <- Heatmap(sampleDistMatrix,
                                   name = "Distance",
                                   top_annotation = ha,
                                   column_title = "Sample-to-sample distance",
                                   column_title_gp = gpar(fontsize = 30,
                                                          fontface = "bold"
                                                          ),
                                   row_title = "Samples",
                                   row_title_gp = gpar(fontsize = 25,
                                                       fontface = "bold"
                                                       ),
                                   cluster_rows = TRUE,
                                   cluster_columns = TRUE,
                                   row_names_gp = gpar(
                                     fontsize = 8
                                     ),
                                   column_names_gp = gpar(fontsize = 8
                                                          ),
                                   
                                   heatmap_legend_param = list(title = "Distance",
                                                               
                                                               title_gp = gpar(fontsize = 20,
                                                                               fontface = "bold"
                                                                               ),
                                                               
                                                               labels_gp = gpar(fontsize = 15
                                                                                ),
                                                               
                                                               grid_height = unit(0.6,
                                                                                  "cm"
                                                                                  ),
                                                               
                                                               grid_width = unit(0.6,
                                                                                 "cm"
                                                                                 )
                                                               )
                                   )


# Draw heatmap
draw(sample_distance_heatmap,
     heatmap_legend_side = "right",
     annotation_legend_side = "right"
     )

pdf("plots/bulk_plot/sample_distance_heatmap.pdf",
    width = 20,
    height = 18
    )

#draw(sample_distance_heatmap,
     #heatmap_legend_side = "right",
     #annotation_legend_side = "right"
     #)

# close all active graphic devices 
while (!is.null(dev.list())) { dev.off() }


# 5.2  run DESeq2  =============================================================================================================================
dds_pairwise <- DESeq(dds_obj)


## post DGE visualisation=======================================================================================================================

# diagnostic dispersion plot 

# Adjust general axis text size and orientation before plotting
# 1. Open the file device FIRST (Set filename and resolution)


# Replace "plots" with the exact name of your existing folder
pdf("plots/bulk_plot/DESeq2_dispersion_plot.pdf", width = 21, height = 18)

par(las = 1, mar = c(5, 5, 4, 2) + 0.1)

plotDispEsts(
  dds_pairwise,
  main = "DESeq2 dispersion estimates",
  xlab = "Mean normalized counts",
  ylab = "Dispersion",
  cex.lab = 1.5,
  cex.axis = 1.0,
  cex.main = 2.0
)

while (!is.null(dev.list())) { dev.off() }




# Results and visualisation ======================================================================================================================


#check comparisons (coefficients)
resultsNames(dds_pairwise)

# retrieve the specific levels results

# positive value means high expression for instance in cDC1 and negative means high expression in B cells

colData(dds_pairwise)$cell_types

# retrieve cDC1 vs B as an example, then convert this into function to get other comparisons 

res_B_cDC1 <- results(dds_pairwise, contrast = c("cell_types", "B", "cDC1"))

res_B_cDC1
# add symbols 
res_B_cDC1$symbol <- rowData(dds_pairwise)[ rownames(res_B_cDC1), "symbol" ]
res_B_cDC1$genes <- rownames(res_B_cDC1)

# Now convert to df 
res_B_cDC1_df <- as.data.frame(res_B_cDC1)
#remove NAs
res_B_cDC1_df <- res_B_cDC1_df[!is.na(res_B_cDC1_df$padj), ]

head(res_B_cDC1_df)


# create a ggMAplot for B cells against cDC1 cells 
b_cdc1_ma_plot <- ggplot(res_B_cDC1_df, aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 2) +
  scale_color_manual(
    values       = c("FALSE" = "gray60", "TRUE" = "red"),
    labels       = c("Not significant", "Significant"),
    na.translate = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(
    ylim = c(-10, 10),
    xlim = range(log10(res_B_cDC1_df$baseMean), na.rm = TRUE)
  ) +
  labs(
    x     = "log10(mean normalized counts)",
    y     = "log2 Fold Change",
    title = "MA Plot (cDC1 vs B)",
    color = "Significance"
  ) +
  theme_bw(base_size = 14) +
  theme(
    # Centred plot title 
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    
    # Axis titles 
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    
    # Axis tick labels 
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    
    # Legend title + labels 
    legend.title    = element_text(size = 15, face = "bold"),
    legend.text     = element_text(size = 12, face = "bold"),
    
    # Legend positioned bottom-right inside plot
    legend.position      = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.background    = element_rect(fill = "white", color = "gray80")
  )
b_cdc1_ma_plot 



# volcano plot data processing  
# prepare the data for a volcano plot 
# colour the direction of DGE 

res_B_cDC1_noNA_df <- res_B_cDC1_df %>% 
  
  # create new columns named as direction and Col
  mutate(
    direction = case_when(
      
      # significant with positive fold change above threshold
      padj <= 0.05 & log2FoldChange >  1.5 ~ "high",
      
      # significant with negative fold change below threshold
      padj <= 0.05 & log2FoldChange < -1.5 ~ "low",
      
      # everything else is not significant
      TRUE ~ "NotSig"
    ))


top25cDC1_B <- res_B_cDC1_noNA_df%>%
  # first, keep only significant genes
  filter(padj < 0.01) %>%   
  # then sort
  arrange(desc(abs(log2FoldChange)), padj) %>%  
  # then take the top 25
  head(25) %>%   
  # get only the symbol
  pull(symbol)                                               


# create plotting groups 
res_B_cDC1_noNA_df <- res_B_cDC1_noNA_df %>%
  mutate(
    plot_group = case_when(
      symbol %in% top25cDC1_B ~ "Top DEG",
      direction == "high" ~ "high",
      direction == "low" ~ "low",
      TRUE ~ "NotSig"
    )
  )



# generate the volcano plot 

b_cdc1_volcano_plot <- ggplot(res_B_cDC1_noNA_df, 
                              aes(log2FoldChange, 
                                  -log10(padj)
                                  )
                              ) +
  
  geom_point(aes(color = plot_group), 
             alpha = 0.75, 
             size = 3.5) +
  
  scale_color_manual(values = c(high = "red",
                                low = "blue",
                                NotSig = "gray",
                                "Top DEG" = "darkgreen"
                                )
                     ) +
  
  # add thresholds
  geom_vline(xintercept = c(-1.5, 1.5), 
             linetype = "dashed", 
             linewidth = 0.6) +
  
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", 
             linewidth = 0.6) +
  
  #add top gene labels 
  geom_text_repel(data = subset(res_B_cDC1_noNA_df, 
                                symbol %in% top25cDC1_B),
                  aes(label = symbol),
                  size = 7,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.3
  ) +
  # add titles 
  labs(title = "Volcano Plot: cDC1 vs B cells Differential Expression",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value",
       color = "Expression levels"
       ) +
  
  # add theme 
  theme_classic(base_size = 14) +
  
  theme(plot.title = element_text(size = 25,
                                  face = "bold",
                                  hjust = 0.5
                                  ),
  # axis titles
  axis.title.x = element_text(size = 20,
                              face = "bold"),
  axis.title.y = element_text(size = 20, 
                              face = "bold"),
  
  # axis tick labels
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
  
  # legend
  legend.title = element_text(size = 18, face = "bold"),
  legend.text = element_text(size = 15),
  
  # legend position
  legend.position = "right"
  )

b_cdc1_volcano_plot




## 5.2.3 create functions for similar analyses ----------------------------------------------------------------------------------------------------------------------------------------------------------
# convert codes from 5.2 into a functions 

# This function get result Deseq result object for each specific comparison 
get_DESeq2_results_df <- function(
    dds_obj,
    contrast,
    symbol_column = "symbol"
    ) {
  
  # Run DESeq2 contrast
  res <- results(dds_obj,
                 contrast = contrast
                 )
  
  # Add gene symbols from rowData
  res$symbol <- rowData(dds_obj)[
    rownames(res),
    symbol_column
  ]
  
  # Add gene IDs
  res$genes <- rownames(res)
  
  # Convert to dataframe
  res_df <- as.data.frame(res)
  
  # Remove genes with missing adjusted p-values
  res_df <- res_df[!is.na(res_df$padj), ]
  
  return(res_df)
}


#  This function generate a MA plot

plot_MA_DESeq2_fun <- function(res_df,
                           title = "MA Plot",
                           padj_cutoff = 0.05
                           ) {
  
  ggplot(res_df,
         aes(x = log10(baseMean),
             y = log2FoldChange
             )
         ) +
    geom_point(aes(color = padj < padj_cutoff),
               alpha = 0.6,
               size = 3
               ) +
    
    scale_color_manual(
      values = c("FALSE" = "gray60",
                 "TRUE" = "red"
                 ),
      labels = c("Not significant",
                 "Significant"
                 ),
      na.translate = FALSE
      ) +
    
    geom_hline(yintercept = 0,
               linetype = "dashed"
               ) +
    
    coord_cartesian(
      ylim = c(-12, 12),
      xlim = range(log10(res_df$baseMean), 
                   na.rm = TRUE)
      ) +
    
    labs(
      x = "log10(mean normalized counts)",
      y = "log2 Fold Change",
      title = title,
      color = "Significance"
    ) +
    
    theme_bw(base_size = 14) +
    
    theme(
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = 0.5
        ),
      
      axis.title.x = element_text(size = 20,
                                  face = "bold"
                                  ),
      
      axis.title.y = element_text(size = 20,
                                  face = "bold"
                                  ),
      
      axis.text.x = element_text(size = 16,
                                 face = "bold"
                                 ),
      
      axis.text.y = element_text(size = 16,
                                 face = "bold"
                                 ),
      
      legend.title = element_text(size = 20,
                                  face = "bold"
                                  ),
      
      legend.text = element_text( size = 18,
                                   face = "bold"
                                  ),
      
      
      legend.position = c(0.95, 0.05),
      
      legend.justification = c(1, 0),
      
      legend.background = element_rect(
        fill = "white",
        color = "gray80"
      )
    )
}


# This function prepares the volcano plot data 

prepare_volcano_data_fun <- function(res_df, 
                                 padj_cutoff = 0.05,
                                 lfc_cutoff = 1.5,
                                 top_n = 25,
                                 top_padj_cutoff = 0.01,
                                 gene_column = "symbol"
                                 ) {
  
  # add DEG direction
  res_volcano_df <- res_df %>%
    mutate(
      direction = case_when(
        
        padj <= padj_cutoff & log2FoldChange > lfc_cutoff ~ "high",
        
        padj <= padj_cutoff & log2FoldChange < -lfc_cutoff ~ "low",
        
        TRUE ~ "NotSig"
      )
    )
  
   
  # Select top significant DEGs based on absolute fold change
  top_DEGs <- res_volcano_df %>%
    filter(padj < top_padj_cutoff) %>%
    arrange(desc(abs(log2FoldChange)), padj) %>%
    head(top_n) %>%
    pull(all_of(gene_column))
  
  
  # create final plotting groups
  res_volcano_df <- res_volcano_df %>%
    mutate(
      plot_group = case_when(
        .data[[gene_column]] %in% top_DEGs ~ "Top DEG",
        direction == "high" ~ "high",
        direction == "low" ~ "low",
        TRUE ~ "NotSig"
      )
    )
  
  
  return(res_volcano_df)
}


# This functions  displays volcano plot 

plot_volcano_DESeq2 <- function(
    volcano_df,
    title = "Volcano Plot: Differential Expression",
    lfc_cutoff = 1.5,
    padj_cutoff = 0.05,
    gene_column = "symbol"
    ) {
  
  ggplot(volcano_df,
         aes(x = log2FoldChange,
             y = -log10(padj)
             )
         ) +
    
    geom_point(aes(color = plot_group),
               alpha = 0.75,
               size = 3.5
               ) +
    
    scale_color_manual(
      values = c(high = "red",
                 low = "blue",
                 NotSig = "gray",
                 "Top DEG" = "darkgreen"
                 )
      ) +
    
    # add fold-change thresholds
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed",
               linewidth = 0.6
               ) +
    
    # Add significance threshold
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed",
               linewidth = 0.6
               ) +
    
    # Label top DEGs
    geom_text_repel(
      data = subset(volcano_df, plot_group == "Top DEG"),
      aes(label = .data[[gene_column]]),
      size = 7,
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3
    ) +
    
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Expression levels"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      plot.title = element_text(
        size = 25,
        face = "bold",
        hjust = 0.5
      ),
      
      axis.title.x = element_text(
        size = 20,
        face = "bold"
      ),
      
      axis.title.y = element_text(
        size = 20,
        face = "bold"
      ),
      
      axis.text.x = element_text(
        size = 15
      ),
      
      axis.text.y = element_text(
        size = 15
      ),
      
      legend.title = element_text(
        size = 18,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 15
      ),
      
      legend.position = "right"
    )
}





## 5.2.4  call the DESeq functions in 5.2.1 -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# call the DESeq2 comparison and result retrieval function 

res_B_cDC1_df <- get_DESeq2_results_df(dds_obj = dds_pairwise,
                                       contrast = c("cell_types", "B", "cDC1")
                                       )
head(res_B_cDC1_df)


# call the MA plot display function 
b_cdc1_ma_plot <- plot_MA_DESeq2_fun(res_df = res_B_cDC1_df,
                                     title = "MA Plot (B vs cDC1)"
                                     )

b_cdc1_ma_plot


# call the volcano data preparation function 
res_B_cDC1_volcano_df <- prepare_volcano_data_fun(res_df = res_B_cDC1_df,
                                                  gene_column = "symbol"
                                                  )
head(res_B_cDC1_volcano_df)


# call the volcano plot generation function 
b_cdc1_volcano_plot <- plot_volcano_DESeq2(volcano_df = res_B_cDC1_volcano_df,
                                           title = "Volcano Plot: B vs cDC1 Differential Expression"
                                           )

b_cdc1_volcano_plot


# call the save_plot_with_ggsave function 

# save the MA plot 
save_plot_with_ggsave(plot_obj = b_cdc1_ma_plot,
                      filename = "b_cdc1_ma_plot.pdf",
                      width = 15,
                      height = 15
                      )

# save volcano plot 
save_plot_with_ggsave(plot_obj = b_cdc1_volcano_plot,
                      filename = "b_cdc1_volcano_plot.pdf",
                      width = 20,
                      height = 20
                      )







# 5.3 Downstream analysis ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## create functions for downstream analyses 

#This function retrieve the all comparisons 
# The function also returned only the significant genes ( specified padj value and log2fold change threshold)
get_pairwise_DEGs <- function(dds_object, 
                              cell_type_col = "cell_types", 
                              log2FC_threshold = 1,
                              padj_threshold = 0.05) {
  
 
  #library(dplyr)
  #library(tibble)
  
  # get all unique cell types
  cts <- levels(colData(dds_object)[[cell_type_col]])
  
  # generate all pairwise combinations
  pairs <- combn(cts, 2, simplify = FALSE)
  
  # prepare storage lists 
  pairwise_results <- list()
  pairwise_summary <- list()
  
  # loop through pairs
  for (pair in pairs) {
    
    cell1 <- pair[1]
    cell2 <- pair[2]
    comparison_name <- paste0(cell1, "_vs_", cell2)
    
    message("Processing ", comparison_name)
    
    # DESeq2 results
    #first cell type  > cell1, second cell type > cell2 
    #DEG > cell1 vs cell2 
    #extract the results for the design cell types DEG
    
    res <- results(dds_object, contrast = c(cell_type_col, cell1, cell2))
    
    # convert to dataframe  and remove the NAs
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      filter(!is.na(log2FoldChange), !is.na(padj))
    
    # Add gene symbols if  not available
    if ("symbol" %in% colnames(rowData(dds_object))) {
      res_df$symbol <- rowData(dds_object)[res_df$gene_id, "symbol"]
    }
    
    # apply DEG filtering both cells 
    
    #first cell (cell1)
    sig_cell1 <- res_df %>%
      filter(log2FoldChange > log2FC_threshold,
             padj < padj_threshold) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      mutate(direction = paste0("higher_in_", cell1))
    
    # second cell (cell2)
    sig_cell2 <- res_df %>%
      filter(log2FoldChange < -log2FC_threshold,
             padj < padj_threshold) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      mutate(direction = paste0("higher_in_", cell2))
    
    # Combine extracted results from both cells 
    combined_results <- bind_rows(sig_cell1, sig_cell2)
    
    # Store all results in a list 
    pairwise_results[[comparison_name]] <- list(
      higher_in_cell1 = sig_cell1,
      higher_in_cell2 = sig_cell2,
      combined_results = combined_results,
      full_results = res_df            
    )
    
    # Summary
    pairwise_summary[[comparison_name]] <- data.frame(
      comparison = comparison_name,
      n_higher_in_cell1 = nrow(sig_cell1),
      n_higher_in_cell2 = nrow(sig_cell2),
      total_DEGs = nrow(combined_results)
    )
  }
  
  # Combine summary
  summary_df <- bind_rows(pairwise_summary)
  
  #  return pairwise list and summary_df 
  return(list(
    results = pairwise_results,
    summary = summary_df
  ))
}




# create a function retrieve cell of interest 
# This function  return data for a specific cell type of which genes are higher in all comparisons 
get_higher_in_cell <- function(pairwise_output, cell_type = "B") {
  
  
  # ensure the input is a list (required structure)
  if (!is.list(pairwise_output)) {
    stop("pairwise_output must be a list")
  }
  
  
  # Some outputs store results inside $results
  #check if the list is not empty 
  if (is.list(pairwise_output$results) && length(pairwise_output$results) > 0) {
    pairwise_output <- pairwise_output$results
  }
  
  
  # intitalise a list for storage 
  
  higher_list <- list()
  
  
  
  #  perform a loop to go through the comparison  and save through as "x"
  for (comp_name in names(pairwise_output)) {
    
    x <- pairwise_output[[comp_name]]
    
    # get comparison name
    # check the comparison name is available 
    
    comp <- if (!is.null(x$comparison_name) && length(x$comparison_name) == 1) {
      x$comparison_name
    } else {
      comp_name
    }
    
    # Split comparison into two cell types
    
    parts <- strsplit(comp, "_vs_")[[1]]
    
    # Skip invalid format
    if (length(parts) != 2) next
    
    
    
    #  Select correct table based on target cell type
    
    df <- NULL
    
    if (parts[1] == cell_type) {
      df <- x$higher_in_cell1
      
    } else if (parts[2] == cell_type) {
      df <- x$higher_in_cell2
      
    } else {
      next
    }
    
    
    
    #  convert all log2FC  to positive     
    # Keep magnitude only (no direction)
    if (!is.null(df) && "log2FoldChange" %in% colnames(df)) {
      df$log2FoldChange <- abs(df$log2FoldChange)
    }
    
    
    
    #  store valid results only 
    
    if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
      higher_list[[comp_name]] <- df
    }
  }
  
  
  
  # no results handling 
  
  if (length(higher_list) == 0) {
    warning("No comparisons found for cell type: ", cell_type)
    return(data.frame())
  }
  
  
  
  # combine the outputs 
  
  combined_df <- dplyr::bind_rows(higher_list, .id = "comparison")
  
  return(combined_df)
}




# This function filters and return genes which most discriminatory based on specified arguments 

top_genes_fn <- function(df_hd_data, top_n = Inf, min_comparisons = 1) {
  df_hd_data %>%
    # Fix: Use consistent column names (assuming DESeq2 standard)
    group_by(gene_id) %>%
    filter(n_distinct(comparison) >= min_comparisons) %>%
    arrange(desc(log2FoldChange), padj) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    ungroup() %>%
    slice_max(order_by = log2FoldChange, n = top_n, with_ties = FALSE)
}







# This function process and return data that can be used to generate a heat map
create_heatmap_data <- function(dds_pairwise, 
                                top_genes_df, 
                                metadata_cols = "cell_types", 
                                use_vst = TRUE) {
  
  
  # check that metadata columns exist
  
  missing_cols <- setdiff(metadata_cols, colnames(colData(dds_pairwise)))
  if(length(missing_cols) > 0) {
    stop("These metadata columns are missing in colData(dds_pairwise): ", 
         paste(missing_cols, collapse=", "))
  }
  
  
  # check that genes exist in the assay
  
  missing_genes <- setdiff(top_genes_df$gene_id, rownames(dds_pairwise))
  if(length(missing_genes) > 0) {
    warning("Some genes in top_genes_df are not found in the assay: ", 
            paste(missing_genes, collapse=", "))
  }
  
  # Keep only genes that exist
  valid_genes <- intersect(top_genes_df$gene_id, rownames(dds_pairwise))
  if(length(valid_genes) == 0) {
    stop("No valid genes found in the assay. Cannot generate heatmap data.")
  }
  
  
  # transform counts
  
  if(use_vst) {
    message("Using VST transformation (faster for large datasets)")
    vsd <- vst(dds_pairwise, blind = FALSE)
  } else {
    message("Using rlog transformation (better for small datasets)")
    vsd <- rlog(dds_pairwise, blind = FALSE)
  }
  
  
  # extract expression matrix for top genes
  
  mat <- assay(vsd)[valid_genes, , drop = FALSE]
  
  # Scale by row (z-score)
  mat_scaled <- t(scale(t(mat)))
  
  
  #  prepare sample annotations
  
  annotation_df <- as.data.frame(colData(dds_pairwise)[, metadata_cols, drop = FALSE])
  
  # ensure rownames match colnames of matrix
  rownames(annotation_df) <- colnames(mat_scaled)
  
  
  # return all components in a list
  
  return(list(
    scaled_matrix = mat_scaled,      # Scaled expression matrix (z-score)
    annotation_df = annotation_df,   # Sample metadata
    top_genes = top_genes_df,        # Original top genes dataframe
    transformation = ifelse(use_vst, "vst", "rlog")  # Transformation method used
  ))
}



# Heatmap generator with unified expression scaling and fixed cell-type palette
generate_heatmap <- function(heatmap_data,
                             main_title = "Top 25 cDC1 cells markers",
                             gene_symbols = NULL,
                             z_min = -2,
                             z_max = 2,
                             legend_breaks = c(-2, -1, 0, 1, 2),
                             cell_colors = c(
                               "cDC1"   = "#E69F00",
                               "cDC2"   = "#D55E00",
                               "DC3"    = "#CC79A7",
                               "pDC"    = "#F0E442",
                               "tDC"    = "#999999",  # Gray
                               "mo1"    = "#0072B2",
                               "mo2"    = "#56B4E9",
                               "moint"  = "#009E73",
                               "B"      = "#E41A1C"   # Red
                             )) {
  
  # 1. Prepare annotation dataframe
  if (is.data.frame(heatmap_data$annotation_df$cell_types)) {
    annotation_df <- heatmap_data$annotation_df$cell_types
    colnames(annotation_df) <- "Cell Types"
  } else {
    annotation_df <- data.frame(`Cell Types` = heatmap_data$annotation_df$cell_types,
                                check.names = FALSE,
                                stringsAsFactors = FALSE)
  }
  
  # 2. Map fixed colors to column annotation
  celltype_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df,
    col = list(`Cell Types` = cell_colors),
    annotation_name_side = "right"
  )
  
  # 3. Set up row labels
  if (!is.null(gene_symbols)) {
    row_labels <- gene_symbols
  } else {
    row_labels <- rownames(heatmap_data$scaled_matrix)
  }
  
  # 4. Define fixed expression color scale (navy -> white -> firebrick3)
  color_fun <- circlize::colorRamp2(
    breaks = seq(z_min, z_max, length.out = 100),
    colors = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  
  
  
  
  
  # 5. Create Heatmap
  ComplexHeatmap::Heatmap(
    matrix = heatmap_data$scaled_matrix,
    top_annotation = celltype_annotation,
    row_labels = row_labels,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = 6),
    col = color_fun,
    column_title = paste(main_title, "-", heatmap_data$transformation),
    name = "expression",
    border = TRUE,
    heatmap_legend_param = list(
      title = "Z-score",
      at = legend_breaks,                   # Locks legend tick locations
      labels = as.character(legend_breaks), # Consistent tick text
      title_position = "topcenter",
      legend_direction = "horizontal",
      title_gp = grid::gpar(fontsize = 10, fontface = "bold")
    )
  )
}


  


# 5.3.1 call the downstream analysis functions ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# call the pairwise function (get_pairwise_DEGs) and wrangle the data 
# shows which comparison is be process at a time 
# it returns 2 outputs, resuslt and a summary table of up-and-down regulated genes for each comparsions  
pairwise_output <- get_pairwise_DEGs(
  dds_object = dds_pairwise,
  cell_type_col = "cell_types",
  log2FC_threshold = 0
  )

# check the output of the pairwise function 
names(pairwise_output)
pairwise_output$summary


# convert the summary from get_pairwise function into table 
summary_tbl <- pairwise_output$summary %>%
  arrange(comparison)


colnames(pairwise_output$summary)

# rename the table columns  
summary_tbl1 <- summary_tbl %>%
  rename(
    Comparison = comparison,
    `Up-regulated in Cell1` = n_higher_in_cell1,
    `Up-regulated in Cell2` = n_higher_in_cell2,
    `Total DEGs` = total_DEGs
  )

summary_tbl1

#create a coloured table 
# customise the table 
gt_tbl_pub <- summary_tbl1 %>%
  gt() %>%
  
  # Left-align all columns
  cols_align(
    align = "left",
    columns = everything()
  ) %>%
  
  # Shrink columns to fit content
  cols_width(
    everything() ~ "auto"
  ) %>%
  
  # Apply gradient colours to key columns only
  data_color(
    columns = "Up-regulated in Cell1",
    fn = col_numeric(
      palette = c("white", "#1f77b4"),
      domain = NULL
    )
  ) %>%
  
  data_color(
    columns = "Up-regulated in Cell2",
    fn = col_numeric(
      palette = c("white", "#ff7f0e"),
      domain = NULL
    )
  ) %>%
  
  data_color(
    columns = "Total DEGs",
    fn = col_numeric(
      palette = c("white", "steelblue"),
      domain = NULL
    )
  ) %>%
  
  # Bold Comparison column
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = "Comparison")
  ) %>%
  
  # Grey header + bold column names
  tab_style(
    style = list(
      cell_fill(color = "grey90"),
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels()
  ) %>%
  
  # Subtle borders
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(2),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table_body.hlines.color = "grey80"
  ) %>%
  
  # row striping
  opt_row_striping()

# display the table
gt_tbl_pub



# save the plot 
gtsave(
  gt_tbl_pub,
  "plots/bulk_plot/table_up_down_regulated_genes.pdf"
)

# call higher function to get higher (up-regulated) genes per cell 

# call get_higher_in_cell function 

names(pairwise_output)
names(pairwise_output$results)
names(pairwise_output$summary)

B_up_df <- get_higher_in_cell(pairwise_output, cell_type = "B")
cDC1_up_df <- get_higher_in_cell(pairwise_output, cell_type = "cDC1")
cDC2_up_df <- get_higher_in_cell(pairwise_output, cell_type = "cDC2")
DC3_up_df <- get_higher_in_cell(pairwise_output, cell_type = "DC3")
pDC_up_df <- get_higher_in_cell(pairwise_output, cell_type = "pDC")
tDC_up_df <- get_higher_in_cell(pairwise_output, cell_type = "tDC")
mo1_up_df <- get_higher_in_cell(pairwise_output, cell_type = "mo1")
mo2_up_df <- get_higher_in_cell(pairwise_output, cell_type = "mo2")
moint_up_df <- get_higher_in_cell(pairwise_output, cell_type = "moint")


# check the number of comparisons in which the this were identified from   
B_up_df %>%
  group_by(gene_id) %>%
  summarise(n_comparisons = n_distinct(comparison)) %>%
  arrange(desc(n_comparisons))


# call the function that would return most differential expressed genes 
top25_B_genes <- top_genes_fn(df_hd_data = B_up_df , top_n = 25, min_comparisons = 8) 
top25_cDC1_genes <- top_genes_fn(df_hd_data = cDC1_up_df , top_n = 25, min_comparisons = 8) 

top25_cDC2_genes <- top_genes_fn(df_hd_data = cDC2_up_df , top_n = 25, min_comparisons = 7) 
top25_DC3_genes <- top_genes_fn(df_hd_data = DC3_up_df , top_n = 25, min_comparisons = 7) 
top25_pDC_genes <- top_genes_fn(df_hd_data = pDC_up_df , top_n = 25, min_comparisons = 8) 
top25_tDC_genes <- top_genes_fn(df_hd_data = tDC_up_df , top_n = 25, min_comparisons = 8) 
top25_mo1_genes <- top_genes_fn(df_hd_data = mo1_up_df , top_n = 25, min_comparisons = 8) 
top25_mo2_genes <- top_genes_fn(df_hd_data = mo2_up_df , top_n = 25, min_comparisons = 7) 
top25_moint_genes <- top_genes_fn(df_hd_data = moint_up_df , top_n = 25, min_comparisons = 7)


# combine the identified top genes per cell into a one single object 

all_cell_gene_table <- bind_rows(
  tibble(Cell = "B",       Gene = top25_B_genes),
  tibble(Cell = "cDC1",    Gene = top25_cDC1_genes),
  tibble(Cell = "cDC2",    Gene = top25_cDC2_genes),
  tibble(Cell = "DC3",     Gene = top25_DC3_genes),
  tibble(Cell = "pDC",     Gene = top25_pDC_genes),
  tibble(Cell = "tDC",     Gene = top25_tDC_genes),
  tibble(Cell = "mo1",     Gene = top25_mo1_genes),
  tibble(Cell = "mo2",     Gene = top25_mo2_genes),
  tibble(Cell = "moint",   Gene = top25_moint_genes)
)


all_cell_gene_table 





cell_gene_table <- bind_rows(
  tibble(Cell = "B",       Gene = top25_B_genes$symbol),
  tibble(Cell = "cDC1",    Gene = top25_cDC1_genes$symbol),
  tibble(Cell = "cDC2",    Gene = top25_cDC2_genes$symbol),
  tibble(Cell = "DC3",     Gene = top25_DC3_genes$symbol),
  tibble(Cell = "pDC",     Gene = top25_pDC_genes$symbol),
  tibble(Cell = "tDC",     Gene = top25_tDC_genes$symbol),
  tibble(Cell = "mo1",     Gene = top25_mo1_genes$symbol),
  tibble(Cell = "mo2",     Gene = top25_mo2_genes$symbol),
  tibble(Cell = "moint",   Gene = top25_moint_genes$symbol)
)

cell_gene_table

# pivot wide the combine top genes data; 
cell_gene_horizontal <- cell_gene_table %>%
  group_by(Cell) %>%
  summarise(
    Genes = paste(Gene, collapse = ", "),
    .groups = "drop"
  )

cell_gene_horizontal



# format and customise the the pivot wide table 

marker_gene_gt <- cell_gene_horizontal %>%
  gt() %>%
  
  # title of the table 
  tab_header(title = "Top 25 Marker Genes Identified for Each Cell Type"
             ) %>%
  
  # rename columns
  cols_label(Cell = "Cell Type",
             Genes = "Top 25 Marker Genes"
             ) %>%
  
  # centre and bold the column labels 
  tab_style(style = cell_text(align = "center",
                              weight = "bold"
                              ),
    locations = cells_column_labels()
    ) %>%
  
  # bold row names
  tab_style(style = cell_text(weight = "bold",
                              size = px(18)),
            locations = cells_body(columns = Cell)
            ) %>%
  
  # column labels
  tab_style(style = list(cell_fill(color = "grey90"),
                         cell_text(weight = "bold",
                                   size = px(18))
                         ),
  locations = cells_column_labels()
  ) %>%

  # increasing body text 
  tab_style(style = cell_text(size = px(15)),
            locations = cells_body()
            ) %>%
  
  # table options
  tab_options(
    table.font.size = px(15),
    data_row.padding = px(12),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table_body.hlines.color = "grey80"
  ) %>%
  
  # alternate row shading
  opt_row_striping()

marker_gene_gt


# then save the table for the top markers 
gtsave(
  marker_gene_gt,
  filename = "plots/bulk_plot/top25_marker_genes_table.pdf"
)



# again create another list of the top genes from all cell types picking only the top 3 genes from each cell
# then make a visual table 
cell_dfs <- list("cDC1" = top25_cDC1_genes,
                 "cDC2" = top25_cDC2_genes,
                 "DC3"  = top25_DC3_genes,
                 "Bcell" = top25_B_genes,
                 "mo1"  = top25_mo1_genes,
                 "mo2"  = top25_mo2_genes,
                 "moint"  = top25_moint_genes,
                 "pDC" = top25_pDC_genes,
                 "tDC" = top25_tDC_genes
                 )




## shared genes 

shared_genes <- list(
  cDC1 = top25_cDC1_genes$symbol,
  cDC2 = top25_cDC2_genes$symbol,
  DC3  = top25_DC3_genes$symbol,
  Bcell = top25_B_genes$symbol,
  mo1 = top25_mo1_genes$symbol,
  mo2 = top25_mo2_genes$symbol,
  moint = top25_moint_genes$symbol,
  pDC = top25_pDC_genes$symbol,
  tDC = top25_tDC_genes$symbol
)


intersect(top25_DC3_genes$symbol, top25_cDC2_genes$symbol)


# generate UPsetR plot 

shared_genes_plot<- upset(
  fromList(shared_genes ),
  sets = c("cDC1", "cDC2",
           "DC3","Bcell",
           "mo1","mo2", 
           "moint", "pDC", 
           "tDC"),
  keep.order = TRUE,
  nintersects = 10,
  order.by = "freq",
  sets.bar.color = "grey30",
  main.bar.color = "steelblue",
  matrix.color = "black",
  text.scale = c(3,2.2,2.2,2.2,2.5,3),
  point.size = 5,
  line.size = 1.5,
  show.numbers = "yes",
  mainbar.y.label = "Number of shared genes",
  sets.x.label = "Number of genes per cell type",
  mb.ratio = c(0.65,0.35)
)

shared_genes_plot


grid.text(
  "Shared Differentially Expressed Genes",
  y = unit(0.98, "npc"),
  gp = gpar(fontsize = 30, fontface = "bold")
)


shared_genes_plot

pdf(
  "plots/bulk_plot/shared_genes_upset.pdf",
  width = 20,
  height = 20
  )




# pick the top 3 for each cell type and  make a table   
top3_stats_summary <- lapply(names(cell_dfs), function(name) {
  cell_dfs[[name]] %>%
    slice_head(n = 3) %>%
    mutate(cell_type = name) %>%
    
 # Selecting and renaming for clarity if necessary
    dplyr::select(cell_type, symbol, log2FoldChange, padj)
  }) %>%
  bind_rows()

top3_stats_summary


# format the table for display
formatted_top3 <- top3_stats_summary %>%
  rename(
    `Cell Type` = cell_type,
    `Gene Symbol` = symbol,
    `log2 Fold Change` = log2FoldChange,
    `Adjusted p-value` = padj
  ) %>%
  mutate(
    `log2 Fold Change` = round(`log2 Fold Change`, 2),
    `Adjusted p-value` = formatC(`Adjusted p-value`, format = "e", digits = 2)
  )

# view the result
print(formatted_top3)


# Create a customised gt table
# Create a customised gt table
top3_DEGs_per_celltype <- formatted_top3 %>%
  gt() %>%
  
  # Table title
  tab_header(
    title = "Top 3 Differentially Expressed Genes Across Myeloid and Lymphoid Subsets"
  ) %>%
  
  # Format table title
  tab_style(
    style = cell_text(
      weight = "bold",
      size = px(14),
      align = "center"
    ),
    locations = cells_title()
  ) %>%
  
  # Column alignment
  cols_align(
    align = "left",
    columns = everything()
  ) %>%
  
  # Bold first column (cell type)
  tab_style(
    style = cell_text(
      weight = "bold"
    ),
    locations = cells_body(columns = 1)
  ) %>%
  
  # Italicize gene names (second column)
  tab_style(
    style = cell_text(
      style = "italic"
    ),
    locations = cells_body(columns = 2)
  ) %>%
  
  # Column header styling
  tab_style(
    style = list(
      cell_fill(color = "#1F4E79"),
      cell_text(
        weight = "bold",
        color = "white",
        size = px(12),
        align = "center"
      )
    ),
    locations = cells_column_labels()
  ) %>%
  
  # Add biological row groups
  tab_row_group(
    label = "Dendritic Cells",
    rows = 1:9
  ) %>%
  
  tab_row_group(
    label = "Lymphoid Cells",
    rows = 10:12
  ) %>%
  
  tab_row_group(
    label = "Monocyte Subsets",
    rows = 13:21
  ) %>%
  
  # Style biological row labels
  tab_style(
    style = list(
      cell_fill(color = "lightblue"),
      cell_text(
        weight = "bold",
        size = px(12),
        color = "black"
      )
    ),
    locations = cells_row_groups()
  ) %>%
  
  # Table formatting options
  tab_options(
    table.font.size = px(12),
    data_row.padding = px(5),
    table.width = pct(60),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table_body.hlines.color = "grey80",
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold"
  ) %>%
  
  # Alternate row shading
  opt_row_striping()


# Display table
top3_DEGs_per_celltype

# save the table 
gtsave(
  top3_DEGs_per_celltype,
  filename = "plots/bulk_plot/top3_DEGs_per_celltype.pdf",)




# 5.5 Heat map visualisation ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# again pick the top 3 genes from each cell type and make global heatmap to show gene cell type association 
top3_combine_genes <- imap_dfr(cell_dfs, ~ { .x %>%
    
    # take first 3 rows
    slice_head(n = 3) %>%  
    
    # remember which df it was from
    mutate(source = .y)
  })                        

top3_combine_genes

# pick on unique gene_id to avoid duplicates check for duplicated gene 
top3_unique <- top3_combine_genes %>%
  distinct(gene_id, .keep_all = TRUE)

# calls the creat_heatmap_data function to prepare heatmap datasets 

# global heatmap data processing 
heatmap_data <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                    top_genes_df = top3_unique, 
                                    metadata_cols = "cell_types", 
                                    use_vst = TRUE)

# B cells heatmap data 
heatmap_data_B <-    create_heatmap_data(dds_pairwise = dds_pairwise, 
                                         top_genes_df = top25_B_genes, 
                                         metadata_cols = "cell_types", 
                                         use_vst = TRUE)

# cDC1 
heatmap_data_cDC1 <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                         top_genes_df = top25_cDC1_genes, 
                                         metadata_cols = "cell_types", 
                                         use_vst = TRUE)
# cDC2

heatmap_data_cDC2 <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                         top_genes_df = top25_cDC2_genes, 
                                         metadata_cols = "cell_types", 
                                         use_vst = TRUE)

# DC3
heatmap_data_DC3 <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                        top_genes_df = top25_DC3_genes, 
                                        metadata_cols = "cell_types", 
                                        use_vst = TRUE)

# pDC
heatmap_data_pDC <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                        top_genes_df = top25_pDC_genes, 
                                        metadata_cols = "cell_types", 
                                        use_vst = TRUE)

# tDC
heatmap_data_tDC <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                        top_genes_df = top25_tDC_genes, 
                                        metadata_cols = "cell_types",
                                        use_vst = TRUE)

# mo1
heatmap_data_mo1 <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                        top_genes_df = top25_mo1_genes, 
                                        metadata_cols = "cell_types", 
                                        use_vst = TRUE)
# mo2 
heatmap_data_mo2 <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                        top_genes_df = top25_mo2_genes, 
                                        metadata_cols = "cell_types", 
                                        use_vst = TRUE)

# mo2 
heatmap_data_moint <- create_heatmap_data(dds_pairwise = dds_pairwise, 
                                          top_genes_df = top25_moint_genes, 
                                          metadata_cols = "cell_types", 
                                          use_vst = TRUE)


# heatmaps of the different cell types plus a combine cell type top 3 genes 
# top3 combine genes (global) heatmap 
# global (combine cell type) heatmap 
heatmap_plot <- generate_heatmap(heatmap_data = heatmap_data,
                  main_title = "Top 25 combine markers",
                                gene_symbols = heatmap_data$top_genes$symbol
                                )
heatmap_plot

# B cells heat map 

b_cell_heatmap <- generate_heatmap(heatmap_data = heatmap_data_B,
                                 main_title = "Top B cell markers",
                                 gene_symbols = heatmap_data_B$top_genes$symbol
                                 )

b_cell_heatmap
# cDC1 heat map 
cDC1_heatmap <- generate_heatmap(heatmap_data = heatmap_data_cDC1,
                                 main_title = "Top cDC1 markers",
                                 gene_symbols = heatmap_data_cDC1$top_genes$symbol
                                 )
cDC1_heatmap
# cDC2 heat map 
cDC2_heatmap <- generate_heatmap(heatmap_data = heatmap_data_cDC2,
                                 main_title = "Top cDC2 markers",
                                 gene_symbols = heatmap_data_cDC2$top_genes$symbol
                                 )
# DC3 heat map 
DC3_heatmap <- generate_heatmap(heatmap_data = heatmap_data_DC3 ,
                                 main_title = "Top DC3 markers",
                                 gene_symbols = heatmap_data_DC3$top_genes$symbol
                                )
# pDC heat map 
pDC_heatmap <- generate_heatmap(heatmap_data = heatmap_data_pDC,
                                 main_title = "Top pDC markers",
                                 gene_symbols = heatmap_data_pDC$top_genes$symbol
                                )
# tDC heat map 
tDC_heatmap <- generate_heatmap(heatmap_data = heatmap_data_tDC,
                                 main_title = "Top tDC markers",
                                 gene_symbols = heatmap_data_tDC$top_genes$symbol
                                )
# mo1 heat map 
mo1_heatmap <- generate_heatmap(heatmap_data = heatmap_data_mo1,
                                 main_title = "Top mo1 markers",
                                 gene_symbols = heatmap_data_mo1$top_genes$symbol
                                )
# mo2 heat map 
mo2_heatmap <- generate_heatmap(heatmap_data = heatmap_data_mo2,
                                 main_title = "Top mo2 markers",
                                 gene_symbols = heatmap_data_mo2$top_genes$symbol
                                )
# moint heat map 
moint_heatmap <- generate_heatmap(heatmap_data = heatmap_data_moint,
                                 main_title = "Top moint markers",
                                 gene_symbols = heatmap_data_moint$top_genes$symbol
                                 )


# make a list of all the heatmaps 
heatmap_list <- list(
  global_heatmap = heatmap_plot,
  Bcell = b_cell_heatmap,
  cDC1 = cDC1_heatmap,
  cDC2 = cDC2_heatmap,
  DC3 = DC3_heatmap,
  pDC = pDC_heatmap,
  tDC = tDC_heatmap,
  mo1  = mo1_heatmap,
  mo2  = mo2_heatmap,
  moint = moint_heatmap
)

# create a heatmap saving function 
save_complex_heatmap <- function(heatmap_obj,
                                 filename,
                                 width = 8,
                                 height = 6) {
  pdf(
    filename,
    width = width,
    height = height
  )
  
  ComplexHeatmap::draw(heatmap_obj)
  
  dev.off()
}


# call the heat map saving funtion using a for loop
for (i in names(heatmap_list)) {
  
  save_complex_heatmap(heatmap_obj = heatmap_list[[i]],
    filename = paste0("plots/bulk_plot/", i, "_heatmap.pdf"),
    width = 8,
    height = 6
  )
}





