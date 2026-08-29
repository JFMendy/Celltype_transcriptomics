

# ============================================= Single cell Analysis =====================================================================================

# 0.1 load the libraries  and data -----------------------------------------------------------------------------------------------------------------------
# libraries 
{library(Seurat)
  library(presto)
  library(gt)
  library(dplyr)
  library(ggplot2)
  library(knitr)
  library(kableExtra)
  library(patchwork)
  library(SingleCellExperiment)
  library(scDblFinder)
}


  
## read in the scRNA-seq  data
single_cell_data <-readRDS("data/scRNASeq_annotated/hsobj_20221209_annotated.rds")

top_bulk_genes <-readRDS("results/cell_gene_table.rds")
top_bulk_genes_list <-readRDS("results/all_cell_gene_table.rds")

# 1 preprocess the data ==================================================================================================================================

# make copy of the data  
single_cell_data2 <- single_cell_data      



rm(single_cell_data)
# 2 scRNA-seq inspection ==============================================================================================================================================

# check the structure of the seurat object 
glimpse(single_cell_data2)  

# using head to return only the first 6 rows  
head(single_cell_data2@meta.data)          

# number of rows and columns 
dim(single_cell_data2)

# check 10 row  names 
head(rownames(single_cell_data2), 10)       



# 3     Quality control and cleaning ======================================================================================================================================

# clean the data by removing: dying cells -> exhibit high mitochondrial contamination 
# by removing genes starting with MT -> mitochondrial genes 
# store mitochondrial percentage in object meta data in "percent.mt" column 

single_cell_data2[["percent.mt"]] <- PercentageFeatureSet(single_cell_data2, 
                                                          pattern = "^MT-")

head(single_cell_data2@meta.data)
head(single_cell_data2$percent.mt)

# 3.1   Pre cleaning visualization ---------------------------------------------------------------------------------------------------------------------------
# Visualize QC Metrics  using the Violin plots
# plot custmoise 
# put all cells into one group 


qc_violin_before <- 
  VlnPlot(
    single_cell_data2,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    pt.size = 0.1,
    ncol = 3
  )

qc_violin_before 



ggsave(
  "plots/scrna_plot/QC_violin_before.pdf",
  plot = qc_violin_before,
  width = 10,
  height = 5,
  units = "in"
)

rm(qc_violin_before)
# Scatter plots used to observed the relationship btw genes or percent_MT vs counts 

scatter_counts_features_before <- FeatureScatter(
  single_cell_data2,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) +
  labs(
    title = "RNA Counts vs Detected Genes",
    x = "RNA Counts",
    y = "Detected Genes"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text = element_text(size = 13),
    legend.position = "none"
  )

scatter_counts_mt_before <- FeatureScatter(
  single_cell_data2,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
) +
  labs(
    title = "RNA Counts vs Mitochondrial RNA",
    x = "RNA Counts",
    y = "Mitochondrial RNA (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.text = element_text(size = 13),
    legend.position = "none"
  )

qc_scatter_before <-
  scatter_counts_features_before +
  scatter_counts_mt_before

qc_scatter_before



ggsave(
  "plots/scrna_plot/qc_scatter_before.pdf",
  plot = qc_scatter_before,
  width = 10,
  height = 15,
  units = "in"
)

rm(qc_scatter_before,scatter_counts_features_before, scatter_counts_mt_before)

# 3.2   post filtering visualisation -------------------------------------------------------------------------------------------------------------------------------------------------

# filter the genes and counts  
# remove low quality cells (dead or dying cells), cell-free RNA, or doublets)
 
# Based on these plots, define thresholds 


# Gene count thresholds 
min_genes <- 200
max_genes <- 6000       
max_mt    <- 20         
min_counts <- 500       

# cleaning the data 
single_cell_data2 <- subset(single_cell_data2,
                           subset = nFeature_RNA > min_genes &
                             nFeature_RNA < max_genes &
                             percent.mt < max_mt &
                             nCount_RNA > min_counts)

 

# cleaning the working space 

rm(min_counts, min_genes, max_genes,max_mt)





# prevent random change in the result 
set.seed(1234)

#convert a seurat object into singleCellexperiement object 
sce <- as.SingleCellExperiment(single_cell_data2)


# this create singlets and doublets and doubletsscores columnns
sce <- scDblFinder(sce,
                   samples = sce$Sample
                   )


# add scDblFinder results back to the Seurat object
single_cell_data2$scDblFinder.score <- sce$scDblFinder.score
single_cell_data2$scDblFinder.class <- sce$scDblFinder.class

colnames(single_cell_data2@meta.data)

table(single_cell_data2$scDblFinder.class)

# Save the threshold as object-level information
single_cell_data2@misc$scDblFinder.threshold <- metadata(sce)$scDblFinder.threshold

summary(single_cell_data2$scDblFinder.score)


single_cell_data2@misc$scDblFinder.threshold


table(single_cell_data2$scDblFinder.class)

mean(single_cell_data2$scDblFinder.class == "doublet") * 100


# generate histogram of doublets 

doublet_score_plot <- ggplot(single_cell_data2@meta.data,
                             aes(x = scDblFinder.score,
                                 fill = scDblFinder.class
                                 )
                             ) +
  
  geom_histogram(bins = 50,
                 alpha = 0.7,
                 position = "identity"
                 ) +
  
  # scDblFinder threshold
  geom_vline(xintercept = single_cell_data2@misc$scDblFinder.threshold,
             linetype = "dashed",
             linewidth = 1.2
             ) +
  
  # threshold value
  annotate("text",
           x = single_cell_data2@misc$scDblFinder.threshold,
           y = Inf,
           label = paste0("Threshold = ",
                          round(single_cell_data2@misc$scDblFinder.threshold, 4)
                          ),
           vjust = 1.5,
           hjust = -0.05,
           size = 6
           ) +
  
  scale_fill_manual(values = c(singlet = "grey70",
                               doublet = "#D73027"
                               )
                    ) +
  
  
  labs(title = "Distribution of scDblFinder scores",
       x = "scDblFinder score",
       y = "Number of cells",
       fill = "Cell classification"
       ) +
  
  scale_x_continuous(limits = c(0, 1.0),
                     breaks = seq(0, 1.0, by = 0.1),
                     expand = c(0.01, 0)
                     ) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))
                     ) +
  
  theme_classic(base_size = 14) +
  
  theme(plot.title = element_text(size = 25,
                                  face = "bold",
                                  hjust = 0.5
                                  ),
  
  axis.title.x = element_text(size = 20,
                              face = "bold"
                              ),
  axis.title.y = element_text(size = 20,
                              face = "bold"
                              ),

  axis.text.x = element_text(size = 15
                             ),
  
  axis.text.y = element_text(
    size = 15
    ),
  
  axis.ticks = element_line(linewidth = 0.7
                            ),
  
  axis.ticks.length = unit(0.2,
                           "cm"
                           ),
  
  legend.title = element_text(size = 15,
                              face = "bold"
                              ),
  
  legend.text = element_text(size = 15
                             )
  )

doublet_score_plot



ggsave(
  "plots/scrna_plot/doublet_score_plot.pdf",
  plot = doublet_score_plot,
  width = 8,
  height = 10,
  units = "in"
)




# generate doublets summary table 


doublet_summary <- tibble(Category = c(rep("Score distribution", 7),
                                       rep("Classification", 5)
                                       ),
Metric = c("Cells assessed",
           "Minimum",
           "1st quartile",
           "Median",
           "Mean",
           "3rd quartile",
           "Maximum",
           "Classification threshold",
           "Singlets",
           "Singlets (%)",
           "Doublets",
           "Doublets (%)"
           ),

Value = c(ncol(single_cell_data2),
          min(single_cell_data2$scDblFinder.score, na.rm = TRUE),
          quantile(single_cell_data2$scDblFinder.score, 0.25, na.rm = TRUE),
          median(single_cell_data2$scDblFinder.score, na.rm = TRUE),
          mean(single_cell_data2$scDblFinder.score, na.rm = TRUE),
          quantile(single_cell_data2$scDblFinder.score, 0.75, na.rm = TRUE),
          max(single_cell_data2$scDblFinder.score, na.rm = TRUE),
          single_cell_data2@misc$scDblFinder.threshold,
          sum(single_cell_data2$scDblFinder.class == "singlet"),
          mean(single_cell_data2$scDblFinder.class == "singlet") * 100,
          sum(single_cell_data2$scDblFinder.class == "doublet"),
          mean(single_cell_data2$scDblFinder.class == "doublet") * 1
          )
)

doublet_summary



doublet_tbl <- doublet_summary %>%
  gt() %>%
  fmt_number(columns = Value,
             rows = Metric %in% c("Minimum",
                                  "1st quartile",
                                  "Median",
                                  "Mean",
                                  "3rd quartile",
                                  "Maximum",
                                  "Classification threshold"
                                  ),
             decimals = 2
             ) %>%
  
  fmt_number(columns = Value,
             rows = Metric %in% c(
               "Singlets (%)",
               "Doublets (%)"
               ),
             
             decimals = 2,
             suffix = "%"
             ) %>%
  
  fmt_number(columns = Value,
             rows = Metric %in% c("Cells assessed",
                                  "Singlets",
                                  "Doublets"
                                  ),
             decimals = 0,
             use_seps = TRUE
             ) %>%
  
  tab_header(title = md(""),
             subtitle = md("**Single-cell double assessment**")
             ) %>%
  
  cols_label(Category = "Category",
             Metric = "Measure",
             Value = "Result"
             ) %>%
  
  tab_style(style = cell_text(size = px(16),
                              weight = "bold"
                              ),
            locations = cells_column_labels()
            ) %>%
  
  tab_style(style = cell_text(size = px(12), weight = "bold"
                              ),
    locations = cells_body()
    ) %>%
  
  tab_options(table.font.size = px(1),
              heading.title.font.size = px(20),
              heading.subtitle.font.size = px(25),
              data_row.padding = px(14),
              column_labels.padding = px(12),
              table.width = pct(80)
              ) %>%
  cols_width(
    Category ~ px(120),
    Metric ~ px(90),
    Value ~ px(70)
  ) %>%
  tab_options(
    table.width = pct(40)
  )
 
doublet_tbl


gtsave(
  doublet_tbl,
  "plots/scrna_plot/doublet_tbl.pdf",
  vwidth = 800,
  vheight = 1000
  )


# cleaning the data 

#subset and create clean subseurat object



single_cell_singlets <- subset(single_cell_data2,
                               subset = scDblFinder.class == "singlet"
                               )




# violin plot 

qc_violin_after <- VlnPlot(single_cell_singlets ,
                           features = c("nFeature_RNA",
                                        "nCount_RNA",
                                        "percent.mt"
                                        ),
                           pt.size = 0.1,
                           ncol = 3
                           )

qc_violin_after 



ggsave(
  "plots/scrna_plot/qc_violin_after .pdf",
  plot = qc_violin_after ,
  width = 10,
  height = 5,
  units = "in"
)



# scatter plots 
scatter_counts_features_after <- FeatureScatter(single_cell_singlets,
                                                feature1 = "nCount_RNA",
                                                feature2 = "nFeature_RNA"
                                                ) +
  
  labs(title = "RNA Counts vs Detected Genes",
       x = "RNA Counts",
       y = "Detected Genes"
       ) +
  
  theme_classic() +
  theme(plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5
                                  ),
    axis.title = element_text(face = "bold",
                              size = 15
                              ),
    axis.text = element_text(size = 13),
    legend.position = "none"
    )


scatter_counts_mt_after <- FeatureScatter(single_cell_singlets,
                                          feature1 = "nCount_RNA",
                                          feature2 = "percent.mt"
                                          ) +
  
  labs(title = "RNA Counts vs Mitochondrial RNA",
       x = "RNA Counts",
       y = "Mitochondrial RNA (%)"
       ) +
  theme_classic() +
  
  theme(plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5
                                  ),
        
    axis.title = element_text(face = "bold",
                              size = 15
                              ),
    
    axis.text = element_text(size = 13),
    legend.position = "none"
    )

qc_scatter_after <-scatter_counts_features_after + scatter_counts_mt_after

qc_scatter_after


# save the plots
ggsave(
  "plots/scrna_plot/qc_scatter_after.pdf",
  plot = qc_scatter_after,
  width = 10,
  height = 15,
  units = "in"
  )



    
#free space 

rm(doublet_summary, doublet_score_plot, doublet_tbl, sce)
rm(single_cell_data2, sce)
# 4 normalization ===========================================================================================================================================

# Normalisation was already done, no need for this step

# perform normalisation 
# single_cell_singlets <- NormalizeData(single_cell_singlets, 
#                                  normalization.method = "LogNormalize", 
#                                scale.factor = 10000)

# run as standalone provide more detail metrics 
single_cell_singlets


# 5 Gene variability calculation ===========================================================================================================================

# calculate the high cell to cell variation 
# In most scRNA-seq experiments only a small proportion of the genes will be informative and biologically variable

# by default 2000 features would be returned
single_cell_singlets <- FindVariableFeatures(single_cell_singlets, 
                                          selection.method = 'vst', 
                                          nfeatures = 2000)


# 6 Scaling ================================================================================================================================================

# next scale the data: Scaling the data will keep highly expressed genes from dominating our analysis -> mean  to 0: variance -> 1 
# scaling can also be used to remove unwanted variations(e.g., cell cycle or percent mitochondrial) 
# Note scaling has also been done. Not required for this analysis 


# 7 Top Variable Gene visualisation ============================================================================================================================

# Identify most highly variable genes

# identify top 10 variable features first
top10 <- head(VariableFeatures(single_cell_singlets), 10)
top10

# create a saving function   

scrna_save_plot_with_ggsave <- function(
    plot_obj,
    filename,
    width,
    height,
    units = "in",
    folder = "plots/scrna_plot/"
) {
  
  # create folder if needed
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

#  plot the most variable genes 
plot1 <- VariableFeaturePlot(single_cell_singlets, pt.size = 4)

# add gene labels for top10
plot2 <- LabelPoints(plot   = plot1,
                     points = top10,
                     repel  = TRUE,
                     xnudge = 0,
                     ynudge = 0,
                     size   = 8
)


# adjust axis labels + title
scrna_top_var_genes <- plot2 + 
  labs(
    x = "Average Expression", 
    y = "Standardized Variance", 
    title = "Top Variable Genes"
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  )

scrna_top_var_genes

# call the save_plot_with_ggsave function 
scrna_save_plot_with_ggsave(
  plot_obj = scrna_top_var_genes,
  filename = "scrna_top_var_genes.pdf",
  width = 12,
  height = 8
)

# clear more space 
rm(plot1,plot2,qc_scatter_before,qc_scatter_after, scatter_counts_features_after, scatter_counts_mt_after, top10)

# 8 PCA Analysis =====================================================================================================================================================

# performing a reduction analysis (PCA)
single_cell_singlets <- RunPCA(single_cell_singlets, 
                            features = VariableFeatures(single_cell_singlets))



# check if PCA results exist  
Reductions(single_cell_singlets)


#calculate standard deviation and variance 
pc_sd <- Stdev(single_cell_singlets,reduction = "pca")

pc_variance <- (pc_sd)^2

head(pc_variance)



pc_percent <- pc_variance / sum(pc_variance) * 100

pc_percent


# generate PC table 

pc_table <- data.frame(PC = paste0("PC", 1:25),
                       Standard_Deviation = round(pc_sd[1:25], 2)
                       )

head(pc_table)




## PC standard  table 
pca_std_tbl <-pc_table %>%
  gt() %>%
  tab_header(
    title = md("**Principal Component Standard Deviations**"),
    subtitle = "First 25 principal components"
  ) %>%
  cols_label(
    PC = "Principal Component",
    Standard_Deviation = "Standard Deviation"
  ) %>%
  fmt_number(
    columns = Standard_Deviation,
    decimals = 2
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_options(
    table.font.size = px(14),
    heading.title.font.size = px(20),
    heading.subtitle.font.size = px(14),
    data_row.padding = px(4)
  )

pca_std_tbl


gtsave(
  pca_std_tbl,
  "plots/scrna_plot/pca_std_tbl.pdf",
  vwidth = 100,
  vheight = 100
)



## PCA Plot using dimplot 

scrrna_pca <- DimPlot(
  single_cell_singlets,
  reduction = "pca",
  group.by = "CellType",
  label = TRUE,
  label.size = 5,
  repel = TRUE
) +
  
  labs(
    x = paste0("PC1 (", round(pc_percent[1], 2), "%)"),
    y = paste0("PC2 (", round(pc_percent[2], 2), "%)"),
    colour = "Cell Type"
  ) +
  
  theme_bw() +
  
  theme(
    axis.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    axis.ticks = element_line(
      linewidth = 1
    ),
    
    axis.ticks.length = unit(0.2, "cm"),
    
    # Legend customisation
    legend.position = "right",
    
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    legend.key.size = unit(0.6, "cm")
  )

scrrna_pca

scrna_save_plot_with_ggsave(
  plot_obj = scrrna_pca,
  filename = "scrrna_pca.pdf",
  width = 8,
  height = 8
)


# to determine the best PCs to choose, use a screeplot/elbowplot

# extract PCA standard deviations manually
# Create custom elbow plot 
scrna_elbow_plot <- ggplot(data.frame(PC = 1:length(pc_sd), 
                                      Stdev = pc_sd), 
                           aes(x = PC, 
                               y = pc_sd)) +
  geom_point(size = 4, 
             alpha = 0.8, 
             color = "steelblue") +
  
  geom_vline(xintercept = 25, 
             linetype = "dashed", 
             color = "red", 
             linewidth = 1.0) +
  
  labs(x = "Principal Components", 
       y = "Standard Deviation", 
       title = "Elbow Plot with chosen PC cutoff") +
  theme(
    plot.title = element_text(hjust = 0.5, 
                              size = 20, 
                              face = "bold"),
    
    axis.title = element_text(size = 16, 
                              face = "bold"),
    
    axis.text = element_text(size = 16, 
                             face = "bold"), 
    
    axis.line = element_line(linewidth = 1.2, 
                             color = "black"), 
    
    panel.background = element_blank(),
    
    panel.grid.major = element_line(color = "grey90", 
                                    linewidth = 0.3),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(0, 50, 5), 
                     expand = c(0, 0)) +
  
  scale_y_continuous(expand = c(0, 0)) +
  
  coord_cartesian(xlim = c(0, 55), 
                  ylim = c(0, max(pc_sd)*1.05))

scrna_elbow_plot

# save the plot 
scrna_save_plot_with_ggsave(
  plot_obj = scrna_elbow_plot,
  filename = "scrna_elbow_plot.pdf",
  width = 10,
  height = 8
)




# 9 clustering ======================================================================================================================================================

# perform clustering 
# Decide “PC = n_pcs” based on elbow, variance, dataset size
# Find neighbors
n_pcs <- 25
single_cell_singlets<- FindNeighbors(single_cell_singlets, dims = n_pcs)

# find clusters - try different resolutions
# lower resolution = fewer clusters, higher = more clusters
# test different resolutions (0.2–1.0)

resolutions = c(0.15, 0.2, 0.5, 0.8,1)

single_cell_singlets <- FindClusters(single_cell_singlets, resolution = resolutions)
 
 
 
 ## run the umap or tSNE
single_cell_singlets <- RunUMAP(single_cell_singlets,
                                 dims = 1:25)
 
 
# number clusters generated by each resolution 
cluster_numbers <- sapply(resolutions,
                           function(x) {length(unique(single_cell_singlets[[paste0("RNA_snn_res.", x)]][, 1]))
                             }
                           )
 

 # generate df  of resolutions 
data.frame(resolution = resolutions,
            n_clusters = cluster_numbers
            )

colnames(single_cell_singlets@meta.data)
 
 #single_cell_singlets$CellType <- as.character(
  # single_cell_singlets$CellType
# )
 
 ## custmoise the plot 
 celltype_colors <- c( "DC3"      = "#984EA3",
                       "act DC3"  = "#D8B7E8",
                       "cDC2"     = "#377EB8",
                       "act cDC2" = "#9EC9E8",
                       "pDC"      = "#666666",
                       "act pDC"  = "#B3B3B3",
                       "mo1"      = "#4DAF4A",
                       "act mo"   = "#A6D854",
                       "moint/2"  = "#F781BF",
                       "cDC1"     = "#E41A1C",
                       "tDC"      = "#66C2A5"
                       )

 
# generate each cell count  
 
# how cell types (clusters) are present 
levels(Idents(single_cell_singlets))
 # the cell counts  per cell type 
table(Idents(single_cell_singlets)) 

cell_counts <- table(single_cell_singlets$CellType)
head(cell_counts)
 

# check cluster names, choose cluster of choice
Idents(single_cell_singlets) <- "CellType"

# Check different clustering results
head(single_cell_singlets@meta.data)


# generate a umap of the cell cluster and counts 

umap_celltype <- DimPlot(single_cell_singlets,
                         reduction = "umap",
                         group.by = "CellType",
                         label = TRUE,
                         label.size = 5,
                         repel = TRUE,
                         cols = celltype_colors
                         ) +
  
  scale_color_manual(values = celltype_colors,
                     labels = paste0(names(cell_counts),
                                     " (n = ",
                                     as.numeric(cell_counts),
                                     ")"
                                     )
                     ) +
  
  labs(title = "UMAP of Single-Cell RNA-seq Data",
       x = "UMAP 1",
       y = "UMAP 2",
       color = "Cell type"
       ) +
  
  
   theme_classic() +
  
   theme(plot.title = element_text(size = 20,
                                   face = "bold",
                                   hjust = 0.5
                                   ),
         
         axis.title = element_text(size = 18,
                                   face = "bold"
                                   ),
         
         axis.text = element_text(size = 14
                                  ),
         
         legend.title = element_text(size = 14,
                                     face = "bold"
                                     ),
         legend.text = element_text(size = 10, 
                                    face = "bold"
                                    )
         )
 
umap_celltype
 
# save  the plot 
scrna_save_plot_with_ggsave(
  plot_obj = umap_celltype,
  filename = "umap_cell_types.pdf",
  width = 10,
  height = 8
)


# 10 perform scRNA-seq differential gene expression (DGE)  ==================================================================================================================

# perform pairs DGE
# use the group_name/cluster for clustering 

# extract the data matrix (layer "data" corresponds to normalized counts in seruat object)
expr_matrix <- GetAssayData(single_cell_singlets, assay = "RNA", layer = "data")

# Run wilcoxauc on the matrix and the cell type labels
de_results <- wilcoxauc(expr_matrix, 
                        y = single_cell_singlets$CellType)

# View top 6 rows of the result 
head(de_results)



# 11 scRNA-seq downstream analysis ==================================================================================================================================

# filter out the non significant genes, log2FC < zero (0) s and AUC < 0.5 

sig_genes <- de_results %>%
  filter(padj < 0.05,
         logFC > 0,
         auc > 0.5)


# 12 create scRNA_seq downstream functions ========================================================================================================================== 

# make a function that will filter for the top genes 
# The function returns df of top genes of interest 

# list containing each cluster dataframe
cluster_top_genes_fun <- function(cluster_list,   
                                  # name of cluster in the list
                                  cluster,   
                                  
                                  # log2 fold change threshold
                                  min_logFC = 0, 
                                  
                                  # AUC threshold
                                  min_auc = 0.0,    
                                  
                                  # % cells expressing gene in cluster
                                  min_pct_in = 1,   
                                  
                                  # % cells expressing gene outside cluster
                                  max_pct_out = 1,  
                                  
                                  # number of genes
                                  top_n = 25) {         
  
  # extract selected cluster
  df <- cluster_list[[cluster]]
  
  # filter and rank genes
  filtered <- df %>%
    filter(
      logFC > min_logFC,
      auc > min_auc,
      pct_in > min_pct_in,
      pct_out < max_pct_out
    ) %>%
    arrange(
      desc(pct_in),
      desc(logFC),
      desc(auc),
      pct_out
    ) %>%
    slice_head(n = top_n)
  
  return(filtered)
}


# create a violin plot function 
plot_violin_genes_plot <- function(
    seurat_obj,
    gene_list,
    main_title = "Gene Expression by Cluster",
    pt_size = 0.1,
    title_size = 18,
    ...
) {
  
  # Create separate violin plots without individual legends
  vln_list <- VlnPlot(
    object = seurat_obj,
    features = gene_list,
    pt.size = pt_size,
    combine = FALSE,
    ...
  )
  
  # Suppress individual legends & remove x axis title
  vln_list <- lapply(vln_list, function(p) {
    p +
      theme(
        legend.position = "none",
        
        # remove "Identity" label
        axis.title.x = element_blank()  
      )
  })
  
  # Combine plots and collect the legend
  combined <- wrap_plots(vln_list, ncol = length(vln_list)) +
    plot_layout(guides = "collect") +
    plot_annotation(title = main_title) &
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = title_size),
      legend.position = "right"
    )
  
  return(combined)
}


# call scRNA_seq downstream functions 

# check column names of filtered data
colnames(sig_genes)

# extract a group of interest 
cDC1_markers <- sig_genes %>%
  filter(group == "cDC1")


head(cDC1_markers)

# split objects into pieces according levels or factor - the clusters (groups)
# create one df for each group 
cluster_list <- split(sig_genes,
                      sig_genes$group
                      )


# view  clusterlist 
glimpse(cluster_list)

# check list names 
names(cluster_list)

# call the filtering function  for each cluster/group (cell type)
# get top 25 marker genes for cluster named "cDC1"

top25_scrna_cDC1 <- cluster_top_genes_fun(cluster_list = cluster_list,
                                          cluster = "cDC1",
                                          min_logFC = 0.5,
                                          min_auc = 0.6,
                                          min_pct_in = 70,
                                          max_pct_out = 27,
                                          top_n = 25
                                          )


head(top25_scrna_cDC1)

top25_scrna_cDC2 <- cluster_top_genes_fun(cluster_list = cluster_list,
                                          cluster = "cDC2",
                                          min_logFC =  0.1,
                                          min_auc =  0.60,
                                          min_pct_in =  60,
                                          max_pct_out = 41
                                          )


top25_scrna_act_cDC2 <- cluster_top_genes_fun(cluster_list = cluster_list,
                                              cluster = "act cDC2",
                                              min_logFC =  0.1,
                                              min_auc =  0.62,
                                              min_pct_in =  60,
                                              max_pct_out = 41
                                              )

top25_scrna_DC3 <- cluster_top_genes_fun(cluster_list = cluster_list,
                                         cluster = "DC3",
                                         min_logFC =  0.2,
                                         min_auc =  0.60,
                                         min_pct_in =  60,
                                         max_pct_out = 41
                                         )

top25_scrna_act_DC3 <- cluster_top_genes_fun(cluster_list = cluster_list,
                                             cluster = "act DC3",
                                             min_logFC =  0.5,
                                             min_auc =  0.70,
                                             min_pct_in =  75,
                                             max_pct_out = 17
                                             )

top25_scrna_pDC <- cluster_top_genes_fun(cluster_list = cluster_list,
                                         cluster = "pDC",
                                         min_logFC =  0.5,
                                         min_auc =  0.70,
                                         min_pct_in =  70,
                                         max_pct_out = 20
                                         )

top25_scrna_act_pDC <- cluster_top_genes_fun(cluster_list = cluster_list,
                                             cluster = "act pDC",
                                             min_logFC =  0.5,
                                             min_auc =  0.70,
                                             min_pct_in =  70,
                                             max_pct_out = 20
                                             )

top25_scrna_tDC <- cluster_top_genes_fun(cluster_list = cluster_list,
                                         cluster = "tDC",
                                         min_logFC =  0.3,
                                         min_auc =  0.68,
                                         min_pct_in =  67,
                                         max_pct_out = 35
                                         )

top25_scrna_mo1 <- cluster_top_genes_fun(cluster_list = cluster_list,
                                         cluster = "mo1",
                                         min_logFC =  0.1,
                                         min_auc =  0.65,
                                         min_pct_in =  62,
                                         max_pct_out = 44
                                         )

top25_scrna_moint_2 <- cluster_top_genes_fun(cluster_list = cluster_list,
                                             cluster = "moint/2",
                                             min_logFC =  0.4,
                                             min_auc =  0.72,
                                             min_pct_in =  0.69,
                                             max_pct_out = 29
                                             )

top25_scrna_act_mo <- cluster_top_genes_fun(cluster_list = cluster_list,
                                            cluster = "act mo",
                                            min_logFC =  0.4,
                                            min_auc =  0.66,
                                            min_pct_in =  0.69,
                                            max_pct_out = 34
                                            )


# combine all the top 25 genes into single tibble 
top25_scrna_all <- bind_rows(top25_scrna_cDC1 %>% 
                               mutate(Cell_type = "cDC1"),
                             top25_scrna_cDC2 %>% 
                               mutate(Cell_type = "cDC2"),
                             top25_scrna_act_cDC2 %>% 
                               mutate(Cell_type = "act cDC2"),
                             top25_scrna_DC3 %>% 
                               mutate(Cell_type = "DC3"),
                             top25_scrna_act_DC3 %>% 
                               mutate(Cell_type = "act DC3"),
                             top25_scrna_pDC %>% 
                               mutate(Cell_type = "pDC"),
                             top25_scrna_act_pDC %>% 
                               mutate(Cell_type = "act pDC"),
                             top25_scrna_tDC %>% 
                               mutate(Cell_type = "tDC"),
                             top25_scrna_mo1 %>% 
                               mutate(Cell_type = "mo1"),
                             top25_scrna_moint_2 %>% 
                               mutate(Cell_type = "moint/2"),
                             top25_scrna_act_mo %>% 
                               mutate(Cell_type = "act mo")
                             )
# preview 
head(top25_scrna_all)


top25_scrna_all <- top25_scrna_all %>%
  dplyr::select(Cell_type,
         feature,
         logFC,
         auc,
         pct_in,
         pct_out,
         padj,
         everything()
         )
# preview 
head(top25_scrna_all)


# select column of interest 
top25_scrna_all <- top25_scrna_all %>%
  dplyr::select(Cell_type:padj)

head(top25_scrna_all)


# extract only gene symbols 
# convert long data to wide 

marker_table <- top25_scrna_all %>%
  arrange(Cell_type, desc(logFC)) %>%
  group_by(Cell_type) %>%
  summarise(
    Genes = paste(feature, collapse = ", "),
    .groups = "drop"
  )

marker_table

# make a table of this genes 

marker_table_gt <- marker_table %>%
  gt() %>%
  tab_header(
    title = md("**Marker Genes by Cell Type**"),
    subtitle = "Genes ranked by log fold-change"
  ) %>%
  cols_label(
    Cell_type = "Cell Type",
    Genes = "Marker Genes"
  ) %>%
  
  # Bold column headers
  tab_style(
    style = cell_text(
      weight = "bold",
      size = px(18)
    ),
    locations = cells_column_labels(everything())
  ) %>%
  
  # Bold Cell Type column
  tab_style(
    style = cell_text(
      weight = "bold",
      size = px(16)
    ),
    locations = cells_body(
      columns = Cell_type
    )
  ) %>%
  
  # Body text
  tab_style(
    style = cell_text(weight = "bold",
      size = px(14)
    ),
    locations = cells_body(
      columns = everything()
    )
  ) %>%
  
  tab_options(
    heading.align = "center",
    table.font.size = px(14),
    heading.title.font.size = px(20),
    heading.subtitle.font.size = px(16)
  )

marker_table_gt




gtsave(
  marker_table_gt,
  "plots/scrna_plot/scrna_marker_table_gt.pdf",
  vwidth = 100,
  vheight = 100
)



# extract the top 3 gene for each cell type 

top3_markers_scrna <-  top25_scrna_all %>%
  group_by(Cell_type) %>%
  slice_max(order_by = logFC,
            n = 3
            ) %>%
  arrange(Cell_type, 
          desc(logFC))

head(top3_markers_scrna)



# remove the adjusted p-values 
top3_markers_scrna <- top3_markers_scrna %>% select(-c(padj))

head(top3_markers_scrna) 

top3_markers_scrna_gt <- top3_markers_scrna %>%
  gt(
    rowname_col = "feature",
    groupname_col = "Cell_type"
  ) %>%
  tab_header(
    title = md("**Top 3 Marker Genes by Cell Type**"),
   # subtitle = "Wilcoxon rank-sum differential expression analysis"
  ) %>%
  cols_label(
    feature = "Gene",
    logFC = "logFC",
    auc = "AUC",
    pct_in = "% In",
    pct_out = "% Out"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  fmt_number(
    columns = c(logFC, auc, pct_in, pct_out),
    decimals = 2
  ) %>%
 
  tab_options(
    heading.align = "center",
    table.font.size = px(14), # Reduced slightly from 16
    heading.title.font.size = px(18), # Reduced slightly from 20
    heading.subtitle.font.size = px(13),
    # --- REDUCE VERTICAL HEIGHT HERE ---
    data_row.padding = px(0.3), # Default is px(6)
    column_labels.padding = px(1),
    row_group.padding = px(1),
    heading.padding = px(1)
  )



top3_markers_scrna_gt


gtsave(
  top3_markers_scrna_gt,
  "plots/scrna_plot/top3_markers_scrna_gt.pdf",
  vwidth = 100,
  vheight = 100
)

# find common genes btw the technologies 

# the "cell_gene_table" is obtained from bulk RNA-seq 
# it contains top25 genes for each cell type 

# the top 25 genes from bulk RNA-seq 
head(top_bulk_genes) 

# top 25 genes from scrna RNA-seq  
head(top25_scrna_all) 

compare_gene_lists <- function(list_bulk, list_scrna) {
  
  # find common genes
  common_genes <- intersect(list_bulk, list_scrna)
  
  cat("Genes common between bulk and scRNA:\n")
  print(common_genes)
  
  # genes in bulk not in scRNA
  not_in_scrna <- setdiff(list_bulk, list_scrna)
  cat("\nGenes in bulk but not in scRNA:\n")
  print(not_in_scrna)
  
  # genes in scRNA not in bulk
  not_in_bulk <- setdiff(list_scrna, list_bulk)
  cat("\nGenes in scRNA but not in bulk:\n")
  print(not_in_bulk)
  
  # return a list with all three sets for further use
  invisible(list(common = common_genes,
                 bulk_only = not_in_scrna,
                 scrna_only = not_in_bulk))
}







comparison_results <- list()

cell_types <- intersect(unique(top_bulk_genes$Cell),
                        unique(top25_scrna_all$Cell_type)
                        )


# loop through the cell types 
for (cell_type in cell_types) {
  
  bulk_genes <- top_bulk_genes %>%
    filter(Cell == cell_type) %>%
    pull(Gene)
  
  scrna_genes <- top25_scrna_all %>%
    filter(Cell_type == cell_type) %>%
    pull(feature)
  
  comparison_results[[cell_type]] <- compare_gene_lists(
    list_bulk = bulk_genes,
    list_scrna = scrna_genes
  )
}



## create a table 

comparison_table <- purrr::map_dfr(
  names(comparison_results),
  function(cell_type) {
    
    result <- comparison_results[[cell_type]]
    
    tibble(Cell_type = cell_type,
           Common = paste(result$common, collapse = ", "),
           Bulk_only = paste(result$bulk_only, collapse = ", "),
           scRNA_only = paste(result$scrna_only, collapse = ", ")
    )
    }
  )


comparison_table


# make the table publishable 
comparison_table_gt <- comparison_table %>%
  gt() %>%
  
  tab_header(
    title = md("**Comparison of Top Marker Genes Between Bulk RNA-seq and scRNA-seq**"),
    subtitle = "Comparison within corresponding cell types"
  ) %>%
  
  cols_label(
    Cell_type = "Cell Type",
    Common = "Common Genes",
    Bulk_only = "Bulk RNA-seq Only",
    scRNA_only = "scRNA-seq Only"
  ) %>%
  
  # Bold column headers
  tab_style(
    style = cell_text(
      weight = "bold",
      size = px(14)
    ),
    locations = cells_column_labels(everything())
  ) %>%
  
  # Bold cell type
  tab_style(
    style = cell_text(
      weight = "bold",
      size = px(12)
    ),
    locations = cells_body(columns = Cell_type)
  ) %>%
  
  # Body text
  tab_style(
    style = cell_text(weight = "bold",
                      size = px(10)
    ),
    locations = cells_body(
      columns = c(Common, Bulk_only, scRNA_only)
    )
  ) %>%
  
  tab_options(
    table.font.size = px(18),
    heading.align = "center",
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(14),
    table.width = pct(70)
  )

comparison_table_gt

gtsave(
  comparison_table_gt,
  filename = "plots/scrna_plot/comparison_table_gt.pdf",
  vwidth = 25000,
  vheight = 1000
)


# shared gens table only 
# combine the bulk RNA-seq and scrna 
overlap_genes <- inner_join(top_bulk_genes,
                            top25_scrna_all,
                            by = c(
                              "Cell" = "Cell_type",  # column name of the bulk RNA-seq table 
                              "Gene" = "feature"     # column names of scrna table 
                              )
                            )


overlap_genes

# create a list of bulk-RNA-seq top 25 genes per cell 
# check colnames of the overlap_genes object 
colnames(overlap_genes)

# extract the cell type and gene symbol columns 
overlap_genes_final_list <- overlap_genes %>%
  dplyr::select(Cell:Gene)

# pivot the overlap_gene final list wide to generate a table 
overlap_gene_table_wide <- overlap_genes_final_list%>%
  group_by(Cell) %>%
  summarise(
    Shared_genes = paste(Gene, collapse = ", "),
    Number_of_genes = n(),
    .groups = "drop"
    ) %>%
  arrange(Cell)

# Format  and customise the gt table
gene_overlap_gt <- overlap_gene_table_wide %>%
  gt() %>%
  
  # rename column headers
  cols_label(
    Cell = "Cell Type",
    Shared_genes = "Shared Genes",
    Number_of_genes = "Number of Shared Genes"
    ) %>%
  
  # left align text
  cols_align(
    align = "left",
    columns = everything()
    ) %>%
  
  # bold cell type column
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Cell)
    ) %>%
  
  # header formatting
  tab_style(
    style = list(
      cell_fill(color = "grey90"),
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels()
    ) %>%
  
  # table options
  tab_options(
    table.font.size = 18,
    data_row.padding = px(4),
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    table_body.hlines.color = "grey80"
    ) %>%
  
  # Alternating row shading
  opt_row_striping()


gene_overlap_gt

gtsave(
  gene_overlap_gt,
  filename = "plots/scrna_plot/bulk_vs_scRNA_shared_genes.pdf",
  vwidth = 2000,
  vheight = 1500
)




# function call 
# genes to be plotted 
genes = c("CLEC9A","FCER1A","JCHAIN")

vln_shared_right <- plot_violin_genes_plot(seurat_obj = single_cell_singlets,
                                           gene_list = genes,
                                           main_title = "Bulk Marker Gene Expression",
                                           pt_size = 0.1
)

print(vln_shared_right)

scrna_save_plot_with_ggsave(
  plot_obj = vln_shared_right,
  filename = "Bulk_marker_gene_violin_plot.pdf",
  width = 15,
  height = 6
)

# 9. statistical performance of the genes ----

#statistical performance of each genes

stat_table_fun <- function(df, 
                           group,
                           min_auc = 0,
                           min_pct_in = 0,
                           max_pct_out = 1) {
  
  # Filter by criteria & cluster
  stat_info <- df %>%
    filter(auc > min_auc) %>%           # select high AUC markers :contentReference[oaicite:1]{index=1}
    filter(pct_in > min_pct_in) %>%
    filter(pct_out < max_pct_out) %>%
    filter(group == group) %>%
    arrange(pct_out, auc, logFC)
  
  # Select columns of interest
  stat_table <- stat_info %>%
    dplyr::select(
      feature,
      logFC,
      auc,
      padj,
      pct_in,
      pct_out
    ) %>%
    head(25)
  
  return(stat_table)
  }


#call the function 

cell_table_stat <- stat_table_fun(df = sig_genes,
                                  group = "cDC1",
                                  min_auc = 0.741,
                                  min_pct_in = 60,
                                  max_pct_out = 25
                                  )

print(cell_table_stat)


#scatter plot of AUC vs pct_in 
pct_in_cut <- 15
auc_cut    <- 0.78

cDC1_scrna_gene_performance <- ggplot(cell_table_stat, aes(x = pct_in, y = auc)) +
  geom_point(aes(color = pct_out), 
             size = 4) +
  geom_vline(xintercept = pct_in_cut, 
             linetype = "dashed", 
             color = "blue",
             size = 1.5) +
  geom_hline(yintercept = auc_cut,     
             linetype = "dashed", 
             color = "red",
             size = 1.5) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    x = "% of cells in cluster expressing gene (pct_in)",
    y = "AUC (marker discrimination)",
    color = "% Expression in other clusters",
    title = "Top cDC1 Gene Perfamnce"
    ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20,
                              hjust = 0.5, 
                              face = "bold"),
    axis.title   = element_text(face = "bold", 
                                size = 16),       
    axis.text    = element_text(face = "bold", 
                                size = 12)       
    ) +
  geom_text(
    data = subset(cell_table_stat, pct_in >= pct_in_cut & auc >= auc_cut),
    aes(label = feature),
    size = 5,
    vjust = -1
    )

cDC1_scrna_gene_performance 


# save the figure 

scrna_save_plot_with_ggsave(plot_obj = cDC1_scrna_gene_performance,
                      filename = "cDC1_scrna_gene_performance.pdf",
                      width = 15,
                      height = 15)





make_styled_stat_table <- function(stat_df,
                                   highlight_condition,
                                   caption = "",
                                   digits = 2,
                                   bootstrap_options = c("striped", "hover"),
                                   latex_options = c("striped", "hold_position"),
                                   highlight_bg = "#DFF0D8") {
  
  # generate row indices matching the highlight condition
  # stat_df must be the same data used for table
  row_idx <- which(eval(substitute(highlight_condition), stat_df, parent.frame()))
  
  # make and style kable
  stat_df %>%
    kable(
      caption = "<b> CDC1 Marker Performance </b>",
      digits = digits,
      format = "html"
      ) %>%
    kable_styling(
      full_width = FALSE,
      bootstrap_options = bootstrap_options,
      latex_options = latex_options
    ) %>%
    row_spec(
      row_idx, 
      background = highlight_bg
    )
  }


make_styled_stat_gt <- function(stat_df,
                                highlight_condition,
                                caption = "Marker Statistics Performance For Each Cell",
                                digits = 2,
                                highlight_bg = "#DFF0D8") {
  
  # identify rows to highlight
  row_idx <- which(
    eval(substitute(highlight_condition), stat_df, parent.frame())
  )
  
  # create gt table
  gt_table <- stat_df %>%
    gt() %>%
    
    # table title
    tab_header(
      title = caption
    ) %>%
    
    # round numeric columns
    fmt_number(
      columns = where(is.numeric),
      decimals = digits
    ) %>%
    
    # highlight selected rows
    tab_style(
      style = list(
        cell_fill(color = highlight_bg)
      ),
      locations = cells_body(
        rows = row_idx
      )
    ) %>%
    
    # bold column headers
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    
    # table options
    tab_options(
      table.font.size = 12,
      table.border.top.color = "black",
      table.border.bottom.color = "black",
      table_body.hlines.color = "grey80"
    )
  
  return(gt_table)
  }



# call the function 
styled_table <- make_styled_stat_table(stat_df = cell_table_stat,
                                       highlight_condition = auc > 0.84 & pct_out < 12,
                                       caption = "Marker statistics Performance For Each Cell",
                                       digits = 2,
                                       highlight_bg = "#DFF0D8"
                                       )

styled_table





# call the make_styled_stat_gt function 
styled_table <- make_styled_stat_gt(
  stat_df = cell_table_stat,
  highlight_condition = auc > 0.84 & pct_out < 12,
  caption = "cDC1 Marker Performance",
  digits = 2
  )

# print the table 
styled_table


# the statistics performance 
gtsave(
  styled_table,
  filename = "plots/scrna_plot/cDC1_marker_statistics.pdf"
  )









