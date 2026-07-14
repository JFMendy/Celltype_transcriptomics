

# ============================================= Single cell Analysis =========================================================================

# 0.1 load the libraries  and data -----
# libraries 
{library(Seurat)
  library(presto)
  library(gt)
  library(dplyr)
}


## read in the scRNA-seq  data
single_cell_data <- readRDS( file = "C:/Users/user/Desktop/MSc_Theses/DIT_MSc_thesis/rawData/scRNASeq_annotated/hsobj_20221209_annotated.rds")


# 1 preprocess the data -------

single_cell_data2 <- single_cell_data      # make copy of the data  
#

# 2 scRNA-seq inspection ------
glimpse(single_cell_data2)                 # check the structure of the seurat object 
head(single_cell_data2@meta.data)          #using head to return only the first 6 rows  

dim(single_cell_data)
head(rownames(single_cell_data), 10)       # check 10 row  names 


# 3 Quality control and cleaning ----
# clean the data by removing: dying cells -> exhibit high mitochondrial contamination 
# by removing genes starting with MT -> mitochondrial genes 
# store mitochondrial percentage in object meta data in "percent.mt" column 

single_cell_data2[["percent.mt"]] <- PercentageFeatureSet(single_cell_data2, pattern = "^MT-")

## 3.1 Pre cleaning visualization ---- 
# Visualize QC Metrics  using the Violin plots
VlnPlot(single_cell_data2, features = c("nFeature_RNA", 
                                        "nCount_RNA",
                                        "percent.mt"),
        ncol = 3, 
        pt.size = 0.1) 



# Scatter plots used to observed the relationship btw genes or percent_MT vs counts 
FeatureScatter(single_cell_data2, 
               feature1 = "nCount_RNA", 
               feature2 = "nFeature_RNA")

FeatureScatter(single_cell_data2, 
               feature1 = "nCount_RNA", 
               feature2 = "percent.mt") 



# data filtering 
#filter the genes and counts  
#remove low quality cells (dead or dying cells), cell-free RNA, or doublets)

# Based on these plots, define thresholds — e.g.:
min_genes <- 200
max_genes <- 6000       # adjust depending on dataset
max_mt    <- 20         # or lower if appropriate
min_counts <- 500       # optional lower bound

# cleaning the data 
single_cell_data2<- subset(single_cell_data2,
                           subset = nFeature_RNA > min_genes &
                             nFeature_RNA < max_genes &
                             percent.mt < max_mt &
                             nCount_RNA > min_counts)




## 3.2 Post cleaning Visualisation  ----
#filtering verification after preprocessing 
#post filtering violin plot 
VlnPlot(single_cell_data2, features = c("nFeature_RNA",
                                        "nCount_RNA", 
                                        "percent.mt"))


FeatureScatter(single_cell_data2,
               "nCount_RNA",
               "percent.mt") 


# 4 normalization ----

# Normalisation was already done, no need for this step
#perform normalisation 
#single_cell_data2 <- NormalizeData(single_cell_data2, 
#                                  normalization.method = "LogNormalize", 
#                                scale.factor = 10000)

#number for genes
nrow(single_cell_data2) # number of genes
ncol(single_cell_data2)  # number of cells 
single_cell_data2

# 5 Gene variability calculation -----
# calculate the high cell to cell variation 
# In most scRNA-seq experiments only a small proportion of the genes will be informative and biologically variable
# by default 2000 features would be returned
single_cell_data2 <- FindVariableFeatures(single_cell_data2, 
                                          selection.method = 'vst', 
                                          nfeatures = 2000)

# 8 Scaling 
# next scale the data: Scaling the data will keep highly expressed genes from dominating our analysis -> mean  to 0: variance -> 1 
# scaling can also be used to remove unwanted variations(e.g., cell cycle or percent mitochondrial) 
# Note scaling has also been done. Not required for this analysis 

## 5.1 Top Variable Gene visualisation ----
# Identify most highly variable genes
# identify top 10 variable features first
top10 <- head(VariableFeatures(single_cell_data2), 10)
top10

## 5.2 Top variability visualisation  ----

scrna_save_plot_with_ggsave <- function(
    plot_obj,
    filename,
    width,
    height,
    units = "in",
    folder = "plots/scrna_plot/"
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

# basic VariableFeaturePlot
plot1 <- VariableFeaturePlot(single_cell_data2, pt.size = 4)

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

# 6 PCA Analysis ----

#performing a reduction analysis (PCA)
single_cell_data2 <- RunPCA(single_cell_data2, 
                            features = VariableFeatures(single_cell_data2))
# to determine the best PCs to choose, use a screeplot/elbowplot
# Visualize variance explained to choose PCs
# Visualize variance explained to choose PCs - LARGER DOTS

# Extract PCA standard deviations manually
pca_stdev <- Stdev(object = single_cell_data2, reduction = "pca")

# Create custom elbow plot with BIGGER DOTS - BOLD AXES
scrna_elbow_plot <- ggplot(data.frame(PC = 1:length(pca_stdev), 
                                      Stdev = pca_stdev), 
                           aes(x = PC, 
                               y = Stdev)) +
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
                  ylim = c(0, max(pca_stdev)*1.05))

scrna_elbow_plot

# save the plot 
scrna_save_plot_with_ggsave(
  plot_obj = scrna_elbow_plot,
  filename = "scrna_elbow_plot.pdf",
  width = 10,
  height = 8
)

# 7 clustering ----

# perform clustering 
# Decide “n_pcs” based on elbow, variance, dataset size
# Find neighbors
n_pcs <- 25
single_cell_data2<- FindNeighbors(single_cell_data2, dims = n_pcs)

# find clusters - try different resolutions
# lower resolution = fewer clusters, higher = more clusters
# test different resolutions (0.2–1.0)
#single_cell_data2 <- FindClusters(single_cell_data2, resolution = c(0.2, 0.5, 0.8,1))

single_cell_data2 <- FindClusters(single_cell_data2,
                                  , resolution = 0.15)  

# Check different clustering results
head(single_cell_data2@meta.data)

# single_cell_data2<- RunUMAP(single_cell_data2, dims = 1:25)
# check cluster names, choose cluster of choice
Idents(single_cell_data2) <- "CellType"
# how cell types are present 
levels(Idents(single_cell_data2))
# the cell counts  per cell type 
table(Idents(single_cell_data2))

# Visualise the UMAP clusters 
umap_plot <- DimPlot(single_cell_data2,
                     reduction = 'umap',
                     label = TRUE,
                     label.size = 6,
                     pt.size = 0.01) + 
  theme_bw() +
  labs(
    x = "UMAP 1",        # custom x-axis label
    y = "UMAP 2"         # custom y-axis label
    # title = "Optional Title"   # uncomment if you want a title
  ) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),   # axis labels (UMAP 1, UMAP 2)
    axis.text = element_text(size = 14, face = "bold")     # tick numbers
  )
umap_plot
scrna_save_plot_with_ggsave(
  plot_obj = umap_plot,
  filename = "umap_cell_types.pdf",
  width = 10,
  height = 8
)

# 13 scRNA-seq differential gene expression (DGE)

#perform pairs DGE
#use the group_name/cluster for clustering 
# Extract the data matrix (layer "data" corresponds to normalized counts)
expr_matrix <- GetAssayData(single_cell_data2, assay = "RNA", layer = "data")

# Run wilcoxauc on the matrix and the cell type labels
de_results <- wilcoxauc(expr_matrix, 
                        y = single_cell_data2$CellType)
# View top 6 genes
head(de_results)

# 8 scRNA-seq downstream analysis ----
#select the significant genes  with an AUC greater 0.5
sig_genes <- de_results %>%
  filter(padj < 0.05,
         logFC > 0,
         auc > 0.5)

## 8.1 scRNA_seq functions -----
# make a function that will filter for the top genes 
# The function returns df of top genes of interest 
cluster_top_genes_fun <- function(cluster_list,        # list containing each cluster dataframe
                                  cluster,             # name of cluster in the list
                                  min_logFC = 0,        # log2 fold change threshold
                                  min_auc = 0.0,        # AUC threshold
                                  min_pct_in = 1,       # % cells expressing gene in cluster
                                  max_pct_out = 1,      # % cells expressing gene outside cluster
                                  top_n = 25) {         # number of genes
  
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

# create a volin plot function 
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
        axis.title.x = element_blank()  # remove "Identity" label
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

## 8.2 scRNA_seq function call ----- 

colnames(sig_genes)

# extract a group of interest 
cDC1_markers <- sig_genes %>%
  filter(group == "cDC1")


cluster_list <- split(
  sig_genes,
  sig_genes$group
)
# view  clusterlist 
glimpse(cluster_list)
# check list names 
names(cluster_list)

# call the filtering function 
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


top25_scrna_all <- top25_scrna_all %>%
  select(Cell_type,
         feature,
         logFC,
         auc,
         pct_in,
         pct_out,
         padj,
         everything()
         )

head(top25_scrna_all)

top25_scrna_all <- top25_scrna_all %>%
  select(Cell_type:padj)

head(top25_scrna_all)

# find common genes btw the technologies 

cell_gene_table  # the bulk RNA-seq table 
top25_scrna_all   # scrna RNA-seq table 

# combine the bulk RNA-seq and scrna 
overlap_genes <- inner_join(cell_gene_table,
                            top25_scrna_all,
                            by = c(
                              "Cell" = "Cell_type",  # column name of the bulk RNA-seq table 
                              "Gene" = "feature"     # column names of scrna table 
                              )
                            )

# create a list of bulk-RNA-seq top 25 genes per cell 
# check colnames of the overlap_genes object 
colnames(overlap_genes)

# extract the cell type and gene symbol columns 
overlap_genes_final_list <- overlap_genes %>%
  select(Cell:Gene)

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
    table.font.size = 12,
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

vln_shared_right <- plot_violin_genes_plot(seurat_obj = single_cell_data2,
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
  filename = "plots/scrna_plot/cDC1_marker_statistics.png"
  )






