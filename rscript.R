
##### package load #########################################################################################################
{library(DESeq2)
 library(biomaRt)
  
library(tidyverse)
library(factoextra)
library(lolliplot)
library(gridExtra) 
library(grid)
  
  
  
  }

##### read in the data #####################################################################################################
mydata_bulk_RNA <- readRDS(file = "C:/Users/user/Desktop/MSc_Theses/DIT_MSc_thesis/rawData/bulkRNASeq/se.rds" )


##### inspect the data and extract the count data set ##################################################################### 
# make a copy of original data 
mydata_bulk_RNA2 <- mydata_bulk_RNA

# check the assays available in the summerisedExperiment 
assayNames(mydata_bulk_RNA2)  

# colData column  and rowData names 
names(colData(mydata_bulk_RNA2)) 

names(rowData(mydata_bulk_RNA2))

# check the column names  
colData(mydata_bulk_RNA)

# check the row names 
rowData(mydata_bulk_RNA)   

# show top 6 rows 
head(colData(mydata_bulk_RNA2)) 

# dimensions should match original
dim(mydata_bulk_RNA2) 

# summary of data file 
glimpse(mydata_bulk_RNA2) 

# detail structure of the data file 
str(mydata_bulk_RNA2)


#### data wrangling and clean ##############################################################################################

#check the sample information  
colData(mydata_bulk_RNA2)
# remove  columns of no interest 
colData(mydata_bulk_RNA2)$...7 <- NULL
colData(mydata_bulk_RNA2)$We.have.4.donors..Donor.252..253..254..255 <- NULL

#check for duplicated ids 
any(duplicated(rowData(mydata_bulk_RNA2)$id))

# extract the duplicates and view the duplicates and view 
duplicates <- (rowData(mydata_bulk_RNA2)$id[duplicated(rowData(mydata_bulk_RNA2)$id)])  
duplicates 

# remove the duplicate
mydata_bulk_RNA2 <- mydata_bulk_RNA2[!duplicated(rowData(mydata_bulk_RNA2)$id), ]

#confirm duplicate removal 
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

# double check again 
str(colData(mydata_bulk_RNA2))


# remove genes with less than 10 row sum count and sample less than 3 
counts <- assay(mydata_bulk_RNA2, "counts")

keep <- rowSums(counts >= 10) >= 3
mydata_bulk_RNA2 <- mydata_bulk_RNA2[keep, ]

# check if the row and column names matach
all(colnames(assay(mydata_bulk_RNA2, "counts")) == 
      rownames(colData(mydata_bulk_RNA2)))


#### exploratory analysis  and plot generation  ########################################################################## 

# extract data of interest to generate demograpgic data and gene count profile per sample 
meta <- as.data.frame(colData(mydata_bulk_RNA2))
head(meta)

#tabulate study participants
table(meta$Patient, meta$Timepoint.after.vaccination)

#number of samples per cell 
table(meta$Cell.population)

#gene count per total sample 
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

#print the table 
tail(demo_table)

# generate a vector called totals  to store the row totals 
totals <- demo_table %>%
  dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
  dplyr::mutate(Patient = "Total") %>%
  dplyr::select(Patient, dplyr::everything())
totals

# join the row totals to the demo_table  
demo_table_total <- dplyr::bind_rows(demo_table, totals)
demo_table_total 

#generate a more polish demogragraphic table 


#  convert standard table to a grob for plotting
#load the library(gridExtra) library(grid)

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

#save the plot 
# wrap the plot in grid_arrange for ggplot to open it easily. 
# Create a grob object for saving
sample_demo_table_plot2 <- gridExtra::arrangeGrob(sample_demo_table)

# Save the table as PDF
ggsave(
  filename = "plots/sample_demo_table.pdf",
  plot = sample_demo_table_plot2,
  width = 10,
  height = 6
)

# extract the gene counts information 
## generate gene count plot
# Use the gene_count_per_sample table 
# order samples by total_gene_counts for nicer visualization
# add a column called total_gene_count 
meta$total_gene_counts <- colSums(counts)
#create table: sampleName, cell type, total counts
gene_count_per_sample <- meta %>%
  dplyr::select(sampleName, Cell.population, total_gene_counts) %>%    #  needed columns 
  dplyr::arrange(Cell.population, total_gene_counts) %>%                 # order base on  total_gene_counts
  dplyr::mutate(sampleName = factor(sampleName, levels = sampleName))    #  create new column name for sample

head(gene_count_per_sample)

#remove any missing value or NAs
# prepare the data for plotting 

gene_count_clean <- gene_count_per_sample %>%
  filter(!is.na(total_gene_counts) & !is.na(Cell.population))

#generate a plot for the gene count per sample 
gene_count_plot <- ggplot(gene_count_clean, 
       aes(x = sampleName, 
           y = total_gene_counts, 
           color = Cell.population)) +
  
  # lollipop sticks
  geom_segment(aes(xend = sampleName, y = 0, yend = total_gene_counts),
               linewidth = 1.2) +
  
  # lollipop heads
  geom_point(size = 3) +
  
  # titles and axis labels
  labs(
    title = "Gene Counts per Sample per Cell Type",
    x = "Sample ID",
    y = "Total Gene Counts per Sample",
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
  
  #  select theme: subtle grids, boxed axes
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    
    axis.title = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    
    # subtle gridlines
    panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor = element_line(color = "gray90", linewidth = 0.15),
    
    # boxed axes
    axis.line = element_line(color = "black", linewidth = 1),
    
    # legend styling
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 12, face = "bold"),
    legend.key = element_rect(fill = "white", color = "black"),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.y = unit(0.3, "cm")
  ) +
  
  #customise the legend symbols
  guides(
    color = guide_legend(
      override.aes = list(size = 4),  # bigger legend points
      ncol = 1                        # number of columns
    )
  )

gene_count_plot

#save the plot 

ggsave(
  filename = "plots/gene_count_plot.pdf",
  plot = gene_count_plot,
  width = 35,
  height = 20,
  units = "in"
)

# arrange liliplot in ascending order 
# extract the count matrixs 
df_counts <- assay(mydata_bulk_RNA2, "counts")


# turn into a tibble in long format
df_long <- as_tibble(df_counts, rownames = "gene") %>% 
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "count"
  )



# Compute QC metrics per sample
df_detected_qc_metrics <- df_long %>%
  group_by(sample) %>%
  summarize(library_size = sum(count),               # total counts per sample
            n_genes = sum(count > 0),                # number of genes detected
            n_zero = sum(count == 0),                # number of zero-count genes
            pct_zero = mean(count == 0) * 100,       # percent of genes with zero count
            mean_count = mean(count),                # mean count per gene
            median_count = median(count),            # median count per gene
            .groups = "drop"
  )

#Library Size Barplot: Samples with lower total counts may indicate low sequencing depth or a problem with that sample
df_detected_qc_metrics




# convert the sample column to a factor levels 
df_detected_qc_metrics <- as.data.frame(df_detected_qc_metrics)

df_detected_qc_metrics$sample<- as.factor(df_detected_qc_metrics$sample
)
str(df_detected_qc_metrics)

# load the forcats library 
# order the plots based on gene counts
df_detected_qc_metrics$sample <- fct_reorder(df_detected_qc_metrics$sample, 
                                             df_detected_qc_metrics$n_genes)



#plot the data as Lollipop plot 
df_detected_qc_metrics_plot <- ggplot(df_detected_qc_metrics, 
                                      aes(x = sample, 
                                          y = n_genes)) +
  
  # Lollipop sticks - thinner and lighter
  geom_segment(aes(xend = sample,  
                   y = 0, 
                   yend = n_genes),
               color = "black",        # line colours
               linewidth = 0.2) +        # line thickness
  
  # Points - gradient colors with viridis scale
  geom_point(aes(fill = n_genes),  # Color by gene count value
             size = 2,                # size of the point 
             shape = 21,              # filled circles with outline
             color = "black",
             stroke = 0.2)  +         # outline thickness 
  #background 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,  
                                   hjust = 1,
                                   vjust = 1, 
                                   size = 3, 
                                   face = 'bold'),
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = 25,
                                  hjust = 0.5, # horizontal centering 
                                  margin = margin(b = 10)),
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
        plot.margin = margin(15, 15, 15, 15)) + # margin round the plot 
  
  #y-axis scaling 
  scale_y_continuous(
    limits = c(0, 
               max(df_detected_qc_metrics$n_genes) * 1.1),
    expand = c(0, 0),
    labels = scales::comma  # Format numbers with commas
  )+
  
  # Labels
  labs(
    title = "Gene Detection Rate per Sample",
    x = "Samples (Ordered by Gene Count)",
    y = "Number of Detected Genes",
    caption = "Yellow Fever: Bulk RNA-Seq Dataset"
  )

df_detected_qc_metrics_plot 

# save the plot 

ggsave(
  filename = "plots/df_detected_qc_metrics_plot.pdf",
  plot = df_detected_qc_metrics_plot,
  width = 30,
  height = 20,
  units = "in"
)


####### Principle component analysis ###########################################


#Normalisation 
# extract TPM matrix
tpms <- assay(mydata_bulk_RNA2, 'tpms')
# filter out zero values 
tpms <- tpms[rowSums(tpms != 0) > 0, ]
sum(is.nan(tpms))

# log-transform + constant value 
mat <- log2(tpms + 0.55)  

# transpose so samples are rows
mat_t <- t(mat)
# perform PCA
pca_result <- prcomp(mat_t)  # already scaled
summary(pca_result)

# extract PC scores + join metadata
pca_df <- as.data.frame(pca_result$x)  # one row per sample
pca_df$sampleName <- rownames(pca_df)

#meta data extraction 
metadata_df <- as.data.frame(colData(mydata_bulk_RNA2))

pca_df <- merge(pca_df, metadata_df, by = "sampleName")  # metadata_df holds sample annotations

# Plot PC1 vs PC2
#round up the % values 
percentVar <- round(100 * (pca_result$sdev^2 / sum(pca_result$sdev^2)), 1)

#provide colors for the cells 
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "black", "#AD7700", "#666666")



pca_cell_population_plot <- ggplot(pca_df, 
                               aes(x = PC1, 
                                   y = PC2, 
                                   colour = Cell.population )) +
  geom_point(size = 3) +
  
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  geom_point(size = 3) +
  
  stat_ellipse(type = "t",         # Add ellipses with multivariate t-distribution
               level = 0.95,       # 95% confidence level      
               linewidth = 1,       # Line thickness        
               linetype = "solid",  # Line type
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

#save the plot 
ggsave(
  filename = "plots/pca_cell_population_plot.pdf",
  plot = pca_cell_population_plot,
  width = 30,
  height = 20,
  units = "in"
)





#scree plot 
get_pca_variance_flags <- function(pca_obj,
                                   prop_threshold = 5,       # per-PC: minimum % variance
                                   cum_threshold = 80,       # cumulative: minimum % for retained PCs
                                   max_N = NULL,             # optionally limit to first N PCs
                                   percent = TRUE) {
  # pca_obj: result of prcomp()
  # prop_threshold: minimum % variance explained by a PC to be flagged (numeric, e.g. 5 means 5%)
  # cum_threshold: cumulative % variance threshold for retention (numeric, e.g. 80)
  # max_N: if not NULL, only consider first max_N PCs
  # percent: whether to return % (default = TRUE) or fraction (FALSE)
  
  sdev <- pca_obj$sdev
  vars  <- sdev^2
  prop  <- vars / sum(vars)
  if (percent) prop <- prop * 100
  cum   <- cumsum(prop)
  
  # If max_N specified, subset
  if (!is.null(max_N)) {
    prop <- prop[1:max_N]
    cum  <- cum[1:max_N]
  }
  
  df <- data.frame(
    PC                = seq_along(prop),
    var_explained     = prop,
    cumvar_explained  = cum,
    flag_prop         = prop >= prop_threshold,
    flag_cum          = cum <= cum_threshold
  )
  
  # another useful flag: keep PCs up to the first where cumulative variance exceeds the threshold
  df$flag_cum_cutoff <- df$PC <= min(which(cum >= cum_threshold))
  
  return(df)
}

# extract the tpms 
tpms <- assay(mydata_bulk_RNA2, 'tpms')
# filter out zero values 
tpms <- tpms[rowSums(tpms != 0) > 0, ]
sum(is.nan(tpms))

# log-transform + constant value 
mat <- log2(tpms + 0.55)  

# transpose so samples are rows
mat_t <- t(mat)

# perform PCA
pca_res <- prcomp(mat_t)



# function call
var_table <- get_pca_variance_flags(pca_res,
                                    prop_threshold = 5,
                                    cum_threshold = 75,
                                    max_N = 20)


# var_explained: the variance explained by each PC (in percent if percent = TRUE)

# cumvar_explained: the cumulative variance explained up to current PC

#flag_prop: logical — TRUE if that PC explains ≥ prop_threshold % of total variance.

#flag_cum: logical — TRUE for PCs whose cumulative variance is still ≤ cum_threshold, i.e. before the “cut-off” (less conservative).

#flag_cum_cutoff: logical — TRUE for all PCs up to (and including) the first PC where cumulative variance ≥ cum_threshold. — this often gives a minimal set of PCs that achieves the desired cumulative variance.

print(var_table)




# Scree plot 
eigen_plot <- fviz_eig(pca_result,
                       choice  = "variance",
                       ncp     = 15,                 # only first 15 PCs
                       geom    = c("bar", "line"),
                       barfill = "navy",
                       barcolor= "black",
                       linecolor = "red",
                       addlabels = TRUE,
                       hjust   = -0.3,
                       size = 10 ) +            # adjust label position
  ylim(0, 20) +                         # y-axis limit: up to 50% variance
  geom_hline(yintercept = 5,            # horizontal line at 5% threshold
             linetype = "dashed", color = "gray50",
             size = 1) +
  labs(
    title = "Scree Plot — Variance Explained",
    x = "Principal Components (PC)",
    y = "Percent Variance Explained"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    axis.title   = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 45, hjust = 1, 
                                face = "bold", color = "black", size = 12 ),
    axis.text.y = element_text(face = "bold", color = "black", size = 11)
    
  )

print(eigen_plot)





######### Differential gene expression analysis ################################ 
# pairwise analysis was performed 1 vs 1
# change coldata column name 
names(colData(mydata_bulk_RNA2))
names(rowData(mydata_bulk_RNA2))


# design sumerisedEperiment object  (dds)

# Rename the 'Cell.population' column to 'cell_types'
colnames(colData(mydata_bulk_RNA2))[colnames(colData(mydata_bulk_RNA2)) == "Cell.population"] <- "cell_types"

#convert to factor
mydata_bulk_RNA2$cell_types <- factor(mydata_bulk_RNA2$cell_types)

#extract celltypes 
levels(colData(mydata_bulk_RNA2)$cell_types)

#build DESeqDataSet
dds_pairwise <- DESeqDataSet(mydata_bulk_RNA2, 
                             design = ~ Patient + cell_types)



# drop PBMCs 
# Filter so only A and B remain
dds_pairwise <- dds_pairwise[, colData(dds_pairwise)$cell_types != "PBMC"]

#drop unused levels 
colData(dds_pairwise) <- droplevels(colData(dds_pairwise))
levels(colData(dds_pairwise)$cell_types)


#filter out total count less 10 per row

dds_pairwise <- dds_pairwise[rowSums(counts(dds_pairwise)) >= 10, ]


# remove the NAs from the data 
# extract the rownames (genes) from dds_pairwise
ensembl_ids <- rownames(dds_pairwise)

# clean IDs (strip version suffix)
ensembl_ids_clean <- sub("\\.\\d+$", "", ensembl_ids)

# query biomaRt for annotation
# connect to the database (bioMart) via the R package (biomaRt)

mart_obj <- useEnsembl(biomart = "genes",
                       dataset = "hsapiens_gene_ensembl",
                       version = 109  )

annot <- getBM(attributes = c("ensembl_gene_id", 
                              "hgnc_symbol", 
                              "gene_biotype"),
               filters = "ensembl_gene_id",
               values = ensembl_ids_clean,
               mart = mart_obj)


# process annotation
annot_clean <- annot %>%
  filter(gene_biotype == "protein_coding",
         hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# create named vectors for lookup
symbol_lookup <- setNames(annot_clean$hgnc_symbol, annot_clean$ensembl_gene_id)
biotype_lookup <- setNames(annot_clean$gene_biotype, annot_clean$ensembl_gene_id)

# create keep vector
matched_genes <- ensembl_ids_clean %in% names(symbol_lookup)

# subset
dds_obj <- dds_pairwise[matched_genes, ]

# add annotations
rowData(dds_obj)$hgnc_symbol <- symbol_lookup[ensembl_ids_clean[matched_genes]]
rowData(dds_obj)$gene_biotype <- biotype_lookup[ensembl_ids_clean[matched_genes]]
rowData(dds_obj)

##run DESeq 
dds_pairwise <- DESeq(dds_obj)

#check comparisons (coeffeicients)
resultsNames(dds_pairwise)


#retrieve the specific levels results

# positive value means high expression in cDC1 and negative means high expression in B cells
#colData(dds_pairwise)$cell_types
res_B_cDC1 <- results(dds_pairwise, contrast = c("cell_types", "B", "cDC1"))

# add symbols 
res_B_cDC1$symbol <- rowData(dds_pairwise)[ rownames(res_B_cDC1), "symbol" ]
res_B_cDC1$genes <- rownames(res_B_cDC1)



#Now convert to df 
res_B_cDC1_df <- as.data.frame(res_B_cDC1)
#remove NAs
res_B_cDC1_df <- res_B_cDC1_df[!is.na(res_B_cDC1_df$padj), ]

res_B_cDC1_df


# create a ggMAplot for B cells against cDC1 cells 

ggplot(res_B_cDC1_df, aes(x = log10(baseMean), y = log2FoldChange)) +
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




