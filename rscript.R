
##### package load #########################################################################################################
{library(DESeq2)
 library(biomaRt)
  
library(tidyverse)
library(factoextra)
library(lolliplot)
library(gridExtra) 
library(grid)
  library(ggrepel)
  library(gt)
  library(scales)
  
  
  
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


# save the MA plot as a sample 
ggsave(
  filename = "plots/b_cdc1_ma_plot.pdf",
  plot = b_cdc1_ma_plot,
  width = 15,
  height = 10,
  units = "in"
)


## volcano plot 
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
  filter(padj < 0.01) %>%                                     # first, keep only significant genes
  arrange(desc(abs(log2FoldChange)), padj) %>%                # then sort
  head(25) %>%                                                 # then take the top 25
  pull(symbol)                                                # get only the symbol


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


b_cdc1_volcano_plot <- ggplot(res_B_cDC1_noNA_df, aes(log2FoldChange, -log10(padj))) +
  
  geom_point(aes(color = plot_group), alpha = 0.75, size = 3.5) +
  
  scale_color_manual(values = c(high = "red",
                                low = "blue",
                                NotSig = "gray",
                                "Top DEG" = "darkgreen"
  )) +
  
  # add thresholds
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.6) +
  
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
  #add theme 
  theme_classic(base_size = 14) +
  
  theme(plot.title = element_text(size = 25,
                                  face = "bold",
                                  hjust = 0.5
  ),
  #axis titles
  axis.title.x = element_text(size = 20,
                              face = "bold"),
  axis.title.y = element_text(size = 20, 
                              face = "bold"),
  
  # Axis tick labels
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
  
  # Legend
  legend.title = element_text(size = 18, face = "bold"),
  legend.text = element_text(size = 15),
  
  # Legend position
  legend.position = "right"
  )


# save the volcano plot as a sample 
ggsave(
  filename = "plots/b_cdc1_volcano_plot.pdf",
  plot = b_cdc1_volcano_plot,
  width = 15,
  height = 10,
  units = "in"
)


# create functions for repeated codes ###################################

# This function would extraact pair wise analyses 

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


# create a function retrive cell of interest 
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





generate_heatmap <- function(heatmap_data,
                             main_title = "Top 25 cDC1 cells markers",
                             gene_symbols = NULL) {
  
  # Prepare annotation – ensure it's a data frame with one column
  # If heatmap_data$annotation_df$cell_types is already a data frame,
  # we convert it to a one‑column data frame with the desired name.
  # If it's a vector, we wrap it.
  if (is.data.frame(heatmap_data$annotation_df$cell_types)) {
    annotation_df <- heatmap_data$annotation_df$cell_types
    colnames(annotation_df) <- "Cell Types"
  } else {
    annotation_df <- data.frame(`Cell Types` = heatmap_data$annotation_df$cell_types,
                                stringsAsFactors = FALSE)
  }
  
  # Create column annotation
  celltype_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = annotation_df,
    annotation_name_side = "right"
  )
  
  # Row labels
  if (!is.null(gene_symbols)) {
    row_labels <- gene_symbols
  } else {
    row_labels <- rownames(heatmap_data$scaled_matrix)
  }
  
  # Create the heatmap
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
    col = circlize::colorRamp2(
      seq(-2, 2, length.out = 100), 
      colorRampPalette(c("navy", "white", "firebrick3"))(100)
    ),
    column_title = paste(main_title, "-", heatmap_data$transformation),
    name = "expression",
    border = TRUE,
    heatmap_legend_param = list(
      title = "expression",
      title_position = "topcenter",
      legend_direction = "horizontal",
      title_gp = grid::gpar(fontsize = 10, fontface = "bold")
    )
  )
}




save_heatmap <- function(
    plot_object,
    file_path = ".",
    filename = "heatmap_plot",
    width = 10,
    height = 6,
    dpi = 300
) {
  # combine file path and filename
  full_path <- file.path(file_path, filename)
  
  # add .pdf extension if not already present
  if (!grepl("\\.pdf$", full_path, ignore.case = TRUE)) {
    full_path <- paste0(full_path, ".pdf")
  }
  
  # save as PDF
  grDevices::pdf(full_path, width = width, height = height)
  print(plot_object)
  dev.off()
  
  message("Heatmap saved as: ", full_path)
  return(invisible(plot_object)) # return but don't display in the console 
}




######functions call ##################################################
# pairwise call

pairwise_output <- get_pairwise_DEGs(
  dds_object = dds_pairwise,
  cell_type_col = "cell_types",
  log2FC_threshold = 0
)

names(pairwise_output)
pairwise_output$summary


# Arrange alphabetically by comparison name
summary_tbl <- pairwise_output$summary %>%
  arrange(comparison)



summary_tbl <- summary_tbl %>%
  rename(
    Comparison = comparison,
    `Up-regulated in Cell1` = n_higher_in_cell1,
    `Up-regulated in Cell2` = n_higher_in_cell2,
    `Total DEGs` = total_DEGs
  )



#  create compact table with three gradient colors and bold Comparison column
# create compact publication-style table with three gradient colors
gt_tbl_pub <- summary_tbl %>%
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
  
  # Row striping
  opt_row_striping()

# display table
gt_tbl_pub

# display the table
gt_tbl_pub



#save the table of differential expressions

gtsave(
  gt_tbl_pub,
  "plots/table_up_down_regulated_genes.png"
)



# get all genes higher in each cell type across all comparisons
names(pairwise_output)
names(pairwise_output$results)
names(pairwise_output$summary)

B_up_df <- get_higher_in_cell(pairwise_output, cell_type = "B")
# preview
head(B_up_df)



cDC1_up_df <- get_higher_in_cell(pairwise_output, cell_type = "cDC1")
cDC2_up_df <- get_higher_in_cell(pairwise_output, cell_type = "cDC2")
DC3_up_df <- get_higher_in_cell(pairwise_output, cell_type = "DC3")
pDC_up_df <- get_higher_in_cell(pairwise_output, cell_type = "pDC")
tDC_up_df <- get_higher_in_cell(pairwise_output, cell_type = "tDC")
mo1_up_df <- get_higher_in_cell(pairwise_output, cell_type = "mo1")
mo2_up_df <- get_higher_in_cell(pairwise_output, cell_type = "mo2")
moint_up_df <- get_higher_in_cell(pairwise_output, cell_type = "moint")




#check the genes present in all the files 
mo2_up_df %>%
  group_by(gene_id) %>%
  summarise(n_comparisons = n_distinct(comparison)) %>%
  arrange(desc(n_comparisons))

top25_B_genes <- top_genes_fn(df_hd_data = B_up_df , top_n = 25, min_comparisons = 8) 
top25_cDC1_genes <- top_genes_fn(df_hd_data = cDC1_up_df , top_n = 25, min_comparisons = 8) 
top25_cDC2_genes <- top_genes_fn(df_hd_data = cDC2_up_df , top_n = 25, min_comparisons = 7) 
top25_DC3_genes <- top_genes_fn(df_hd_data = DC3_up_df , top_n = 25, min_comparisons = 7) 

top25_pDC_genes <- top_genes_fn(df_hd_data = pDC_up_df , top_n = 25, min_comparisons = 8) 

top25_tDC_genes <- top_genes_fn(df_hd_data = tDC_up_df , top_n = 25, min_comparisons = 8) 
top25_mo1_genes <- top_genes_fn(df_hd_data = mo1_up_df , top_n = 25, min_comparisons = 8) 
top25_mo2_genes <- top_genes_fn(df_hd_data = mo2_up_df , top_n = 25, min_comparisons = 7) 
top25_moint_genes <- top_genes_fn(df_hd_data = moint_up_df , top_n = 25, min_comparisons = 7)




# the top 25 genes for each cell type 
top25_B_genes$symbol

top25_cDC1_genes$symbol
top25_cDC2_genes$symbol
top25_DC3_genes$symbol

top25_pDC_genes$symbol
top25_tDC_genes$symbol

top25_mo1_genes$symbol
top25_mo2_genes$symbol
top25_moint_genes$symbol

#find common between cell types 

intersect(top25_moint_genes$symbol, top25_mo2_genes$symbol)



# create a table of top 3 genes for each cell type 

cell_dfs <- list(
  "cDC1" = top25_cDC1_genes,
  "cDC2" = top25_cDC2_genes,
  "DC3"  = top25_DC3_genes,
  "Bcell" = top25_B_genes,
  "mo1"  = top25_mo1_genes,
  "mo2"  = top25_mo2_genes,
  "moint"  = top25_moint_genes,
  "pDC" = top25_pDC_genes,
  "tDC" = top25_tDC_genes
)



# extract top 3, keep stats, and combine
top3_stats_summary <- lapply(names(cell_dfs), function(name) {
  cell_dfs[[name]] %>%
    slice_head(n = 3) %>%
    mutate(cell_type = name) %>%
    # Selecting and renaming for clarity if necessary
    dplyr::select(cell_type, symbol, log2FoldChange, padj)
}) %>%
  bind_rows()

# print the table
# Formatting the table for display
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

# View the result
print(formatted_top3)


#save the table of differential expressions 
ggsave(
  filename = "plots/formatted_top3.pdf",
  plot = formatted_top3,
  width = 15,
  height = 30,
  units = "in"
)


