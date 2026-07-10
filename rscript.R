
##### package load #########################################################################################################
library(DESeq2)
library(tidyverse)
library(factoextra)
library(lolliplot)
library(gridExtra) 
library(grid)

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


# Task 3, Principle component analysis

