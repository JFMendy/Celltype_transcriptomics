
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
demo_table

#  generate a vector called totals  to store the row totals 
totals <- demo_table %>%
  dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
  dplyr::mutate(Patient = "Total") %>%
  dplyr::select(Patient, dplyr::everything())
totals

# join the row totals to the demo_table  
demo_table <- dplyr::bind_rows(demo_table, totals)
demo_table 

#generate a more polish demogragraphic table 


#  convert standard table to a grob for plotting
#load the library(gridExtra) library(grid)

sample_demo_table <- tableGrob(demo_table, rows = NULL, theme = ttheme_default(
  core = list(fg_params = list(cex = 0.9)),
  colhead = list(fg_params = list(cex = 1, fontface = "bold"))
))


# display table on a blank plot
grid.draw(sample_demo_table)


#save the plot 
# wrap the plot in grid_arrange for ggplot to open it easily. 
ggsave("plots/sample_demo_table.pdf", 
       plot = gridExtra::grid.arrange(sample_demo_table), 
       width = 8, 
       height = 6)



