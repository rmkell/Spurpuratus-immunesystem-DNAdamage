#Analysis of timecourse bulk-RNAseq data comparing control and UVB treated (1000 mJ/cm2) S.purpuratus coelomocytes at 0, 1, 3, 6, and 24 hours recovery in vitro. 

#pOverA data filtering method is exemplified here: https://github.com/echille/Mcapitata_OA_Developmental_Gene_Expression_Timeseries/blob/main/4-Differential-Gene-Expression-Analysis/pHTreatment_RNAseqDE.Rmd

#Updated 2/26/26 RMK

#Import necessary libraries
library("tidyverse")
library("genefilter")
library("DESeq2")
library("RColorBrewer")
library("WGCNA")
library("flashClust")
library("gridExtra")
library("ComplexHeatmap")
library("goseq")
library("dplyr")
library("clusterProfiler")
library("simplifyEnrichment")
library(Rmisc)
library(factoextra)
#install.packages("extrafont")
library(extrafont)
#install.packages("ggsignif")
library(ggsignif)

##Data input, cleaning, and pre-processing

# Read in the combined sample, raw count data 
gcount <- as.data.frame(read.csv("RAW_COUNTS.csv", header=TRUE, row.names="Geneid"))
head(gcount)
dim(gcount) #29,213 genes x 29 samples
summary(gcount) #quick way to see if any have NAN values

#Read in the treatment info file: a metadata .csv
#Columns are Sample ID (matching those names in gcount), time_point (0hr, 1hr, 3hr, 6hr, or 24hr), and treatment (control or zap)
metadata <- read.csv("RNAseq_treatmentinfo.csv", header = TRUE, sep = ",")
#This is a .csv file with 3 columns: sample_id, time_point, and treatment.
#sample_id lists every sample name (Control_0hr_A, Control_0hr_B, etc. There is no Zap_3hr_A due to sequencing failure).
#Time_point is 0hr, 1hr, 3hr, 6hr, or 24hr
# Create the new column 'time_point_treatment'
metadata <- metadata %>%
  mutate(time_point_treatment = paste(time_point, treatment, sep = "_"))
#Treatment is either control or zap
head(metadata)

##Quality-filter gene counts
### Check that there are no genes with 0 counts across all samples: if there are, remove those genes

nrow(gcount) #29,213
gcount<-gcount %>%
  mutate(Total = rowSums(.[, 1:29]))%>% #adds a new column Total, which is the sum of values in the first 29 columns for each row
  filter(!Total==0)%>% #filters out rows where total=0
  dplyr::select(!Total) #removes the total column from the dataframe
nrow(gcount) #29,213
#No change in this case, since rows with all zeros were filtered out previously, in the raw data creation step.

### Conduct data filtering

# The following info taken from Ariana Huffmyer and Erin Chille scripts: 
#   
# *pOverA*: Specifying the minimum count for a proportion of samples for each gene. Here, we are using a pOverA of 0.10. This is because we have 29 samples with a minimum of n=3 samples per bleached per date. Therefore, we will accept genes that are present in 3/29 = 0.10 of the samples because we expect different expression by group. 
# 
# We are further setting the minimum count of genes to 10, such that 10% of the samples must have a gene count of >10 in order for the gene to remain in the data set.  
# 
# Filter in the package "genefilter". Pre-filtering our dataset to reduce the memory size dataframe, increase the speed of the transformation and testing functions, and improve quality of statistical analysis by removing low-coverage counts. Removed counts could represent outliers in the data and removing these improves sensitivity of statistical tests


filt <- filterfun(pOverA(0.10,10))

#create filter for the counts data
gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#How many rows do we have before and after filtering?
nrow(gcount) #Before = 29,213
nrow(gcount_filt) #After = 18,032

# Before filtering, we had 29,213 genes. After filtering for pOverA, we have approximately 18,032 genes. This indicates that there were many genes present in <10% of samples at <10 counts per gene.  
# 
# In order for the DESeq2 algorithms to work, the SampleIDs on the metadata file and count matrices have to match exactly and in the same order. The following R clump will check to make sure that these match. Should return TRUE. 


# Display current order of metadata and gene count matrix.  
metadata$sample_id
colnames(gcount_filt)

#Order metadata the same as the column order in the gene matrix.  
list<-colnames(gcount_filt)
list<-as.factor(list)
list

metadata$sample_id<-as.factor(metadata$sample_id)

# Re-order the levels
metadata$sample_id <- factor(as.character(metadata$sample_id), levels=list)
metadata$sample_id

# Re-order the data.frame
metadata_ordered <- metadata[order(metadata$sample_id),]
metadata_ordered$sample_id 

#Checking that all row and column names match. Should return "TRUE"
all(rownames(metadata$sample_id) %in% colnames(gcount_filt))
all(rownames(metadata$sample_id) == colnames(gcount_filt)) 


### Plot the sequence counts after trimming in each sample to confirm that sequence depth is not different between groups

#Visualize distribution of total read counts across different groups to check for sig differences in sequence depth
# Calculate total read counts per sample
total_counts <- colSums(gcount_filt)

# Create a data frame for visualization
sample_ids <- colnames(gcount_filt)
group_info <- metadata$time_point_treatment  # Adjust this to match your metadata structure
visualization_data <- data.frame(sample_id = sample_ids, total_counts = total_counts, group = group_info)

# Create boxplot
ggplot(visualization_data, aes(x = group, y = total_counts)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Total Read Counts per Sample by Group", x = "Group", y = "Total Read Counts")

#Use the compare_means function to perform pairwise comparisons between the control and zap groups at each timepoint.

library(ggpubr)

# Perform pairwise comparisons
comparisons <- list(
  c("0hr_control", "0hr_zap"),
  c("1hr_control", "1hr_zap"),
  c("3hr_control", "3hr_zap"),
  c("6hr_control", "6hr_zap"),
  c("24hr_control", "24hr_zap")
)

# Create boxplot with actual p-values
ggplot(visualization_data, aes(x = group, y = total_counts)) +
  geom_boxplot() +
  theme_classic() +
  labs(title = "Total Read Counts per Sample by Group", x = "Group", y = "Total Read Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = comparisons, method = "t.test", label = "p.signif") #change label to "p.format" to see actual p values

## Construct DESeq2 data set

Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at groups to test for any differences in gene expression across Time and Treatment.

#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = metadata_ordered,
                               design = ~time_point*treatment)
#NOTE regarding the design: ~time_point*treatment model includes the main effects of time_point and treatment, as well as their interaction.
#The interaction term allows the effect of treatment to vary depending on the time_point and vice versa. Using this bc we suspect the treatment effect may differ at different timepoints
#This model is more flexible and can capture differences in the treatment effect at different time_points.
#If you wanted to test the main effects of treatment and time_point separately, you would use design = ~treatment + time_point

metadata_ordered <- metadata_ordered %>% unite(group, treatment, time_point, remove = FALSE) #combines the treatment and time_point columns into a new column called group. remove=False keeps the original time_point and treatment columns in the dataframe.

#Set alternate DESeq2 design using group
gdds2 <- DESeqDataSetFromMatrix(countData = gcount_filt,
                                colData = metadata_ordered,
                                design = ~group)

# First we are going to log-transform the data using a variance stabilizing transforamtion (VST). This is only for visualization purposes. Essentially, this is roughly similar to putting the data on the log2 scale. It will deal with the sampling variability of low counts by calculating within-group variability (if blind=FALSE). Importantly, it does not use the design to remove variation in the data, and so can be used to examine if there may be any variability do to technical factors such as extraction batch effects.
# 
# To do this we first need to calculate the size factors of our samples. This is a rough estimate of how many reads each sample contains compared to the others. In order to use VST (the faster log2 transforming process) to log-transform our data, the size factors need to be less than 4.
# 
# Chunk should return TRUE if <4.  

### first design (gdds)
SF.gdds <- estimateSizeFactors(gdds) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 for us to use vst

print(sizeFactors(SF.gdds)) #View size factors
all(sizeFactors(SF.gdds)) < 4 #returned TRUE

### second design (gdds2)
SF.gdds2 <- estimateSizeFactors(gdds2) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than 4 for us to use vst

print(sizeFactors(SF.gdds2)) #View size factors
all(sizeFactors(SF.gdds2)) < 4 #returned TRUE

# The size factors of gdds and gdds2 are exactly the same. The size factors are determined solely by the count data and are independent of the design formula used for differential expression analysis. This is why the size factors for gdds and gdds2 are the same, even though the design formulas are different. The design formula affects the modeling of the data for differential expression testing, but not the normalization step.
# All size factors are less than 4, so we can use VST transformation.  

## using the treatment x time design
gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transformation to minimize effects of small counts and normalize wrt library size
head(assay(gvst), 3) #view transformed gene count data for the first three genes in the dataset.

## using the group design
gvst2 <- vst(gdds2, blind=FALSE) #apply a variance stabilizing transformation to minimize effects of small counts and normalize wrt library size
head(assay(gvst2), 3) #view transformed gene count data for the first three genes in the dataset.


#NOTE: design ~group tells DESeq2 to model gene expression based on the group variable, which has 10 levels:
# 
# control_0hr
# zap_0hr
# control_1hr
# zap_1hr
# control_3hr
# zap_3hr
# control_6hr
# zap_6hr
# control_24hr
# zap_24hr
# 
# This is a categorical variable with multiple levels, so DESeq2 will treat one of them as the reference level (usually the first alphabetically or the first factor level), and compare all others to it.
# 
# DESeq2 by default compares the first non-reference level to the reference level of your group factor. For example, if control_0hr is the reference, then:
#   
#   log2FoldChange in res will show the log2 fold change of the first contrast, e.g., zap_0hr vs control_0hr.
# 
# If you want to compare specific groups, use the contrast argument in results():
#   
#   results(DEG2, contrast = c("group", "zap_1hr", "control_1hr"))


#### Differential Gene Expression Analysis
#Run differential expression test using a Wald model. 


#Group model
This design treats each combination of time_point and treatment as a unique group. It simplifies the model by comparing each group directly to the reference group (ctrl0hr).

DEG2 <-DESeq(gdds2)

#Look at the raw counts 
head(DEG2@assays@data$counts)

#Look at the normalized counts 
normalized_counts <- counts(DEG2, normalized = T)
head(normalized_counts)
#These normalized counts can be used to compare gene expression values between samples.



#Run DESeq and Extract results. Note that you have to explicitly call what timepoints you are comparing. If you dont, it will just compare the first "group" to the second one. See note above.

# Load required libraries
library(dplyr)

# Define time points
time_points <- c("0hr", "1hr", "3hr", "6hr", "24hr")

# Load annotation lookup table
lookup <- read.csv("updated_SP_lookuptable_from_GTF.csv", header = TRUE)

# Create a list to store filtered results
filtered_results <- list()

# Loop through each time point
for (tp in time_points) {
  control_group <- paste0("control_", tp)
  zap_group <- paste0("zap_", tp)
  
  # Get DE results for zap vs control at this time point
  res <- results(DEG2, contrast = c("group", zap_group, control_group))
  
  # Add comparison label
  res$comparison <- paste(zap_group, "vs", control_group)
  
  # Convert to data frame and add gene_id column
  res_df <- as.data.frame(res) %>%
    mutate(gene_id = rownames(res))
  
  # Save unfiltered results
  write.csv(res_df, file = paste0("DEG_", zap_group, "_vs_", control_group, "_ALL.csv"), row.names = FALSE)
  
  # Filter for significant DEGs
  res_filtered <- res_df %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1)
  
  # Count up/downregulated genes
  upregulated <- sum(res_filtered$log2FoldChange > 1, na.rm = TRUE)
  downregulated <- sum(res_filtered$log2FoldChange < -1, na.rm = TRUE)
  total_DEGs <- upregulated + downregulated
  
  # Print summary
  cat("\nComparison:", zap_group, "vs", control_group, "\n")
  cat("  Upregulated genes:   ", upregulated, "\n")
  cat("  Downregulated genes: ", downregulated, "\n")
  cat("  Total DEGs:          ", total_DEGs, "\n")
  
  # Annotate with lookup table
  res_annotated <- merge(res_filtered, lookup, by = "gene_id", all.x = TRUE)
  
  # Save filtered annotated results
  write.csv(res_annotated, file = paste0("DEG_", zap_group, "_vs_", control_group, "_FILTERED_ANNOTATED.csv"), row.names = FALSE)
  
  # Store in list
  filtered_results[[tp]] <- res_annotated
}

# Combine all filtered annotated results
combined_filtered <- do.call(rbind, filtered_results)
#Note: gene_ID column will get timepoint added to it, because row names cannot be duplicated. 
#To remove, just use data > sort> delimited = . in excel

nrow(combined_filtered) #3,813 differentially expressed genes over all timepoints

# Save combined results
write.csv(combined_filtered, file = "DEG_all_filtered_zap_vs_control_ANNOTATED.csv", row.names = FALSE)

#Need to pull out the gene IDs from DEG2 corresponding to differentially expressed genes (lfc >1 or <-1, and padj < 0.05)
sig_ids <- data.frame(gene_id = unique(combined_filtered$gene_id))

#Sort alphabetically
sig_ids <- data.frame(gene_id = sort(unique(combined_filtered$gene_id)))
nrow(sig_ids) #Down to 2,581 unique DEGs

#Ensure sig_ids$gene_id is a character vector
sig_gene_ids <- as.character(sig_ids$gene_id)
head(sig_gene_ids)

# Subset DEG2 to only include rows with gene_id that are in sig_gene_ids. 
DEG2_SigGenes_Subset <- DEG2[rownames(DEG2) %in% sig_gene_ids, ]
nrow(DEG2_SigGenes_Subset) #good, still 2,581

#Apply a variance stabilizing transformation to the subset DEG2 dataset, which contains only the 2,581 DE genes, using vst(). Applying a variance stabilizing transformation minimizes the effects of small counts and normalizes with regard to library size.

VST_DEG2 <- vst(DEG2_SigGenes_Subset, blind = TRUE)
head(assay(VST_DEG2)) #view transformed gene count data
dim(assay(VST_DEG2)) ## should be 2,581 DE genes x 29 samples 


# Distance plot of the DEGs using variance stabilizing transformation values
# The distance plot calculates the Euclidean distance between the expression values for each individual sample. This is a way to identify which samples are most closely related in terms of their gene expression patterns.
plotDists = function (VST_DEG2.obj) {
  sampleDists <- dist(t(assay(VST_DEG2.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste(VST_DEG2.obj$condition )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
}
plotDists(VST_DEG2)


#PCA plot using VST_DEG2 values
#VST_DEG2 is the SUBSETTED DEG2 values (including only those 2,581 DE genes, not all genes)
#A principal components plot is another way to observe how diverse the samples are. As always, the samples which are most similar with each other in terms of their gene expression values will be closer to each other on the plot. This plot uses a Cartesian coordinate system, and the axes that are displayed correspond to the top two principal components that explain the majority of the variability in the data.

plot_PCA = function (VST_DEG2.obj) {
  pcaData <- plotPCA(VST_DEG2.obj,  intgroup = c("time_point_treatment"), returnData = T)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=time_point_treatment)) +
    geom_point(size=3) +
    labs(x = paste0("PC1: ",percentVar[1],"% variance"),
         y = paste0("PC2: ",percentVar[2],"% variance"),
         title = "PCA Plot colored by condition") +
    ggrepel::geom_text_repel(aes(label = name), color = "black")
}
plot_PCA(VST_DEG2)



### Plot a PCA
metadata_ordered$time_point <- as.character(metadata_ordered$time_point)


gPCAdata <- plotPCA(VST_DEG2, intgroup = c("time_point", "treatment"), returnData=TRUE) #use ntop to specify top n genes or all genes; 18,032 for all
gPCAdata$time_point <- as.character(gPCAdata$time_point)
metadata$time_point <- as.character(metadata$time_point)
gPCAdata <- gPCAdata %>% dplyr::rename(sample_id = name) 
gPCAdata <- left_join(gPCAdata, metadata_ordered, by = c("sample_id", "time_point", "treatment"))

percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data

library(scales)
library(ggrepel)

# Convert treatment to a factor
gPCAdata$treatment <- as.factor(gPCAdata$treatment)

# Extract the Urchin identifier and add it as a new column
gPCAdata <- gPCAdata %>%
  mutate(Urchin = sub(".*_(.)$", "\\1", sample_id))

# Check the levels of the treatment variable
print(levels(gPCAdata$treatment))

gPCAdata$Urchin <- as.factor(gPCAdata$Urchin)
print(levels(gPCAdata$Urchin))

# Plot ALL SAMPLES in one PCA plot
# Define the fill colors for the urchins
fill_values <- c("A" = "#ffff00", "B" = "#8000ff", "C" = "#ff69b4") # Adjust as needed

# Define the shapes for controls and zaps
shape_values <- c("control" = 21, "zap" = 24) # Adjust the names to match your actual treatment levels

allgenesPCA <- gPCAdata %>% 
  unite(Group, time_point, remove = FALSE) %>%
  ggplot(aes(PC1, PC2, color=time_point, shape=treatment, fill=Urchin, group=treatment)) + 
  geom_point(size=3, stroke=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c("#ea5f94", "#006bff", "#ff7f0e", "#2ca02c", "#d62728")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_fill_manual(values = fill_values, guide = guide_legend(override.aes = list(shape = 21))) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background=element_blank()) +
  theme(legend.text = element_text(size=18), 
        legend.position="right",
        legend.title = element_text(size=20), 
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size=18), 
        title = element_text(size=25, face="bold"), 
        axis.title = element_text(size=18))

print(allgenesPCA)
ggsave(filename="PCA_allDEGgenes_allsamples.jpeg", 
       plot=allgenesPCA, dpi=300, width=8, height=7, units="in")

# Plot PCA split by TIME

# Convert treatment and Urchin to factors if not already
gPCAdata$treatment <- as.factor(gPCAdata$treatment)
gPCAdata$Urchin <- as.factor(gPCAdata$Urchin)

# Reorder the time_point factor levels
gPCAdata$time_point <- factor(gPCAdata$time_point, levels = c("0hr", "1hr", "3hr", "6hr", "24hr"))

# Check the levels of the treatment and Urchin variables
print(levels(gPCAdata$treatment))
print(levels(gPCAdata$Urchin))
print(levels(gPCAdata$time_point))

# Define the fill colors for the urchins
fill_values <- c("A" = "#ffff00", "B" = "#8000ff", "C" = "#ff69b4") # Adjust as needed

# Define the shapes for controls and zaps
shape_values <- c("control" = 21, "zap" = 24) # Adjust the names to match your actual treatment levels

allgenesPCA_tp_split <- gPCAdata %>% 
  unite(Group, time_point, remove = FALSE) %>%
  ggplot(aes(PC1, PC2, shape=treatment, fill=Urchin, group=treatment)) + 
  geom_point(size=3, stroke=1) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_shape_manual(values = shape_values) +
  scale_fill_manual(values = fill_values, guide = guide_legend(override.aes = list(shape = 21))) +
  theme_bw() + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border around each plot
        axis.line = element_line(colour = "black"),
        plot.background=element_blank(),
        strip.text = element_text(size = 20, face = "bold")) + # Increase font size of plot titles
  theme(legend.text = element_text(size=18), 
        legend.position="right",
        legend.title = element_text(size=20), 
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text = element_text(size=18), 
        title = element_text(size=25, face="bold"), 
        axis.title = element_text(size=18)) +
  facet_wrap(~time_point, ncol = 5)

print(allgenesPCA_tp_split)

ggsave(filename="PCA_allgenes_by_time_point.jpeg", 
       plot=allgenesPCA, dpi=300, width=25, height=5, units="in")


##### Create matrix and cluster the 2,581 DEGs using NbClust ####

mat_DEG2 <- assay(VST_DEG2) #Extracts the expression data from the DEG2_vst object
nrow(mat_DEG2) #2,581 genes
mat_DEG2_centered <- mat_DEG2 - rowMeans(mat_DEG2) #Subtracts the mean expression value of each gene from its expression values across all samples. This centers the data by making the mean expression value of each gene zero.
nrow(mat_DEG2_centered)

### Compute optimal number for clustering using input of the 2,821 DEGs stored in DEG2_VST

Find the optimum number of clusters using 30 indexes with the NbClust() package.
NOTE: The NbClust function below only work mat_DEG2, NOT mat_DEG2_centered
```{r}
#install.packages('NbClust')
library(NbClust)

nb <- NbClust(mat_DEG2, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans") #use k-means clustering, consider a minimum of 2 clusters and a max of 10


### Get and save clustering results with DESeq2 results matrix "combined_filtered"
#NOTE: it suggested 3 k means clusters, but this resulted in cluster2 having only 61 genes. Therefore using 2 clusters only.
calc_kmeans <- kmeans(mat_DEG2_centered, 2)
head(calc_kmeans)

#Create a dataframe of the clustering results
cluster_res <- data.frame(gene_id = names(calc_kmeans$cluster), cluster = calc_kmeans$cluster)

#Append the clustering results to the "combined_filtered" dataframe
combined_filtered_clusters_added <- merge(combined_filtered, cluster_res, by = "gene_id")
#Save to .csv
write.csv(combined_filtered_clusters_added, file="FINAL_DEG_LIST_with_annotations_and_clusters.csv")


# Plot a heatmap
#NOTE: I made cluster assignments above using calc_kmeans <- kmeans(mat_DEG2_centered, 2). Do NOT let ComplexHeatmap compute its own clusters, as you will get different results. 


# First, make sure the cluster assignments are in the same order as the rows in mat_DEG2_centered
gene_clusters <- cluster_res[match(rownames(mat_DEG2_centered), cluster_res$gene_id), "cluster"]


# Convert to a factor with desired cluster order and labels
gene_clusters <- factor(gene_clusters, levels = c(1, 2), labels = c("Cluster1", "Cluster2"))
row_title = levels(gene_clusters)


DEG2_clust <- subset(combined_filtered_clusters_added, select = c(gene_id, cluster)) # creates a new data frame DEG2_clust containing only the gene_id and cluster columns from DEG2_list.
nrow(DEG2_clust) #3,813, since not filtered to uniqueness yet

#Creates a newdataframe, hm_ann_row, that contains only the unique gene ids
hm_ann_row <- unique(DEG2_clust) #remove duplicate rows
nrow(hm_ann_row) #now correct, 2,581 unique DEGs

rownames(hm_ann_row) <- hm_ann_row$gene_id #set row names to the gene_id values
hm_ann_row <- subset(hm_ann_row, select=cluster) #keep only the cluster column

#replaces numeric cluster labels with more descriptive names (Cluster1, Cluster2).
hm_ann_row$cluster <- gsub(1,"Cluster1",hm_ann_row$cluster)  
hm_ann_row$cluster <- gsub(2,"Cluster2",hm_ann_row$cluster)


hm_ann_row <- as.matrix(hm_ann_row[rownames(mat_DEG2_centered),]) #converts hm_ann_row to a matrix, ensuring the rows match those in mat_DEG2.

hmGroup <- colData(VST_DEG2)[c("treatment", "time_point")] #Extracts the treatment and time_point columns from VST_DEG2.
hmGroup$time_point <- as.character(hmGroup$time_point)  #Converts the time_point column to character type.
# Ensure correct order of timepoints
hmGroup$time_point <- factor(hmGroup$time_point, levels = c("0hr", "1hr", "3hr", "6hr", "24hr"))

Treatmentcol = c("control"="#006bff", "zap" ="#ff0000") #define color mappings for the treatment and time_point annotations.
Timepointcol = c(
  "0hr" = "#e5f5e0",  # Light green
  "1hr" = "#a1d99b",  # Light-medium green
  "3hr" = "#74c476",  # Medium green
  "6hr" = "#31a354",  # Dark-medium green
  "24hr" = "#006d2c"  # Dark green
)

hm_ann_col <- HeatmapAnnotation(df=hmGroup, #gp = gpar(col = "grey"),
                                col = list(treatment=Treatmentcol, time_point=Timepointcol)) # creates a HeatmapAnnotation object using the hmGroup data frame and the defined color mappings.

nrow(mat_DEG2_centered)
row_title = levels(gene_clusters)
gene_clusters <- factor(gene_clusters, levels = c("Cluster1", "Cluster2"))
levels(gene_clusters)


DEGheatmap <-  Heatmap(mat_DEG2_centered, column_title = "", 
                       name = "expression",
                       show_row_names = FALSE, 
                       top_annotation = hm_ann_col, show_column_names = TRUE, row_dend_side = "left" ,
                       column_split = 3, column_dend_height = unit(0.5, "in"),
                       cluster_row_slices = FALSE, #disable reordering, ensuring cluster 1 is on top
                       row_split = gene_clusters,
                       row_km_repeats = 100, 
                       use_raster = TRUE,
                       #row_title = levels(gene_clusters),
                       row_gap = unit(2.5, "mm"), border = TRUE,
                       column_names_gp =  gpar(fontsize = 10))
#Save to PDF
pdf(file = "DEGheatmap_centered_NEW_2851DEGs.pdf")
print(DEGheatmap)
dev.off()

#How many genes did Heatmap put into each of the n=2 k-means clusters??
# Extract the row clustering (k-means) information
row_clusters <- row_order(DEGheatmap)

# Count how many genes are in each cluster
sapply(row_clusters, length)

Cluster2 Cluster1 
1334     1247

### List of clusters 
cluster1 <- cluster_res %>% filter(cluster == 1) 
cluster2 <- cluster_res %>% filter(cluster == 2) 

nrow(cluster1) #1,247
nrow(cluster2) #1,334

total = nrow(cluster1) + nrow(cluster2)
print(total) #2,581

#Ok NOW the heatmap matches the clusters assigned by kmeans()

#Save lists of gene IDs assigned to each cluster
write.csv(cluster1, file = 'Cluster1_ID_list.csv')
write.csv(cluster2, file = 'Cluster2_ID_list.csv')


Cluster1_ID_list = read.csv('Cluster1_ID_list.csv', header = TRUE)
Cluster1_ID_list

#Annotate the clusterID lists
# Load lookup table (.csv)
lookup = read.csv('updated_SP_lookuptable_from_GTF.csv', header = TRUE)
lookup

#Merge DEG2_list_with_clusters with lookup
Cluster1_ID_list_annotated <- merge(Cluster1_ID_list, lookup, by = "gene_id", all.x = TRUE)

# View the merged dataframe
head(Cluster1_ID_list_annotated)
nrow(Cluster1_ID_list_annotated)

#Save DEG2_list_with_clusters_genedescriptions
write.csv(Cluster1_ID_list_annotated, file="Cluster1_ID_list_annotated.csv")

Cluster2_ID_list = read.csv('Cluster2_ID_list.csv', header = TRUE)
Cluster2_ID_list
#Annotate the clusterID lists
# Load lookup table (.csv)
lookup = read.csv('updated_SP_lookuptable_from_GTF.csv', header = TRUE)
lookup
##Merge with lookup
Cluster2_ID_list_annotated <- merge(Cluster2_ID_list, lookup, by = "gene_id", all.x = TRUE)
# View the merged dataframe
head(Cluster2_ID_list_annotated)
nrow(Cluster2_ID_list_annotated)
#Save DEG2_list_with_clusters_genedescriptions
write.csv(Cluster2_ID_list_annotated, file="Cluster2_ID_list_annotated.csv")


#Retrieve protein FASTA of cluster 1 and cluster2 sequences
#2_4_25 Retrieve protein sequences FASTA given an input of Spurp LOC gene IDs.
#Biostrings is used to read/handle FASTA files.

library(Biostrings)

################################## CLUSTER 1 ##########################
#Read the protein FASTA into R
fasta_file <- "C:/Users/Riss Kellogg/Desktop/scRNAseq/Spur_5.0.pep.all_LOCheaders.fasta"
protein_sequences <- readAAStringSet(fasta_file)

#Load Your List of LOC IDs: Assuming your list of LOC IDs is in a text file or a CSV, you can read it into R using read.table or read.csv.
Cluster1_loc_ids <- read.csv("Cluster1_ID_list.csv", header = TRUE, stringsAsFactors = FALSE)
Cluster1_loc_ids <- Cluster1_loc_ids$gene_id  # Assuming the LOC IDs are in column named "gene_id"
Cluster1_loc_ids
summary(Cluster1_loc_ids) #should be 1247 total

#Match LOC IDs and Extract Sequences: Match the LOC IDs with the names in your protein_sequences object and extract the corresponding sequences.
Cluster1_matched_sequences <- protein_sequences[names(protein_sequences) %in% Cluster1_loc_ids]
summary(Cluster1_matched_sequences) #1165 matched
#Find unmatched IDs
Cluster1_unmatched_ids <- Cluster1_loc_ids[!Cluster1_loc_ids %in% names(protein_sequences)]
# Print the unmatched IDs
print(Cluster1_unmatched_ids) #Many of them

#Input the unmatched IDs to echinobase, try to find protein sequence, add it to the reference protein FASTA. Run matching/protein seq pull again.

#Save the matched sequences AND unmatched sequences to a new FASTA file
writeXStringSet(Cluster1_matched_sequences, "Cluster1_proteinFASTA.fasta")
write.csv(Cluster1_unmatched_ids, file = 'Cluster1_NCRNAs.csv')

################################## CLUSTER 2 ##########################
#Load Your List of LOC IDs: Assuming your list of LOC IDs is in a text file or a CSV, you can read it into R using read.table or read.csv.
Cluster2_loc_ids <- read.csv("Cluster2_ID_list.csv", header = TRUE, stringsAsFactors = FALSE)
Cluster2_loc_ids <- Cluster2_loc_ids$gene_id  # Assuming the LOC IDs are in column named "gene"
Cluster2_loc_ids
summary(Cluster2_loc_ids) #1334 total

#Match LOC IDs and Extract Sequences: You can then match the LOC IDs with the names in your protein_sequences object and extract the corresponding sequences.
Cluster2_matched_sequences <- protein_sequences[names(protein_sequences) %in% Cluster2_loc_ids]
summary(Cluster2_matched_sequences) #1300 total
#Find unmatched IDs
Cluster2_unmatched_ids <- Cluster2_loc_ids[!Cluster2_loc_ids %in% names(protein_sequences)]
# Print the unmatched IDs
print(Cluster2_unmatched_ids) #many

#Input the unmatched IDs to echinobase, try to find protein sequence, add it to the reference protein FASTA. Run matching/protein seq pull again.
#Confirmed that many unmatched LOC IDs hit to noncoding or pseudo RNA (so no protein seq).

#Save the matched sequences to a new FASTA file
writeXStringSet(Cluster2_matched_sequences, "Cluster2_proteinFASTA.fasta")
write.csv(Cluster2_unmatched_ids, file = 'Cluster2_NCRNAs.csv')



#####Gene count plots of the above (VST counts matrix)
GO_meta <- metadata_ordered %>% dplyr::rename(Sample.ID = sample_id) #renames the sample_id column in the metadata_ordered dataframe to Sample.ID.

GOcluster_exp <- as.data.frame(assay(VST_DEG2)) ### Converts the VST-transformed counts matrix (VST_DEG2) to a dataframe.
GOcluster_exp <- GOcluster_exp %>% rownames_to_column(., var = "gene") #change the row names to a new column called "gene"
GOcluster_exp <- GOcluster_exp %>% gather("Sample.ID", "count", 2:30)
GOcluster_exp <- GOcluster_exp %>% left_join(., GO_meta, by = "Sample.ID")


#These lines filter the VST_DEG expression data to include only genes from cluster1 and cluster2.
GOcluster1_exp <- GOcluster_exp[GOcluster_exp$gene %in% cluster1$gene_id,]
GOcluster2_exp <- GOcluster_exp[GOcluster_exp$gene %in% cluster2$gene_id,] 


#These lines calculate summary statistics (e.g., mean, standard error) for the counts in each cluster, grouped by Sample.ID, treatment, and time_point.
GOcluster1_samplemeans <- summarySE(GOcluster1_exp, measurevar = c("count"), groupvars = c("Sample.ID","treatment","time_point"))
GOcluster2_samplemeans <- summarySE(GOcluster2_exp, measurevar = c("count"), groupvars = c("Sample.ID","treatment","time_point"))

GOcluster1_groupmeans <- summarySE(GOcluster1_samplemeans, measurevar = c("count"), groupvars = c("treatment", "time_point"))
GOcluster2_groupmeans <- summarySE(GOcluster2_samplemeans, measurevar = c("count"), groupvars = c("treatment", "time_point"))


## CLUSTER 1 PLOT ##
#font_import()
#font_import(prompt = FALSE)
loadfonts(device = "win")
library(extrafont)

GOcluster1_expmeans.plot <-
  GOcluster1_samplemeans %>% 
  ggplot(., aes(x=time_point, y=count, color=treatment)) +
  geom_jitter(alpha=0.2, width=0.1) +
  theme_bw() + 
  scale_color_manual(values = c("#006bff", "#ff0000"), labels = c("Control", "UVB-dosed")) +  # Change legend labels
  geom_point(data = GOcluster1_groupmeans, aes(x=time_point, y=count), size=3) +
  geom_errorbar(data = GOcluster1_groupmeans, aes(ymin=count-se, ymax=count+se), size=0.8, width=.1) + #error bars calculated using standard error (SE) of the mean. "summarySE" calculates the mean and standard error for the counts in each cluster, grouped by Sample.ID, treatment, and time_point.
  geom_line(data=GOcluster1_groupmeans, aes(group = treatment), size=0.8) +
  labs(title = "Cluster 1", y = "Counts (VST)", x = "Time point", color = "Treatment") +  # Add plot title, y-axis title, x-axis title, and legend title
  scale_x_discrete(limits = c("0hr", "1hr", "3hr", "6hr", "24hr")) +  # Reorder x-axis
  theme(
    text = element_text(family = "Arial"),  # Set all fonts to Arial
    legend.text = element_text(size = 12),  # Set legend text size to 12
    legend.title = element_text(size = 16),  # Set legend title size to 18
    axis.text = element_text(size = 12, color = "black"),  # Set axis tick mark label size to 12 and color to black
    axis.title = element_text(size = 16),  # Set axis titles to 16
    plot.title = element_text(size = 18, face = "bold"),  # Set chart title to size 18 and bold
    strip.text.x = element_text(size = 12, color = "black", face = "bold"),
    strip.text.y = element_text(size = 12, color = "black", face = "bold"),
    strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid")
  ) +
  annotate("text", x = Inf, y = Inf, label = "time_point: p = 0.391(ns)\ntreatment: p = 6.21e-13(***)\ntime_point x treatment: p = 0.143(ns)", 
           hjust = 1.1, vjust = 1.1, size = 3, color = "black")

# Print the plot
print(GOcluster1_expmeans.plot)



#Statistics on the above

aov <-aov(count ~ time_point*treatment, data=GOcluster1_samplemeans)
summary(aov)

#Copy pasted output to excel sheet to save

# Perform Tukey's HSD test
tukey_result <- TukeyHSD(aov)

# Convert the TukeyHSD result to a data frame
tukey_table <- as.data.frame(tukey_result$`time_point:treatment`)

# Print the table
print(tukey_table)

ggsave(filename="GOcluster1_expmeans.jpeg", 
       plot=GOcluster1_expmeans.plot, dpi=300, width=5.5, height=4.5, units="in")

## CLUSTER 2 PLOT ##

GOcluster2_expmeans.plot <-
  GOcluster2_samplemeans %>% 
  ggplot(., aes(x=time_point, y=count, color=treatment)) +
  geom_jitter(alpha=0.2, width=0.1) +
  theme_bw() + 
  scale_color_manual(values = c("#006bff", "#ff0000"), labels = c("Control", "UVB-dosed")) +  # Change legend labels
  geom_point(data = GOcluster2_groupmeans, aes(x=time_point, y=count), size=3) +
  geom_errorbar(data = GOcluster2_groupmeans, aes(ymin=count-se, ymax=count+se), size=0.8, width=.1) +
  geom_line(data=GOcluster2_groupmeans, aes(group = treatment), size=0.8) +
  labs(title = "Cluster 2", y = "Counts (VST)", x = "Time point", color = "Treatment") +  # Add plot title, y-axis title, x-axis title, and legend title
  scale_x_discrete(limits = c("0hr", "1hr", "3hr", "6hr", "24hr")) +  # Reorder x-axis
  theme(
    text = element_text(family = "Arial"),  # Set all fonts to Arial
    legend.text = element_text(size = 12),  # Set legend text size to 12
    legend.title = element_text(size = 16),  # Set legend title size to 18
    axis.text = element_text(size = 12, color = "black"),  # Set axis tick mark label size to 12 and color to black
    axis.title = element_text(size = 16),  # Set axis titles to 16
    plot.title = element_text(size = 18, face = "bold"),  # Set chart title to size 18 and bold
    strip.text.x = element_text(size = 12, color = "black", face = "bold"),
    strip.text.y = element_text(size = 12, color = "black", face = "bold"),
    strip.background = element_rect(color="black", fill="white", linewidth=0.5, linetype="solid")
  ) +
  annotate("text", x = -Inf, y = Inf, label = "time_point: p = 1.08e-4(***)\ntreatment: p = 1.38e-8(***)\ntime_point x treatment: p = 0.38(*)", 
           hjust = -0.1, vjust = 1.1, size = 3, color = "black")


# Print the plot
print(GOcluster2_expmeans.plot)


#Statistics on the above (replace cluster # for 1,2)

aov2 <-aov(count ~ time_point*treatment, data=GOcluster2_samplemeans)
summary(aov2)

# Perform Tukey's HSD test
tukey_result2 <- TukeyHSD(aov2)
tukey_result2

# Convert the TukeyHSD result to a data frame
tukey_table2 <- as.data.frame(tukey_result2$`time_point:treatment`)

# Print the table
print(tukey_table2)

ggsave(filename="GOcluster2_expmeans.jpeg", 
       plot=GOcluster2_expmeans.plot, dpi=300, width=5.5, height=4.5, units="in")


####### Identify the top 100 genes in both clusters using a centroid approach. Centroid is the av. expression of all genes in that cluster. Genes most representative of the cluster have expression profiles over time that match the cluster's overall pattern.

# Use centered, VST matrix
mat <- mat_DEG2_centered

# Get cluster info into mat
cluster_vector <- calc_kmeans$cluster
names(cluster_vector) <- rownames(mat)

# Get top genes defined by best Pearson correlation to cluster centroid
get_top_genes <- function(cluster_number, n = 10) {
  genes_in_cluster <- names(cluster_vector[cluster_vector == cluster_number])
  cluster_mat <- mat[genes_in_cluster, ]
  
  # Compute centroid (mean expression profile)
  centroid <- colMeans(cluster_mat)
  
  # Compute correlation of each gene to centroid
  cors <- apply(cluster_mat, 1, function(x)
    cor(x, centroid)) #Note: cor() function finds Pearson correlation between value X and the centroid. Pearson is used by default, if you dont specify.
  cors_sorted <- sort(cors, decreasing = TRUE)[1:n]

  # Get dataframe of top genes
  top_genes <- data.frame(
    gene = names(cors_sorted),
    correlation_to_centroid = as.numeric(cors_sorted),
    cluster = paste0("Cluster", cluster_number),
    row.names = NULL
  )
  
  return(top_genes)
}

top10_cluster1 <- get_top_genes(1, 10)
top10_cluster2 <- get_top_genes(2, 10)

top10_cluster1
top10_cluster2

top100_cluster1 <- get_top_genes(1, 100)
top100_cluster2 <- get_top_genes(2, 100)

top100_cluster1
top100_cluster2


# Combine top 10 genes into a data frame
top_genes_df <- rbind(top10_cluster1, top10_cluster2)

# Combine top 100 genes into a data frame
top_100genes_df <- rbind(top100_cluster1, top100_cluster2)

# Save to Excel file
write_xlsx(top_genes_df, "Top10_genes_per_cluster.xlsx")
write_xlsx(top_100genes_df, "Top100_genes_per_cluster.xlsx")

######## Create line plots of top 10 genes each cluster #########

#All gene names are LOC IDs. Map them to more informative short gene names. These were looked up manually and added to dictionary below.
gene_map1 <- data.frame(
  short_name = c(
    "RPL18",
    "RPS14",
    "RPL28",
    "L26",
    "RWDD4",
    "RPL7",
    "RPS18",
    "RPS24",
    "RPL23a",
    "RPS25",
    "DDX17",
    "ARF16",
    "TSPAN4",
    "NACA",
    "AGO1",
    "STK17A",
    "Unchar.",
    "ITIH3",
    "JNK",
    "LLGL2"
  ),
  gene = c(
    "LOC373246",
    "rps14",
    "LOC754209",
    "LOC585872",
    "LOC100893494",
    "LOC755731",
    "LOC105436394",
    "LOC592738",
    "LOC588157",
    "LOC589961",
    "LOC578813",
    "LOC115918858",
    "LOC593238",
    "LOC100889842",
    "ago1",
    "LOC582485",
    "LOC100889516",
    "LOC584806",
    "LOC580919",
    "LOC578619"
  ),
  stringsAsFactors = FALSE
)

#Combine top10 genes into 1 list
top_genes_all <- c(top10_cluster1, top10_cluster2)

#Extract expression data
top_gene_exp <- GOcluster_exp %>%
  filter(gene %in% top_genes_all)

#Order by timepoint
top_gene_exp$time_point <- factor(top_gene_exp$time_point,
                                  levels = c("0hr","1hr","3hr","6hr","24hr"))
# Summarize mean +/- SE
top_gene_summary <- top_gene_exp %>%
  group_by(gene, treatment, time_point) %>%
  dplyr::summarise(
    mean_expr = mean(count),
    se = sd(count)/sqrt(length(count)),
    .groups = "drop"
  )

#For plotting, want all cluster 1 genes to appear first, therefore order plot facets by cluster
gene_order <- c(top10_cluster1, top10_cluster2)

top_gene_summary$gene <- factor(top_gene_summary$gene,
                                levels = gene_order)

#PLOT
top_gene_plot <- ggplot(top_gene_summary,
                        aes(x=time_point,
                            y=mean_expr,
                            color=treatment,
                            group=treatment)) +
  geom_line(size=0.8) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_expr-se,
                    ymax=mean_expr+se),
                width=0.1) +
  facet_wrap(~gene, scales="free_y", ncol=5) +
  theme_bw() +
  scale_color_manual(values = c("#006bff", "#ff0000")) +
  labs(title="Top 10 Key Genes per Cluster",
       y="VST expression",
       x="Time point")

print(top_gene_plot)
#Plot is correctly ordered, but need to have short gene names rather than LOC ID

#Merge top_gene_sumary to add short gene names
top_gene_summary <- top_gene_summary %>%
  left_join(gene_map1, by = c("gene" = "gene"))

#Check merge
if(!"short_name" %in% colnames(top_gene_summary)){
  stop("Merge failed: short_name not found.")
}

# Get short names in same order as gene_order
gene_order_short <- gene_map1$short_name[
  match(gene_order, gene_map1$gene)
]

top_gene_summary$gene_label <- factor(
  top_gene_summary$short_name,
  levels = gene_order_short
)

#PLOT, now showing short gene names rather than LOC ID
top_gene_plot <- ggplot(top_gene_summary,
                        aes(x=time_point,
                            y=mean_expr,
                            color=treatment,
                            group=treatment)) +
  geom_line(size=0.8) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_expr-se,
                    ymax=mean_expr+se),
                width=0.1) +
  facet_wrap(~gene_label, scales="free_y", ncol=5) +
  theme_bw() +
  scale_color_manual(values = c("#006bff", "#ff0000")) +
  labs(title="Top 10 Key Genes per Cluster",
       y="VST expression",
       x="Time point")

print(top_gene_plot)



































##Old lookup table is not complete, need to update it using the GTF file (BCF_000002235.5_Spur_5.0_genomic.gtf, from NCBI) by pulling gene_id, db_xref, and product

# Load necessary library
library(dplyr)

# Initialize vectors to store gene_ids, products, and db_gene_ids
gene_ids <- c()
products <- c()
db_gene_ids <- c()

# Read the GTF file
gtf_file <- "GCF_000002235.5_Spur_5.0_genomic.gtf"
con <- file(gtf_file, "r") #opens the GTF file for reading and assigns the connection to the variable con.

while (TRUE) {
  line <- readLines(con, n = 1)
  if (length(line) == 0) break
  
  if (grepl("gene_id", line) && grepl("product", line) && grepl("db_xref \"GeneID:", line)) {
    # Extract gene_id
    gene_id <- sub(".*gene_id \"(LOC[^\"]+)\".*", "\\1", line)
    
    # Extract product
    product <- sub(".*product \"([^\"]+)\".*", "\\1", line)
    
    # Extract db_xref GeneID
    db_gene_id <- sub(".*db_xref \"GeneID:([^\"]+)\".*", "\\1", line)
    
    # Append to vectors
    gene_ids <- c(gene_ids, gene_id)
    products <- c(products, product)
    db_gene_ids <- c(db_gene_ids, db_gene_id)
  }
}

close(con) #close connection

# Create a dataframe from the vectors
lookup_df <- data.frame(gene_id = gene_ids, product = products, db_gene_id = db_gene_ids)


#Some entries in the gene_id column are wrong, too much info. Take out only the gene_id.
#Function to extract gene ID
extract_gene_id <- function(line) {
  match <- regmatches(line, regexpr('gene_id "([^"]+)"', line))
  if (length(match) > 0) {
    return(sub('gene_id "([^"]+)"', '\\1', match))
  } else {
    return(line)
  }
}

lookup_df$gene_id <- sapply(lookup_df$gene_id, extract_gene_id)

# Remove duplicates
lookup_df_dupsremoved <- lookup_df %>% distinct(gene_id, .keep_all = TRUE)

# Save the final lookup dataframe to a CSV file
#write.csv(lookup_df_dupsremoved, "updated_SP_lookuptable_from_GTF.csv", row.names = FALSE)
#I manually added the few genes it missed in excel.

 
# Load lookup table (.csv)
lookup = read.csv('updated_SP_lookuptable_from_GTF.csv', header = TRUE)
lookup

##Merge DEG2_list_with_clusters with lookup
DEG2_list_with_clusters_genedescriptions <- merge(DEG2_list_with_clusters, lookup, by = "gene_id", all.x = TRUE)

# View the merged dataframe
head(DEG2_list_with_clusters_genedescriptions)
#Save DEG2_list_with_clusters_genedescriptions
write.csv(DEG2_list_with_clusters_genedescriptions, file="DEG2_list_annotated.csv")



#DEG2_list <- read.csv("DEG2_list.csv") #to read back in








