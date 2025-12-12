#HTSEQ Counts Data from S.purpuratus UVB bulkRNAseq experiment: manipulate counts files.
#Updated 1/31/25 by RMK
#HTSEQ was performed perviously using the GMGI RHEL server. Output .counts files are used as input here.
#Each counts file was opened manually, and the last 4 lines were deleted, and the file was resaved.

#This R script is simply to read in each .counts file as a dataframe, name the columns, and then merge everything together. 

# Define a function to generate the new dataframe name based on the file name
generate_df_name <- function(file_name) {
  parts <- unlist(strsplit(file_name, "_"))
  replicate <- switch(parts[1],
                      "SP1" = "A",
                      "SP3" = "B",
                      "SP4" = "C")
  condition <- ifelse(grepl("^C", parts[2]), "Control", "Zap")
  time <- gsub("^C|^Z", "", parts[2])
  paste0(condition, "_", time, "hr_", replicate)
}

# List all counts files in the current directory with the pattern *.combined.counts
counts_files <- list.files(pattern = "*.combined.counts")

# Loop through each file, read it in, rename columns, and assign to a new dataframe
for (file in counts_files) {
  df_name <- generate_df_name(file)
  df <- data.frame(read.table(file, sep = " "))
  colnames(df) <- c("target_id", paste0(df_name, "_count"))
  assign(df_name, df)
}

# Verify the individual dataframes
ls(pattern = "Control|Zap")

#Merge all dataframes into one final summary dataframe
# List of all your dataframe names
df_names <- c("Control_0hr_A", "Control_0hr_B", "Control_0hr_C", "Zap_0hr_A", "Zap_0hr_B", "Zap_0hr_C", 
  "Control_1hr_A", "Control_1hr_B", "Control_1hr_C", "Zap_1hr_A", "Zap_1hr_B", "Zap_1hr_C",
  "Control_3hr_A", "Control_3hr_B", "Control_3hr_C", "Zap_3hr_B", "Zap_3hr_C",
  "Control_6hr_A", "Control_6hr_B", "Control_6hr_C", "Zap_6hr_A", "Zap_6hr_B", "Zap_6hr_C",
  "Control_24hr_A", "Control_24hr_B", "Control_24hr_C", "Zap_24hr_A", "Zap_24hr_B", "Zap_24hr_C")

# Create a list of dataframes
df_list <- lapply(df_names, get)

# Merge all dataframes by "target_id"
BulkRNASeq_CountsTable <- Reduce(function(x, y) merge(x, y, by = "target_id", all = TRUE), df_list)

# Verify the merged dataframe
head(BulkRNASeq_CountsTable)

#Write to final output file
write.table(BulkRNASeq_CountsTable, "BulkRNASeq_SP_timecourse_combined_readCount.tsv", sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

#Remove any genes that have 0 counts across all samples
# Load the data
BulkRNASeq_CountsTable <- read.table("BulkRNASeq_SP_timecourse_combined_readCount.tsv", sep = "\t", header = TRUE)

# Filter out rows where all counts across the 29 samples are 0
filtered_counts <- BulkRNASeq_CountsTable[rowSums(BulkRNASeq_CountsTable[, -1] != 0) > 0, ]

# Write the filtered data to a new file
write.table(filtered_counts, "Filtered_BulkRNASeq_SP_timecourse_combined_readCount.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Verify the filtered data
head(filtered_counts)

#Right now filtered_counts has genes as rows and samples as columns.
#Need to transpose for TPM analysis: need genes as columns and samples as rows.

# Transpose the dataframe
transposed_counts <- t(filtered_counts)

#Convert back to dataframe
test <- as.data.frame(transposed_counts)
write.table(test, "TimecourseRNA_noZeros_Transposed.tsv", sep = "\t", quote =FALSE)
#Also write file without row names for downstream analysis pipeline
write.table(test, "TimecourseRNA_noZeros_Transpose_noRowNames.csv", row.names = FALSE, col.names= FALSE, sep = ",", quote = FALSE)




























         

