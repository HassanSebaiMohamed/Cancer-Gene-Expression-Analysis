# Cancer Gene Expression Analysis
wt_seq <- "ATGCGTACGGT"
mut_seq <- "ATGCCTACGGT"

set.seed(123)
counts <- matrix(sample(0:20, 48, replace = TRUE), nrow = 12)
rownames(counts) <- paste0("Gene", 1:12)
colnames(counts) <- paste0("Sample", 1:4)

# Add missing values
counts[3,1] <- NA
counts[8,4] <- NA

cat("Original Data:\n")
print(counts)


#Global Alignment (Simple Similarity) 
wt <- strsplit(wt_seq, "")[[1]]
mut <- strsplit(mut_seq, "")[[1]]

matches <- sum(wt == mut)
similarity <- matches / length(wt)

cat("\nSequence Matches:", matches, "\n")
cat("Similarity Score:", similarity, "\n")


# Filter Genes 
gene_sums <- rowSums(counts, na.rm = TRUE)
filtered_counts <- counts[gene_sums > 10, ]

cat("\nFiltered Data:\n")
print(filtered_counts)


# Handle Missing Values 
na_count <- sum(is.na(filtered_counts))
cat("\nNumber of Missing Values:", na_count, "\n")

# Replace NA with row mean
for(i in 1:nrow(filtered_counts)) {
  row_mean <- mean(filtered_counts[i, ], na.rm = TRUE)
  filtered_counts[i, is.na(filtered_counts[i, ])] <- row_mean
}

cat("\nAfter Replacing NA:\n")
print(filtered_counts)


#  CPM Normalization 
library(edgeR)

cpm_data <- cpm(filtered_counts)

cat("\nCPM Normalized Data:\n")
print(cpm_data)


#Log Transformation 
log_data <- log2(cpm_data + 1)

cat("\nLog Transformed Data:\n")
print(log_data)


# Z-score Scaling 
z_data <- t(scale(t(log_data)))

cat("\nZ-score Scaled Data:\n")
print(z_data)


# Visualization
par(mfrow=c(1,2))

boxplot(filtered_counts,
        main = "Raw Data",
        col = "red")

boxplot(z_data,
        main = "Processed Data",
        col = "green")