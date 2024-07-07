# Load necessary packages if not already installed
# install.packages("readr")
library(readr)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("Rtsne")
library(Rtsne)

# File paths for the four CSV files
file_paths <- c("GSM958226-tbl-1.txt", "GSM958227-tbl-1.txt", "GSM958220-tbl-1.txt", "GSM958221-tbl-1.txt")

# Read the four files
file1 <- read.table("GSM958226-tbl-1.txt", header = FALSE)  # Replace "file1.txt" with your actual file path
file2 <- read.table("GSM958227-tbl-1.txt", header = FALSE)  # Replace "file2.txt" with your actual file path

file3 <- read.table("GSM958220-tbl-1.txt", header = FALSE)  # Replace "file3.txt" with your actual file path
file4 <- read.table("GSM958221-tbl-1.txt", header = FALSE)  # Replace "file4.txt" with your actual file path

# Rename columns
colnames(file1) <- c("ID", "Value1")
colnames(file2) <- c("ID", "Value2")
colnames(file3) <- c("ID", "Value3")
colnames(file4) <- c("ID", "Value4")

# Merge files based on ID
merged_files <- merge(file1, file2, by = "ID", all = TRUE)
merged_files <- merge(merged_files, file3, by = "ID", all = TRUE)
merged_files <- merge(merged_files, file4, by = "ID", all = TRUE)

# Fill missing values with NA
merged_files[is.na(merged_files)] <- 0  # Change 0 to NA if you prefer NA for missing values

# View the merged table
print(merged_files)

data_for_clustering <- merged_files[, c("Value1", "Value2", "Value3", "Value4")]

# Perform dimensionality reduction using t-SNE (t-distributed Stochastic Neighbor Embedding)
tsne_result <- Rtsne(data_for_clustering, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

# Extract the 2D representation of the data
tsne_data <- as.data.frame(tsne_result$Y)
colnames(tsne_data) <- c("Dim1", "Dim2")

# Perform clustering on the data_for_clustering
num_clusters <- 5  # Change this number as needed
clusters <- kmeans(data_for_clustering, centers = num_clusters)

# Add cluster assignments to tsne_data
tsne_data$Cluster <- as.factor(clusters$cluster)

# Plotting the clustered data using ggplot2
ggplot(tsne_data, aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point() +
  labs(title = "t-SNE Visualization of Clustering") +
  theme_minimal()

