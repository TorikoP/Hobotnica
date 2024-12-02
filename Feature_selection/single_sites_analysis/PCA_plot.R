library(dplyr)
library(ggplot2)
library(tidyr)

# Input folder containing all the files
input_folder <- "/home/vpal/hobotnica/All_datasets/single_sites_analysis/h_score_result"

# Output PDF file path
output_pdf <- "/home/vpal/hobotnica/All_datasets/single_sites_analysis/pca_plot.pdf"

# Get the list of all files in the folder
files <- list.files(input_folder, full.names = TRUE)

# Combine all files into a single data frame
combined_data <- lapply(files, function(file) {
  data <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  colnames(data) <- c("Column", "Result")  # Ensure consistent naming
  data$File <- basename(file)  # Add the file name to identify the source
  return(data)
}) %>%
  bind_rows()  # Combine all files into one data frame

# Reshape the data into wide format: Rows = Columns (sites), Columns = Files
pca_data <- combined_data %>%
  select(Column, File, Result) %>%
  pivot_wider(names_from = File, values_from = Result)  # Convert to wide format

# Save 'Column' as a separate column instead of rownames
site_names <- pca_data$Column  # Extract site names
pca_data <- pca_data %>%
  select(-Column)  # Remove the 'Column' column

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# Prepare PCA results for plotting
pca_plot_data <- as.data.frame(pca_result$x)  # PCA-transformed data
pca_plot_data$Site <- site_names  # Add site names explicitly as a column

# Create the PCA plot
pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, label = Site)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.5, size = 2.5) +
  labs(
    title = "PCA of H-Scores for Sites Across All Files",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

# Save the PCA plot to a PDF
pdf(output_pdf, width = 10, height = 7)  # Open PDF device
print(pca_plot)  # Print the plot to the PDF
dev.off()  # Close the PDF device

cat("PCA plot saved to:", output_pdf, "\n")