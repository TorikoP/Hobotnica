library(dplyr)
library(arrow)
library(doSNOW)
library(GADES)
library(parallel)

input_dir <- "/10tb/vpal/datasets/feather_ds_no_age_39_ds_only_interception"
filtered_sites_dir <- "/home/vpal/hobotnica/Feature_selection/diff_methylation/diff_met_LogFC>0.035_and_pval_0.01"
output_dir <- "/home/vpal/hobotnica/Feature_selection/diff_methylation/sig_diff_met_hobotnica"

files <- list.files(input_dir, full.names = TRUE)
filtered_files <- list.files(filtered_sites_dir, full.names = TRUE)

num_cores <- max(1, detectCores() - 10)
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

calculate_h_score <- function(submatrix, annotation) {
  DistMatrix <- mtrx_distance(submatrix, metric = 'kendall', type = 'gpu', sparse = FALSE, write = TRUE)
  Hobotnica(DistMatrix, annotation)
}

Hobotnica <- function(distMatrix, annotation) {
  if (typeof(annotation) == "list") {
    annotation <- as.vector(unlist(annotation))
  } else {
    annotation <- as.vector(annotation)
  }
  rank.m <- as.matrix(distMatrix)
  rank.m[lower.tri(rank.m)] <- rank(rank.m[lower.tri(rank.m)])
  rank.m[upper.tri(rank.m)] <- rank(rank.m[upper.tri(rank.m)])

  inclass_sum <- 0
  classes <- unique(annotation)
  Ns <- vector()

  for (i in seq_along(classes)) {
    clas <- classes[i]
    class_samples <- which(annotation == clas)
    l_tmp <- length(class_samples)
    Ns[i] <- l_tmp
    tmp_sum_inclass <- sum(rank.m[class_samples, class_samples])
    inclass_sum <- inclass_sum + tmp_sum_inclass
  }

  Ns_sum <- sum(Ns)
  biggest_bossible_rank <- Ns_sum * (Ns_sum - 1) / 2
  number_of_unique_inclass_elements <- sum(Ns * (Ns - 1)) / 2
  maximal_value <- number_of_unique_inclass_elements * (2 * biggest_bossible_rank - number_of_unique_inclass_elements + 1)
  minimal_value <- number_of_unique_inclass_elements * (1 + number_of_unique_inclass_elements)

  normalization_factor <- maximal_value - minimal_value
  return(max(0, 1 - (inclass_sum - minimal_value) / normalization_factor))
}


results <- foreach(file_path = files, .packages = c("arrow", "dplyr", "GADES")) %dopar% {
  file_name <- basename(file_path)
  dataset_name <- tools::file_path_sans_ext(file_name)
  output_file <- file.path(output_dir, paste0(dataset_name))
 
  ds <- read_feather(file_path)
  annotation <- ds$Condition

  filtered_file <- filtered_files[grep(dataset_name, filtered_files)]
  filtered_sites <- read.csv(filtered_file, stringsAsFactors = FALSE)$X
  columns_to_process <- intersect(colnames(ds), filtered_sites)

  if (length(columns_to_process) == 0) {
    cat("No matching columns found for dataset:", dataset_name, "\n")
    return(NULL)
  }

  submatrix <- t(as.matrix(ds[, columns_to_process]))
  h_score <- calculate_h_score(submatrix, annotation)

  cat("Processed and saved:", file_name, "\n")
}

stopCluster(cl)

all_results <- do.call(rbind, results)
write.table(all_results, file = file.path(output_dir, "all_h_scores.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
