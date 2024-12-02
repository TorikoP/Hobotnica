library(dplyr)
library(arrow)
library(doSNOW)
library(GADES)
library(parallel)

input_dir <- "/home/vpal/hobotnica/All_datasets/datasets/data_feather_no_age_big_ds"
sites_folder <- "/home/vpal/hobotnica/All_datasets/single_sites_analysis/a_sites_for_clocks"
output_base_dir <- "/home/vpal/hobotnica/All_datasets/single_sites_analysis/clocks_sites_h_score"

num_cores <- max(1, detectCores() - 10)
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

columns_Hobotnica_single <- function(ds, annotation, columns_to_process) {
  pb <- txtProgressBar(min = 0, max = length(columns_to_process), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  H_scores <- foreach(col_idx = seq_along(columns_to_process), .packages = c("GADES"),
                      .combine = rbind, .options.snow = opts, 
                      .export = c("Hobotnica")) %dopar% {

    colname <- columns_to_process[col_idx]
    column_data <- ds[[colname]]
    as_matrix <- t(matrix(column_data, ncol = 1))
    DistMatrix <- mtrx_distance(as_matrix, metric = 'euclidean', type = 'gpu', sparse = FALSE, write = TRUE)
    H_result <- Hobotnica(DistMatrix, annotation)

    data.frame(Column = colname, Result = H_result, stringsAsFactors = FALSE)
  }
  return(H_scores)
}

Hobotnica <- function(distMatrix, annotation){
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
    biggest_possible_rank <- Ns_sum * (Ns_sum - 1) / 2
    number_of_unique_inclass_elements <- sum(Ns * (Ns - 1)) / 2
    maximal_value <- number_of_unique_inclass_elements * (2 * biggest_possible_rank - number_of_unique_inclass_elements + 1)
    minimal_value <- number_of_unique_inclass_elements * (1 + number_of_unique_inclass_elements)

    normalization_factor <- maximal_value - minimal_value
    return(max(0, 1 - (inclass_sum - minimal_value) / normalization_factor))
}

site_files <- list.files(sites_folder, full.names = TRUE)
for (site_file in site_files) {
  site_list <- read.table(site_file, header = FALSE, stringsAsFactors = FALSE)[[1]]
  site_file_name <- tools::file_path_sans_ext(basename(site_file))

  clock_output_dir <- file.path(output_base_dir, site_file_name)
  dir.create(clock_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (file_path in list.files(input_dir, full.names = TRUE)) {
    dataset_name <- tools::file_path_sans_ext(basename(file_path))
    output_file <- file.path(clock_output_dir, paste0(dataset_name))

    if (file.exists(output_file)) {
      cat("Skipping already processed dataset for clock file:", site_file_name, "Dataset:", dataset_name, "\n")
      next
    }

    ds <- read_feather(file_path)
    annotation <- ds$Condition
    columns_to_process <- intersect(site_list, colnames(ds))
    if (length(columns_to_process) == 0) {
      cat("No matching columns found for sites in:", site_file_name, "Dataset:", dataset_name, "\n")
      next
    }

    H_scores <- columns_Hobotnica_single(ds, annotation, columns_to_process)
    
    write.table(H_scores, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    cat("Processed and saved:", dataset_name, "for clock file:", site_file_name, "\n")
  }
}

stopCluster(cl)