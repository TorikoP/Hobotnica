library(dplyr)
library(arrow)
library(doParallel)
library(dplyr)
library(amap)

input_dir <- "/10tb/vpal/datasets/feather_ds_no_age_39_ds_only_interception"
filtered_sites_dir <- "/home/vpal/hobotnica/Feature_selection/diff_methylation/diff_met_LogFC>0.035_and_pval_0.01"
distribution_path <- "/10tb/vpal/diff_methylation/distrib_h_score"
output_results_file <- "/home/vpal/hobotnica/Feature_selection/diff_methylation/diff_met_h_scores_with_pval.csv"

file_list1 <- list.files(path = input_dir, full.names = TRUE)
file_list2 <- list.files(path = filtered_sites_dir, full.names = TRUE)

if (length(file_list1) != length(file_list2)) {
  stop("The number of files in each directory does not match.")
}

new_result_table <- data.frame(
  Dataset_ID = character(),
  H_score = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

kendall_dist <- function(ds_samples) {
    ds_samples <- as.data.frame(ds_samples)
    if ("Condition" %in% colnames(ds_samples) || "Age" %in% colnames(ds_samples)) {
        matrix <- dplyr::select(ds_samples, -Condition, -Age)
    } else {
        matrix <- ds_samples
    }
    distMatrix <- Dist(matrix, method = "kendall", nbproc = 2) 
    return(distMatrix)
}

distrib_and_hobot <- function(distMatrix, annotation, number) {
    h_val <- Hobotnica(distMatrix, annotation)
    distMatrix_2 = as.matrix(distMatrix)
    distribution <- RandomeDistribution(distMatrix_2, annotation, number)
    pval <- (1+ sum(distribution >= h_val))/number
    result <- list(pval = pval, h_val = h_val, random_h_scores = distribution)
    return(result)
}

RandomeDistribution <- function(distMatrix, annotation, nPermutations) {
    if (length(dim(distMatrix)) != 2) {
        stop("The distMatrix dim length should be equal 2, stopping.")
    }
    if (dim(distMatrix)[1] != dim(distMatrix)[2]) {
        stop("distMatrix should be a square matrix, stopping")
    }
    
    H_scores <- foreach (i = 1:nPermutations, .packages = c("dplyr"),
                      .export = c("Hobotnica"),
                      .combine = 'c') %dopar% {
        permutedAnnotation <- sample(annotation, length(annotation), replace=FALSE)
        H_result <- Hobotnica(distMatrix, permutedAnnotation)
        
        return(H_result)
    }

}

Hobotnica <- function(distMatrix, annotation){
    if (typeof(annotation) == "list") {
        annotation <- as.vector(unlist(annotation))
    } else {
        annotation <- as.vector(annotation)
    }
    rank.m <- as.matrix(distMatrix) # transform distance matrix to matrix object
    rank.m[lower.tri(rank.m)] <- rank(rank.m[lower.tri(rank.m)]) # transform distances to ranks
    rank.m[upper.tri(rank.m)] <- rank(rank.m[upper.tri(rank.m)]) #

    inclass_sum <- 0
    classes <- unique(annotation) # unique classes
    Ns <- vector()

    for (i  in 1:length(classes)){

        clas <- classes[i]
        class_samples <- which(annotation == clas)
        l_tmp <- length(class_samples)
        Ns[i] <- l_tmp
        tmp_sum_inclass <- sum(rank.m[class_samples,class_samples]) # sum of ranks, describing in-class distances
        inclass_sum <- inclass_sum + tmp_sum_inclass


    }
    Ns_sum <- sum(Ns)
    biggest_bossible_rank <-  Ns_sum * (Ns_sum - 1)/2
    number_of_unique_inclass_elements <-  sum(Ns * (Ns-1))/2
    maximal_value <- number_of_unique_inclass_elements * (2*biggest_bossible_rank - number_of_unique_inclass_elements + 1)
    minimal_value <- number_of_unique_inclass_elements* (1 + number_of_unique_inclass_elements)

    normalization_factor <- maximal_value - minimal_value
    return (max(0, 1 - (inclass_sum - minimal_value)/normalization_factor ))

}

for (i in seq_along(file_list1)) {
  file1 <- file_list1[i]
  file2 <- file_list2[i]

  dataset_name <- tools::file_path_sans_ext(basename(file1))
  ds <- read_feather(file1)
  annotation <- ds$Condition

  filtered_sites <- read.csv(file2, stringsAsFactors = FALSE)$X
  columns_to_process <- intersect(colnames(ds), filtered_sites)

  submatrix <- ds[, columns_to_process, drop = FALSE]
  submatrix <- as.data.frame(submatrix)
  rownames(submatrix) <- submatrix$index
  submatrix <- submatrix[, !colnames(submatrix) %in% "index"]
  distMatrix <- kendall_dist(submatrix)
  
  distrib_result <- distrib_and_hobot(distMatrix, annotation, 1000)

  distrib_file <- file.path(distribution_path, paste0(dataset_name, "_distrib.txt"))
  write.table(distrib_result$random_h_scores, file = distrib_file, row.names = FALSE, col.names = FALSE)

  new_result_table <- new_result_table %>%
    add_row(Dataset_ID = dataset_name, H_score = distrib_result$h_val, p_value = distrib_result$pval)

  cat("Processed dataset:", dataset_name, "\n")
}

write.csv(new_result_table, output_results_file, row.names = FALSE)