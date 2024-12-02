
library(dplyr)
library(arrow)
library(doSNOW)
library(GADES)
library(parallel)

input_dir <- "/home/vpal/hobotnica/All_datasets/data_feather"
output_dir <- "/home/vpal/hobotnica/All_datasets/H_scores_triplets"
batch_size <- 1000

files <- list.files(input_dir, full.names = TRUE)

num_cores <- max(1, detectCores() - 10)
cl <- makeCluster(num_cores)
registerDoSNOW(cl)

columns_Hobotnica_triplets <- function(ds, annotation) {

  excluded_cols <- c("Condition", "Age", "index")
  columns_to_process <- setdiff(colnames(ds), excluded_cols)
  total_columns <- length(columns_to_process)

  if (total_columns %% 3 != 0) {
    total_columns <- total_columns - (total_columns %% 3)
  }

  H_scores <- foreach(triplet_idx = seq(1, total_columns, by = 3), 
                      .packages = c("GADES"), 
                      .combine = rbind, 
                      .export = c("Hobotnica")) %dopar% {
    
    colname1 <- columns_to_process[triplet_idx]
    colname2 <- columns_to_process[triplet_idx + 1]
    colname3 <- columns_to_process[triplet_idx + 2]

    column_data1 <- ds[[colname1]]
    column_data2 <- ds[[colname2]]
    column_data3 <- ds[[colname3]]

    as_matrix <- t(cbind(column_data1, column_data2, column_data3))
    DistMatrix <- mtrx_distance(as_matrix, metric = 'kendall', type = 'gpu', sparse = FALSE, write = TRUE)

    H_result <- Hobotnica(DistMatrix, annotation)

    data.frame(Column1 = colname1, Column2 = colname2, Column3 = colname3, Result = H_result, stringsAsFactors = FALSE)
  }

  return(H_scores)
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


for (file_path in files) {
  file_name <- basename(file_path)
  output_file <- file.path(output_dir, file_name)

  if (file.exists(output_file)) {
    cat("Skipping already processed file:", file_name, "\n")
    next
  }

  ds <- read_feather(file_path)
  annotation <- ds$Condition

  H_scores <- columns_Hobotnica_triplets(ds, annotation)

  write.table(H_scores, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  cat("Processed and saved:", file_name, "\n")
}

stopCluster(cl)