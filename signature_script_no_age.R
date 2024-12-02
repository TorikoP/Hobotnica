#!/usr/bin/env Rscript
library(dplyr)
library(GADES)
source("/home/vpal/hobotnica/custom_function.R")
path_sig_matrix <- "/home/vpal/hobotnica/PhenoAgeV2_res_imputed_rand_signature/sig_matrix_no_age"
substract_age <- "/home/vpal/hobotnica/substracted_age_PhenoAgeV2_imputed_NA"
output_dir <- "/home/vpal/hobotnica/PhenoAgeV2_res_imputed_rand_signature/distrib_h_score_no_age"
existing_results_path <- "/home/vpal/hobotnica/PhenoAgeV2_res_imputed_rand_signature/PhenoAgeV2_H_scores_no_age.csv"

folder_list <- list.dirs(path_sig_matrix, full.names = TRUE, recursive = FALSE)
file_list_PhenoAge <- list.files(path = substract_age, full.names = TRUE)

result_table <- data.frame(
  Dataset_ID = character(), 
  H_score_no_age = numeric(), 
  p_value_no_age = numeric(), 
  stringsAsFactors = FALSE
)
for (i in 1:length(file_list_PhenoAge)) {
    file1 <- file_list_PhenoAge[i]
    sample_id1 <- basename(file1)
    sample_id_no_ext <- tools::file_path_sans_ext(sample_id1) 
    txt_output_path <- file.path(output_dir, paste0("H_score_distrib_", sample_id_no_ext, ".txt"))
    data <- read.csv(file1, sep = ",", header = TRUE, row.names = 1)
    annotation <- data$Condition
    if ("Condition" %in% colnames(data) || "Age" %in% colnames(data)) {
      matrix <- data %>% select(-Condition, -Age)
    } else {
      matrix <- data
    }

    as_matrix <- t(as.matrix(matrix))
    distMatrix <- mtrx_distance(as_matrix, metric = 'kendall', type='gpu', sparse=F, write=T)
    h_score <- Hobotnica(distMatrix, annotation)
    
    matching_folder <- folder_list[basename(folder_list) == sample_id_no_ext]#check to start calculation on the matrix from the same dataset
    if (length(matching_folder) > 0) {
        folder <- matching_folder[1]
        H_score_distrib <- ParallelHobotnica(folder, annotation)
        pval <- (1 + sum(H_score_distrib >= h_score)) / length(H_score_distrib)
        write.table(H_score_distrib, txt_output_path, row.names = FALSE, col.names = FALSE)
        new_row <- data.frame(Dataset_ID = sample_id_no_ext, H_score_no_age = h_score, p_value_no_age = pval)
        write.table(new_row, existing_results_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    } else {
        warning(paste("No matching folder found for", sample_id_no_ext))
    }
}