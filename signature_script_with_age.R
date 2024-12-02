#!/usr/bin/env Rscript
library(dplyr)
source("/home/vpal/hobotnica/custom_function.R")
path_sig_matrix <- "/home/vpal/hobotnica/PhenoAgeV2_res_imputed_rand_signature/sig_matrix"
pheno_age_ds <- "/home/vpal/hobotnica/datasets_PhenoAgeV2_imputed_NA"
output_dir <- "/home/vpal/hobotnica/PhenoAgeV2_res_imputed_rand_signature/distrib_h_score_initial"

folder_list <- list.dirs(path_sig_matrix, full.names = TRUE, recursive = FALSE)
file_list_PhenoAge <- list.files(path = pheno_age_ds, full.names = TRUE)

for (i in 1:length(file_list_PhenoAge)) {
    file1 <- file_list_PhenoAge[i]
    sample_id1 <- basename(file1)
    data <- read.csv(file1, sep = ",", header = TRUE, row.names = 1)
    annotation <- data$Condition
    matching_folder <- folder_list[basename(folder_list) == tools::file_path_sans_ext(sample_id1)]
    if (length(matching_folder) > 0) {
        folder <- matching_folder[1]
        H_score_distrib <- ParallelHobotnica(folder, annotation)
        txt_output_path <- file.path(output_dir, paste0("H_score_distrib_", tools::file_path_sans_ext(sample_id1), ".txt"))
        write.table(H_score_distrib, txt_output_path, row.names = FALSE, col.names = FALSE)
    }
}