###############################################################################
## INSTALL AND LOAD PACKAGES
###############################################################################

packages <- c("tidyverse","impute","readr", "imputeLCMD")

if(sum(as.numeric(!packages %in% installed.packages())) != 0){
  instalador <-packages[!packages %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(packages, require, character = T) 
} else {
  sapply(packages, require, character = T) 
}


###############################################################################
## DATA LOADING AND PREPROCESSING
###############################################################################

# Load data
dataset <- read_csv("lfq.proteins.csv")

data <- dataset %>%
  # 1) Keep only Top proteins
  filter(Top == TRUE) %>%
  
  # 2) Select Accession and abundance columns
  select(
    Accession,
    `WT-30_1 Area`, `WT-30_2 Area`, `WT-30_3 Area`,
    `WT-22_1 Area`, `WT-22_2 Area`, `WT-22_3 Area`,
    `A1-30_1 Area`, `A1-30_2 Area`, `A1-30_3 Area`,
    `A1-22_1 Area`, `A1-22_2 Area`, `A1-22_3 Area`
  ) %>%
  
  # 3) Remove " Area" from column names
  rename_with(~ str_replace(.x, " Area$", ""))


###############################################################################
## LOG2 TRANSFORMATION AND ZERO HANDLING
###############################################################################

# Columns of interest
columns_to_transform <- c(
  "A1-30_1", "A1-30_2", "A1-30_3",
  "A1-22_1", "A1-22_2", "A1-22_3",
  "WT-30_1", "WT-30_2", "WT-30_3",
  "WT-22_1", "WT-22_2", "WT-22_3"
)

# Replace zeros with NA and apply log2 transformation
df <- data %>%
  mutate(across(all_of(columns_to_transform), ~ ifelse(. == 0, NA, .))) %>%
  mutate(across(all_of(columns_to_transform), ~ log2(.)))


###############################################################################
## KNN IMPUTATION – MATRIX PREPARATION
###############################################################################

# Extract numeric matrix
mat_numeric <- df[, columns_to_transform]
mat_numeric <- as.matrix(mat_numeric)
mode(mat_numeric) <- "numeric"
rownames(mat_numeric) <- df$Accession

# ---- FIX PROTEIN ORDER (CRITICAL FOR MINPROB & KNN REPRODUCIBILITY) ----
mat_numeric <- mat_numeric[order(rownames(mat_numeric)), ]


###############################################################################
## EXPERIMENTAL GROUP DEFINITION
###############################################################################

groups <- list(
  A1_30 = c("A1-30_1","A1-30_2","A1-30_3"),
  A1_22 = c("A1-22_1","A1-22_2","A1-22_3"),
  WT_30 = c("WT-30_1","WT-30_2","WT-30_3"),
  WT_22 = c("WT-22_1","WT-22_2","WT-22_3")
)


###############################################################################
## KNN OPTIMIZATION FUNCTIONS
###############################################################################

# Function to simulate missing values
simulate_missing <- function(mat, prop = 0.1, seed = 1234){
  set.seed(seed)
  mat_sim <- mat
  idx <- which(!is.na(mat), arr.ind = TRUE)
  n_missing <- round(nrow(idx) * prop)
  missing_idx <- idx[sample(1:nrow(idx), n_missing), ]
  
  for(i in 1:nrow(missing_idx)){
    mat_sim[missing_idx[i, "row"], missing_idx[i, "col"]] <- NA
  }
  
  return(list(mat_sim = mat_sim, missing_idx = missing_idx))
}

# Function to calculate RMSE
calc_rmse <- function(original, imputed, missing_idx){
  diffs <- numeric(nrow(missing_idx))
  for(i in 1:nrow(missing_idx)){
    r <- missing_idx[i, "row"]
    c <- missing_idx[i, "col"]
    diffs[i] <- original[r,c] - imputed[r,c]
  }
  sqrt(mean(diffs^2))
}

# Function to test multiple k values with repeated simulations
test_knn_k <- function(mat_original, k_values = c(3, 5, 7,10, 15, 20), 
                       prop_missing = 0.1, n_iter = 100){
  rmse_results <- matrix(NA, nrow = n_iter, ncol = length(k_values))
  colnames(rmse_results) <- k_values
  
  for(j in 1:n_iter){
    sim <- simulate_missing(mat_original, prop = prop_missing, seed = 1234+j)
    mat_sim <- sim$mat_sim
    missing_idx <- sim$missing_idx
    
    for(i in seq_along(k_values)){
      k <- k_values[i]
      imputed <- impute.knn(mat_sim, k = k)$data
      rmse_results[j,i] <- calc_rmse(mat_original, imputed, missing_idx)
    }
  }
  colMeans(rmse_results, na.rm = TRUE)
}


###############################################################################
## OPTIMAL K SELECTION PER GROUP
###############################################################################

mat_numeric_original <- mat_numeric
optimal_k <- list()

for(grp_name in names(groups)){
  cat("Testing KNN for group:", grp_name, "\n")
  
  group_mat <- mat_numeric_original[, groups[[grp_name]], drop = FALSE]
  group_mat <- group_mat[rowSums(!is.na(group_mat)) > 0, ]
  
  rmse_res <- test_knn_k(group_mat, k_values = c(3,5,7,10,15,20), n_iter = 50)
  
  optimal_k[[grp_name]] <- names(rmse_res)[which.min(rmse_res)]
  
  print(rmse_res)
  cat("Optimal k for", grp_name, "=", optimal_k[[grp_name]], "\n\n")
  
  plot(as.numeric(names(rmse_res)), rmse_res, type="b",
       xlab="k", ylab="RMSE", main=paste("Group", grp_name))
}

optimal_k


###############################################################################
## KNN IMPUTATION (PROTEINS WITH 2 OF 3 REPLICATES)
###############################################################################

imputed_groups <- list()

for (grp_name in names(groups)) {
  
  group_mat <- mat_numeric[, groups[[grp_name]], drop = FALSE]
  n_obs <- rowSums(!is.na(group_mat))
  
  idx_2of3 <- which(n_obs == 2)
  if (length(idx_2of3) == 0) next
  
  mat_knn <- group_mat[n_obs >= 2, , drop = FALSE]
  
  knn_res <- impute.knn(
    as.matrix(mat_knn),
    k = 15
  )$data
  
  prot_2of3 <- rownames(group_mat)[idx_2of3]
  
  for (prot in prot_2of3) {
    na_pos <- is.na(group_mat[prot, ])
    group_mat[prot, na_pos] <- knn_res[prot, na_pos]
  }
  
  imputed_groups[[grp_name]] <- group_mat[prot_2of3, , drop = FALSE]
}

# Merge KNN-imputed values into full matrix
df_final <- mat_numeric
for (grp_name in names(imputed_groups)) {
  df_final[
    rownames(imputed_groups[[grp_name]]),
    groups[[grp_name]]
  ] <- imputed_groups[[grp_name]]
}


###############################################################################
## MINPROB IMPUTATION
###############################################################################

set.seed(1234)

# Function to apply MinProb imputation within each group
impute_minprob_group <- function(df_final, group_cols, q = 0.01) {
  
  df_tmp <- df_final[, group_cols, drop = FALSE]
  group_absent <- apply(df_tmp, 1, function(x) all(is.na(x)))
  
  df_imputed <- df_tmp
  df_imputed[!group_absent, ] <- impute.MinProb(
    as.matrix(df_tmp[!group_absent, ]),
    q = q
  )
  
  df_imputed[group_absent, ] <- NA
  return(df_imputed)
}

data_minprob <- df_final

data_minprob[, groups$A1_30] <- impute_minprob_group(df_final, groups$A1_30)
data_minprob[, groups$A1_22] <- impute_minprob_group(df_final, groups$A1_22)
data_minprob[, groups$WT_30] <- impute_minprob_group(df_final, groups$WT_30)
data_minprob[, groups$WT_22] <- impute_minprob_group(df_final, groups$WT_22)

write.csv(data_minprob, "data_refined_KNN+MinProb.csv", row.names = TRUE)


###############################################################################
## STUDENT'S T-TEST
###############################################################################

calc_ttest_student_p <- function(mat, group1_cols, group2_cols) {
  apply(mat[, c(group1_cols, group2_cols)], 1, function(x) {
    
    g1 <- x[group1_cols]
    g2 <- x[group2_cols]
    
    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]
    
    if (length(g1) < 2 || length(g2) < 2) return(NA)
    
    tryCatch(
      t.test(g1, g2, var.equal = TRUE)$p.value,
      error = function(e) NA
    )
  })
}


###############################################################################
## LOG2 FOLD CHANGE
###############################################################################

calc_log2fc <- function(mat, group1_cols, group2_cols) {
  apply(mat[, c(group1_cols, group2_cols)], 1, function(x) {
    
    g1 <- x[group1_cols]
    g2 <- x[group2_cols]
    
    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]
    
    if (length(g1) == 0 || length(g2) == 0) return(NA)
    
    mean(g1) - mean(g2)
  })
}


###############################################################################
## APPLY STATISTICAL ANALYSES
###############################################################################

data_minprob_df <- as.data.frame(data_minprob)

data_minprob_df$pvalue_A1_22_vs_A1_30 <- calc_ttest_student_p(
  data_minprob_df, groups$A1_22, groups$A1_30
)

data_minprob_df$pvalue_WT_22_vs_WT_30 <- calc_ttest_student_p(
  data_minprob_df, groups$WT_22, groups$WT_30
)

data_minprob_df$log2FC_A1_22_vs_A1_30 <- calc_log2fc(
  data_minprob_df, groups$A1_22, groups$A1_30
)

data_minprob_df$log2FC_WT_22_vs_WT_30 <- calc_log2fc(
  data_minprob_df, groups$WT_22, groups$WT_30
)


###############################################################################
## SIGNIFICANCE FILTERING
###############################################################################

data_minprob_df$Significance_A1_22_vs_A1_30 <-
  with(data_minprob_df,
       !is.na(pvalue_A1_22_vs_A1_30) &
         pvalue_A1_22_vs_A1_30 < 0.05 &
         abs(log2FC_A1_22_vs_A1_30) >= 0.26
  )

data_minprob_df$Significance_WT_22_vs_WT_30 <-
  with(data_minprob_df,
       !is.na(pvalue_WT_22_vs_WT_30) &
         pvalue_WT_22_vs_WT_30 < 0.05 &
         abs(log2FC_WT_22_vs_WT_30) >= 0.26
  )


###############################################################################
## DIRECTION (UP / DOWN / NS)
###############################################################################

data_minprob_df$Direction_A1_22 <- with(
  data_minprob_df,
  ifelse(
    Significance_A1_22_vs_A1_30 & log2FC_A1_22_vs_A1_30 >= 0.26, "Up",
    ifelse(
      Significance_A1_22_vs_A1_30 & log2FC_A1_22_vs_A1_30 <= -0.26, "Down",
      "NS"
    )
  )
)

data_minprob_df$Direction_WT_22 <- with(
  data_minprob_df,
  ifelse(
    Significance_WT_22_vs_WT_30 & log2FC_WT_22_vs_WT_30 >= 0.26, "Up",
    ifelse(
      Significance_WT_22_vs_WT_30 & log2FC_WT_22_vs_WT_30 <= -0.26, "Down",
      "NS"
    )
  )
)


###############################################################################
## SUMMARY COUNTS
###############################################################################

n_sig_A1 <- sum(data_minprob_df$Significance_A1_22_vs_A1_30, na.rm = TRUE)
n_sig_WT <- sum(data_minprob_df$Significance_WT_22_vs_WT_30, na.rm = TRUE)

cat("A1_22 vs A1_30 - significant proteins:", n_sig_A1, "\n")
cat("WT_22 vs WT_30 - significant proteins:", n_sig_WT, "\n")

cat("A1_22 vs A1_30\n")
print(table(data_minprob_df$Direction_A1_22))

cat("\nWT_22 vs WT_30\n")
print(table(data_minprob_df$Direction_WT_22))


###############################################################################
## EXPORT FINAL DATASET
###############################################################################
data_minprob_df <- data.frame(
  Accession = rownames(data_minprob_df),
  data_minprob_df,
  row.names = NULL,
  check.names = FALSE
)

data_minprob_df <- data_minprob_df %>%
  left_join(
    dataset %>%
      select(
        Accession,
        `Protein Group`,
        `Avg. Mass`,
        Description,
        `Coverage (%)`,
        `#Peptides`,
        `#Unique`,
        PTM
      ),
    by = "Accession"
  )


write.csv(
  data_minprob_df,
  "LFQ_ouput.csv",
  row.names = FALSE
)

