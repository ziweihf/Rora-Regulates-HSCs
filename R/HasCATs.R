# Title: HasCATs Framework for TF Contribution Analysis
# Description: This script calculates transcription factor (TF) contributions to aging
# using single-cell RNA-seq data, TF activity, and regulatory relationships.

# Load required libraries
library(Seurat)
library(dplyr)
library(Matrix)

#' Load and preprocess single-cell data
#' @param rds_path Path to Seurat RDS file
#' @return List containing expression and TF activity matrices
load_sc_data <- function(rds_path) {
  if (!file.exists(rds_path)) stop("RDS file not found: ", rds_path)
  sc_data <- readRDS(rds_path)
  expr_matrix <- as.matrix(GetAssayData(sc_data, slot = "data", assay = "RNA"))
  tf_act_matrix <- as.matrix(GetAssayData(sc_data, slot = "data", assay = "AUC"))
  
  # Ensure common TFs between expression and activity matrices
  common_tfs <- intersect(rownames(expr_matrix), rownames(tf_act_matrix))
  if (length(common_tfs) == 0) stop("No common TFs found between matrices")
  
  list(
    expr = expr_matrix[common_tfs, , drop = FALSE],
    tf_act = tf_act_matrix[common_tfs, , drop = FALSE]
  )
}

#' Calculate TF expression-activity product
#' @param expr_matrix Gene expression matrix
#' @param tf_act_matrix TF activity matrix
#' @return Matrix of TF expression × activity
calc_tf_expr_activity <- function(expr_matrix, tf_act_matrix) {
  tf_expr_activity <- expr_matrix * tf_act_matrix
  return(tf_expr_activity)
}

#' Load and preprocess TF regulatory data
#' @param reg_path Path to TF regulation file
#' @return Data frame with adjusted mode of regulation (mor)
load_tf_regulation <- function(reg_path) {
  if (!file.exists(reg_path)) stop("Regulation file not found: ", reg_path)
  tf_reg <- read.table(reg_path, header = TRUE, sep = "\t")
  
  # Adjust mor based on senescence effect
  tf_reg$mor <- ifelse(tf_reg$Senescence.Effect == "Induces",
                       tf_reg$mor,
                       -1 * tf_reg$mor)
  return(tf_reg)
}

#' Calculate edge weights with sigmoid transformation
#' @param tf_reg Data frame with TF-target relationships
#' @param expr_matrix Gene expression matrix
#' @param tf_list List of TFs to process
#' @return Matrix of sigmoid-transformed edge weights
calc_edge_weights <- function(tf_reg, expr_matrix, tf_list) {
  # Filter valid TF-target pairs
  tf_targets <- tf_reg %>%
    filter(source %in% tf_list, target %in% rownames(expr_matrix)) %>%
    select(source, target, mor) %>%
    distinct()
  
  if (nrow(tf_targets) == 0) stop("No valid TF-target pairs found")
  
  # Initialize result matrix
  result <- matrix(0,
                   nrow = length(unique(tf_targets$source)),
                   ncol = ncol(expr_matrix),
                   dimnames = list(unique(tf_targets$source), colnames(expr_matrix)))
  
  # Vectorized edge weight calculation
  for (tf in rownames(result)) {
    tf_data <- tf_targets %>% filter(source == tf)
    if (nrow(tf_data) == 0) next
    expr_subset <- expr_matrix[tf_data$target, , drop = FALSE]
    weights <- t(tf_data$mor) %*% expr_subset
    result[tf, ] <- weights
  }
  
  # Apply sigmoid transformation
  sigmoid <- function(x) 1 / (1 + exp(-x))
  result <- sigmoid(result)
  
  return(result)
}

#' Calculate final TF contribution scores
#' @param tf_expr_activity TF expression-activity matrix
#' @param edge_weights Edge weight matrix
#' @return Matrix of final TF contribution scores
calc_tf_contribution <- function(tf_expr_activity, edge_weights) {
  common_tfs <- intersect(rownames(tf_expr_activity), rownames(edge_weights))
  if (length(common_tfs) == 0) stop("No common TFs for final score calculation")
  
  tf_expr_activity <- tf_expr_activity[common_tfs, , drop = FALSE]
  edge_weights <- edge_weights[common_tfs, , drop = FALSE]
  
  final_scores <- tf_expr_activity * edge_weights
  return(final_scores)
}

#' Identify aging-related TFs
#' @param sc_data Seurat object
#' @param final_scores Final TF contribution scores
#' @return List of aging-inducing and aging-resisting TFs
identify_aging_tfs <- function(sc_data, final_scores) {
  sc_data[["TFscore"]] <- CreateAssayObject(counts = final_scores)
  Idents(sc_data) <- "type"
  DefaultAssay(sc_data) <- "TFscore"
  
  deg <- FindMarkers(sc_data,
                     ident.1 = "elderly",
                     ident.2 = "young",
                     min.pct = 0,
                     logfc.threshold = 0,
                     test.use = "wilcox")
  
  # Categorize TFs
  aging_inducing <- deg %>%
    filter(p_val < 0.05, avg_log2FC > 0) %>%
    rownames()
  aging_resisting <- deg %>%
    filter(p_val < 0.05, avg_log2FC < 0) %>%
    rownames()
  
  list(
    inducing = aging_inducing,
    resisting = aging_resisting
  )
}

#' Main function to run HasCATs analysis
#' @param rds_path Path to Seurat RDS file
#' @param reg_path Path to TF regulation file
#' @return List containing final scores and aging TFs
run_hascats <- function(rds_path, reg_path) {
  # Step 1: Load data
  cat("Loading single-cell data...\n")
  sc_data_list <- load_sc_data(rds_path)
  
  # Step 2: Calculate TF expression × activity
  cat("Calculating TF expression × activity...\n")
  tf_expr_activity <- calc_tf_expr_activity(sc_data_list$expr, sc_data_list$tf_act)
  
  # Step 3: Load TF regulation data
  cat("Loading TF regulation data...\n")
  tf_reg <- load_tf_regulation(reg_path)
  
  # Step 4: Calculate edge weights
  cat("Calculating edge weights...\n")
  edge_weights <- calc_edge_weights(tf_reg, sc_data_list$expr, rownames(tf_expr_activity))
  
  # Step 5: Calculate final TF contribution
  cat("Calculating final TF contributions...\n")
  final_scores <- calc_tf_contribution(tf_expr_activity, edge_weights)
  
  # Step 6: Identify aging-related TFs
  cat("Identifying aging-related TFs...\n")
  aging_tfs <- identify_aging_tfs(readRDS(rds_path), final_scores)
  
  list(
    final_scores = final_scores,
    aging_tfs = aging_tfs
  )
}

# Run analysis
rds_path <- "./data/HSC_OldYoung_merge.rds"
reg_path <- "./RoraScript/results/supple/TF and Target for senescence.xls"
results <- run_hascats(rds_path, reg_path)
cat("Analysis complete.")