setwd('Desktop/Research/binchenlab/disease_signature/')

##### Stratify case samples - for one case #####
case <- readRDS('output_files/archs_exp_data/cases/adjacent_normals/Endometrial-Neoplasms_Uterus_GSE146889.rds')
log_case <- log2(case + 1)
var_genes <- apply(log_case, 2, var)
log_case_filtered <- log_case[, var_genes > 0]

# pca
pca <- prcomp(log_case_filtered, center = TRUE, scale. = TRUE)
plot(pca$x[,1:2], col = "darkred", pch = 19, main = "PCA of Samples")
text(pca$x[,1], pca$x[,2], labels = rownames(log_case), pos = 3, cex = 0.6)

# tsne
library(Rtsne)
tsne <- Rtsne(pca$x[,1:10], perplexity = 5)  # tweak perplexity based on sample size
plot(tsne$Y, col = "blue", pch = 19, main = "t-SNE of Samples")
text(tsne$Y[,1], tsne$Y[,2], labels = rownames(case), pos = 3, cex = 0.6)

# umap
library(uwot)
umap <- umap(pca$x[, 1:min(30, ncol(pca$x))], n_neighbors = 15, min_dist = 0.1)
plot(umap, col = "blue", pch = 19, main = "umap of Samples")
text(umap[,1], umap[,2], labels = rownames(case), pos = 3, cex = 0.6)

dist_matrix <- dist(pca$x[,1:10])
hc <- hclust(dist_matrix)
plot(hc, labels = rownames(case), main = "Hierarchical Clustering of Cases")
#####

#### Stratify case samples - All #####
case <- readRDS('output_files/archs_exp_data/cases/adjacent_normals/Endometrial-Neoplasms_Uterus_GSE146889.rds')
controls <- readRDS('output_files/archs_exp_data/controls/adjacent_normals/Endometrial-Neoplasms_Uterus_GSE146889.rds')
gtex <- readRDS('output_files/gtex_exp_data/adjacent_normals/Endometrial-Neoplasms_Uterus_GSE146889.rds')

common_genes <- Reduce(intersect, list(colnames(case), colnames(controls), colnames(gtex)))
case_sub <- case[, common_genes]
controls_sub <- controls[, common_genes]
gtex_sub <- gtex[, common_genes]

combined <- rbind(case_sub, controls_sub, gtex_sub)
log_combined <- log2(combined + 1)
var_genes <- apply(log_combined, 2, var)
log_combined_filtered <- log_combined[, var_genes > 0]
pca <- prcomp(log_combined_filtered, center = TRUE, scale. = TRUE)
group <- c(rep("Case", nrow(case_sub)), rep("Control", nrow(controls_sub)), rep("GTEx", nrow(gtex_sub)))

library(ggplot2)
pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$Group <- factor(group)
pca_df$Sample <- rownames(log_combined_filtered)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) + geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 2.5) + theme_minimal() + ggtitle("PCA of Samples (Case, Control, GTEx)") +
  theme(plot.title = element_text(hjust = 0.5))

# Function
visualize_sample_stratification <- function(case_dir, control_dir, gtex_dir, plot_output_dir, method = c("pca", "tsne")) {
  method <- match.arg(method)
  if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
  if (method == "tsne" && !requireNamespace("Rtsne", quietly = TRUE)) install.packages("Rtsne")
  
  case_files <- list.files(case_dir, full.names = TRUE)
  case_names <- list.files(case_dir)
  
  for (i in seq_along(case_files)) {
    fname <- case_names[i]
    case_file <- case_files[i]
    control_file <- file.path(control_dir, fname)
    gtex_file <- file.path(gtex_dir, fname)
    
    if (!file.exists(control_file) || !file.exists(gtex_file)) {
      cat("Skipping", fname, "- matching control or GTEx not found.\n")
      next
    }
    
    case <- readRDS(case_file)
    controls <- readRDS(control_file)
    gtex <- readRDS(gtex_file)
    
    common_genes <- Reduce(intersect, list(colnames(case), colnames(controls), colnames(gtex)))
    case_sub <- case[, common_genes]
    controls_sub <- controls[, common_genes]
    gtex_sub <- gtex[, common_genes]
    
    log_combined <- log2(rbind(case_sub, controls_sub, gtex_sub) + 1)
    log_case <- log2(case_sub + 1)
    
    var_comb <- apply(log_combined, 2, var)
    var_case <- apply(log_case, 2, var)
    
    log_combined_filtered <- log_combined[, var_comb > 0]
    log_case_filtered <- log_case[, var_case > 0]
    
    # PCA preprocessing
    pca_combined <- prcomp(log_combined_filtered, center = TRUE, scale. = TRUE)
    pca_case <- prcomp(log_case_filtered, center = TRUE, scale. = TRUE)
    num_comb_pcs <- min(30, nrow(pca_combined$x) - 1, ncol(pca_combined$x))
    num_case_pcs <- min(30, nrow(pca_case$x) - 1, ncol(pca_case$x))
    red_input_comb <- pca_combined$x[, 1:num_comb_pcs, drop = FALSE]
    red_input_case <- pca_case$x[, 1:num_case_pcs, drop = FALSE]
    
    set.seed(42)
    # coords_comb <- switch(method, pca = red_input_comb[, 1:2], tsne = Rtsne::Rtsne(red_input_comb, perplexity = min(5, floor(nrow(red_input_comb) / 3)))$Y)
    # coords_case <- switch(method, pca = red_input_case[, 1:2], tsne = Rtsne::Rtsne(red_input_case, perplexity = min(5, floor(nrow(red_input_case) / 3)))$Y)
    coords_comb <- switch(method, pca = red_input_comb[, 1:2], tsne = {
      n_comb <- nrow(red_input_comb)
      max_perp <- floor((n_comb - 1) / 3)
      if (max_perp < 2) {
        warning(paste("Too few samples for t-SNE in combined data:", fname))
        matrix(NA, nrow = n_comb, ncol = 2)} else {
          tsne_perp_comb <- min(5, max_perp)
          Rtsne::Rtsne(red_input_comb, perplexity = tsne_perp_comb)$Y
        }
    }
    )
    
    coords_case <- switch(method, pca = red_input_case[, 1:2], tsne = {
      n_case <- nrow(red_input_case)
      max_perp <- floor((n_case - 1) / 3)
      if (max_perp < 2) {
        warning(paste("Too few samples for t-SNE in case-only data:", fname))
        matrix(NA, nrow = n_case, ncol = 2)} else {
          tsne_perp_case <- min(5, max_perp)
          Rtsne::Rtsne(red_input_case, perplexity = tsne_perp_case)$Y
        }
    }
    )
    
    group <- c(rep("Case", nrow(case_sub)), rep("Control", nrow(controls_sub)), rep("GTEx", nrow(gtex_sub)))
    sample_labels_comb <- rownames(log_combined_filtered)
    sample_labels_case <- rownames(log_case_filtered)
    
    # Hierarchical clustering on case
    pc_count <- min(10, ncol(pca_case$x))
    dist_matrix <- dist(pca_case$x[, 1:pc_count, drop = FALSE])
    hc <- hclust(dist_matrix)
    out_file <- file.path(plot_output_dir, paste0(tools::file_path_sans_ext(fname), "_", method, ".pdf"))
    pdf(out_file)
    
    # Plot 1: All samples
    colors <- c("Case" = "darkred", "Control" = "grey", "GTEx" = "green")
    point_colors <- colors[group]
    # plot(coords_comb, col = point_colors, pch = 19, main = paste(toupper(method), "- All Samples"), xlab = paste0(toupper(method), "1"), ylab = paste0(toupper(method), "2"))
    # text(coords_comb[,1], coords_comb[,2], labels = sample_labels_comb, pos = 3, cex = 0.6)
    # legend("topright", legend = names(colors), col = colors, pch = 19, cex = 0.8, box.lwd = 0.5)
    if (all(is.finite(coords_comb))) {
      plot(coords_comb, col = point_colors, pch = 19,
           main = paste(toupper(method), "- All Samples"),
           xlab = paste0(toupper(method), "1"), ylab = paste0(toupper(method), "2"))
      text(coords_comb[,1], coords_comb[,2], labels = sample_labels_comb, pos = 3, cex = 0.6)
      legend("topright", legend = names(colors), col = colors, pch = 19, cex = 0.8, box.lwd = 0.5)
    } else {
      plot.new()
      title(main = paste("t-SNE Skipped (Too Few Samples) - All Samples"))
    }
    
    # Plot 2: Case only
    # plot(coords_case, col = "darkred", pch = 19, main = paste(toupper(method), "- Case Only"), xlab = paste0(toupper(method), "1"), ylab = paste0(toupper(method), "2"))
    # text(coords_case[,1], coords_case[,2], labels = sample_labels_case, pos = 3, cex = 0.6)
    if (all(is.finite(coords_case))) {
      plot(coords_case, col = "darkred", pch = 19, main = paste(toupper(method), "- Case Only"), xlab = paste0(toupper(method), "1"), ylab = paste0(toupper(method), "2"))
      text(coords_case[,1], coords_case[,2], labels = sample_labels_case, pos = 3, cex = 0.6)
    } else {
      plot.new()
      title(main = paste("t-SNE Skipped (Too Few Samples) - Case Only"))
    }
    
    # Plot 3: Clustering
    plot(hc, labels = sample_labels_case, main = paste("Hierarchical Clustering - Case Only"))
    dev.off()
    cat("Saved:", out_file, "\n")
  }
}

# Normal
case_dir <- "output_files/archs_exp_data/cases/normals/"
control_dir <- "output_files/archs_exp_data/controls/normals/"
gtex_dir <- "output_files/gtex_exp_data/normals/"

# PCA
plot_output_dir <- "output_files/plots/sample_stratification/normals/"
visualize_sample_stratification(case_dir, control_dir, gtex_dir, plot_output_dir, method = "pca")

# TSNE
plot_output_dir <- "output_files/plots/sample_stratification/normals/"
visualize_sample_stratification(case_dir, control_dir, gtex_dir, plot_output_dir, method = "tsne")

# Adjacent normal
case_dir <- "output_files/archs_exp_data/cases/adjacent_normals/"
control_dir <- "output_files/archs_exp_data/controls/adjacent_normals/"
gtex_dir <- "output_files/gtex_exp_data/adjacent_normals/"

# PCA
plot_output_dir <- "output_files/plots/sample_stratification/adjacent_normals/"
visualize_sample_stratification(case_dir, control_dir, gtex_dir, plot_output_dir, method = "pca")

# TSNE
plot_output_dir <- "output_files/plots/sample_stratification/adjacent_normals/"
visualize_sample_stratification(case_dir, control_dir, gtex_dir, plot_output_dir, method = "tsne")
#####

#####
case <- readRDS('output_files/archs_exp_data/cases/adjacent_normals/Endometrial-Neoplasms_Uterus_GSE146889.rds')
log_case <- log2(case + 1)
var_genes <- apply(log_case, 2, var)
log_case_filtered <- log_case[, var_genes > 0]

# pca
pca <- prcomp(log_case_filtered, center = TRUE, scale. = TRUE)
# plot(pca$x[,1:2], col = "darkred", pch = 19, main = "PCA of Samples")
# text(pca$x[,1], pca$x[,2], labels = rownames(log_case), pos = 3, cex = 0.6)

# Gap Statistic (to test optimal number of clusters)
library(cluster)
set.seed(42)
gap <- clusGap(pca$x[, 1:10], FUN = kmeans, K.max = 6, B = 100)  # B = bootstrap samples
plot(gap, main = "Gap Statistic (Optimal Clusters in Case Samples)")
best_k <- which.max(gap$Tab[, "gap"])

k_range <- 2:6
avg_sil <- sapply(k_range, function(k) {
  km <- kmeans(pca$x[, 1:10], centers = k, nstart = 10)
  sil <- silhouette(km$cluster, dist(pca$x[, 1:10]))
  mean(sil[, 3])
})
plot(k_range, avg_sil, type = "b", xlab = "Number of Clusters (k)", ylab = "Average Silhouette Width", main = "Silhouette Analysis - Case Samples")
#####

##### Loop for all the cases #####
cases_file_names <- list.files('output_files/archs_exp_data/cases/adjacent_normals/')
control_files_names <- list.files('output_files/archs_exp_data/controls/adjacent_normals/')
final_cases_files <- intersect(cases_file_names, control_files_names)

case_file_path <- 'output_files/archs_exp_data/cases/adjacent_normals/'
length(final_cases_files)

library(cluster)

# # # For normals
# # case_file_path <- "output_files/archs_exp_data/cases/normals/"
# # final_cases_files <- intersect(list.files("output_files/archs_exp_data/cases/normals/"), list.files("output_files/archs_exp_data/controls/normals/"))
# 
# # # For adjacent normals
# # case_file_path <- "output_files/archs_exp_data/cases/adjacent_normals/"
# # final_cases_files <- intersect(list.files("output_files/archs_exp_data/cases/adjacent_normals/"), list.files("output_files/archs_exp_data/controls/adjacent_normals/"))

summary_df <- data.frame(File = character(), Num_Samples = integer(), Gap_k = integer(), Silhouette_k = integer(), Sample_Clusters = character(), stringsAsFactors = FALSE)

for (file_name in final_cases_files) {
  file_path <- file.path(case_file_path, file_name)
  cat("Processing:", file_name, "\n")
  
  # Load and preprocess
  case <- readRDS(file_path)
  log_case <- log2(case + 1)
  var_genes <- apply(log_case, 2, var)
  log_case_filtered <- log_case[, var_genes > 0]
  
  num_samples <- nrow(log_case_filtered)
  
  if (num_samples < 3) {
    warning(paste("Too few samples to cluster:", file_name))
    summary_df <- rbind(summary_df, data.frame(
      File = file_name,
      Num_Samples = num_samples,
      Gap_k = NA,
      Silhouette_k = NA,
      Sample_Clusters = NA
    ))
    next
  }
  
  # PCA
  pca <- prcomp(log_case_filtered, center = TRUE, scale. = TRUE)
  pc_mat <- pca$x[, 1:min(10, ncol(pca$x))]
  
  K.max <- min(6, nrow(pc_mat) - 1)
  if (K.max < 2) {
    warning(paste("Too few samples for multiple clusters:", file_name))
    summary_df <- rbind(summary_df, data.frame(
      File = file_name,
      Num_Samples = num_samples,
      Gap_k = NA,
      Silhouette_k = NA,
      Sample_Clusters = NA
    ))
    next
  }
  
  # Gap Statistic
  set.seed(42)
  gap <- clusGap(pc_mat, FUN = kmeans, K.max = K.max, B = 100)
  gap_k <- which.max(gap$Tab[, "gap"])
  
  # Silhouette
  sil_k_range <- 2:K.max
  avg_sil <- sapply(sil_k_range, function(k) {
    km <- kmeans(pc_mat, centers = k, nstart = 10)
    sil <- silhouette(km$cluster, dist(pc_mat))
    mean(sil[, 3])
  })
  sil_k <- sil_k_range[which.max(avg_sil)]
  
  # Run k-means with best gap_k and generate sample:cluster string
  km_assign <- kmeans(pc_mat, centers = gap_k, nstart = 10)
  sample_ids <- rownames(pc_mat)
  cluster_assignments <- paste0(sample_ids, ":", km_assign$cluster)
  sample_clusters_str <- paste(cluster_assignments, collapse = ", ")
  
  # Append to summary_df
  summary_df <- rbind(summary_df, data.frame(
    File = file_name,
    Num_Samples = num_samples,
    Gap_k = gap_k,
    Silhouette_k = sil_k,
    Sample_Clusters = sample_clusters_str,
    stringsAsFactors = FALSE
  ))
}
write.csv(summary_df, "output_files/plots/sample_stratification/cluster_validation_summary_adjacent_normal.csv", row.names = FALSE)
write.csv(summary_df, "output_files/plots/sample_stratification/cluster_validation_summary_normal.csv", row.names = FALSE)
#####

##### Relation between clustering metrics with correlation #####
# # # For normal
# # corr <- read.csv('output_files/plots/corr_logFC_geo_normals/corr_logFC_summary_geo_normals.csv')
# # clustering <- read.csv('output_files/plots/sample_stratification/cluster_validation_summary_normal.csv')
# 
# # # For adjacent normal
# # corr <- read.csv('output_files/plots/corr_logFC_geo_adjacent_normals/corr_logFC_summary_geo_adjacent_normals.csv')
# # clustering <- read.csv('output_files/plots/sample_stratification/cluster_validation_summary_adjacent_normal.csv')

df_merged <- merge(corr, clustering, by.x = 'Filename', by.y = 'File')
head(df_merged, 3)

df_merged$Pearson_Correlation <- as.numeric(df_merged$Pearson_Correlation)
df_merged$N_Genes <- as.numeric(df_merged$N_Genes)
df_merged$Num_Samples <- as.numeric(df_merged$Num_Samples)
df_merged$Gap_k <- as.numeric(df_merged$Gap_k)
df_merged$Silhouette_k <- as.numeric(df_merged$Silhouette_k)
numeric_vars <- df_merged[, c("Pearson_Correlation", "N_Genes", "Num_Samples", "Gap_k", "Silhouette_k")]
cor_results <- cor(numeric_vars, use = "complete.obs")
cor_results["Pearson_Correlation", ]

library(ggplot2)
cor_data <- data.frame(Variable = names(cor_results["Pearson_Correlation", -1]), Correlation = cor_results["Pearson_Correlation", -1])
ggplot(cor_data, aes(x = reorder(Variable, Correlation), y = Correlation)) + geom_bar(stat = "identity") +
  coord_flip() + theme_minimal() + labs(title = "Correlation with Pearson Correlation - Normal", x = "", y = "Pearson Correlation")
ggplot(cor_data, aes(x = reorder(Variable, Correlation), y = Correlation)) + geom_bar(stat = "identity") +
  coord_flip() + theme_minimal() + labs(title = "Correlation with Pearson Correlation - Adjacent Normal", x = "", y = "Pearson Correlation")


library(reshape2)
cor_results <- cor(numeric_vars, use = "complete.obs")
cor_results[lower.tri(cor_results)] <- NA
melted_cor <- melt(cor_results, na.rm = TRUE)

ggplot(melted_cor, aes(x = Var1, y = Var2, fill = value)) + geom_tile() +
  scale_fill_gradient2(low = "grey", high = "darkred", limit = c(0,1)) + 
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  theme_minimal() + labs(title = "Correlation Heatmap", fill = "Pearson\nCorrelation")
#####