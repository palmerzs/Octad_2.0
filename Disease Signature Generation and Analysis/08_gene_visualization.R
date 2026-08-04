setwd('Desktop/Research/binchenlab/disease_signature/')

##### Required libraries #####
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(readr)
library(ggplot2)


library(metafor)
#####

##### Explore how can see gene or gene-set enrichment #####
sig_path <- 'output_files/signature_for_drug_discovery/signatures_new/'
disease <- list.files(sig_path)
# head(disease); length(unique(disease))

read_signature <- function(file, disease_name) {
  df <- read.csv(file)
  df <- df %>% dplyr::select(identifier, log2FoldChange, padj) %>% mutate(disease = disease_name)
  return(df)
}
disease_dirs <- list.files(sig_path, full.names = TRUE)
all_results <- lapply(disease_dirs, function(d_dir) {
  files <- list.files(d_dir, full.names = TRUE)
  disease_name <- basename(d_dir)
  lapply(files, function(f) read_signature(f, disease_name))
})
all_results <- do.call(rbind, unlist(all_results, recursive = FALSE))
# all_results <- all_results %>% mutate(log2FoldChange = ifelse(padj < 0.05, log2FoldChange, NA))
gene_disease_matrix <- all_results %>% group_by(identifier, disease) %>% summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = disease, values_from = log2FoldChange) %>% mutate(across(-identifier, ~replace_na(., 0)))
# head(gene_disease_matrix)

mapping <- read.csv('../gene_info_biomart.csv')
gene_disease_matrix <- gene_disease_matrix %>% 
  left_join(mapping %>% group_by(gene_id) %>% summarise(gene_name = first(gene_name), .groups = "drop"), by = c("identifier" = "gene_id")) %>%
  relocate(gene_name, .after = identifier)
# saveRDS(gene_disease_matrix, file = 'output_files/signature_for_drug_discovery/gene_enrichment/gene_disease_matrix.rds')

gene_disease_matrix <- readRDS('output_files/signature_for_drug_discovery/gene_enrichment/gene_disease_matrix.rds')
dim(gene_disease_matrix)
head(gene_disease_matrix)[1:5]

gene_profile_plot <- function(data, gene_of_interest, fc_threshold = 1, type = c("bar", "trend"), order = TRUE) {
  type <- match.arg(type)
  df <- data %>% filter(gene_name == gene_of_interest)
  if (nrow(df) == 0) {
    message("Gene named ", '"', gene_of_interest, '"', ' not found')
    return(invisible(NULL))
  }
  disease_order <- names(df)[sapply(df, is.numeric)]
  df_long <- df %>% select(gene_name, all_of(disease_order)) %>% pivot_longer(cols = -gene_name, names_to = "disease", values_to = "log2FoldChange") %>% filter(abs(log2FoldChange) >= fc_threshold)
  if (order) {df_long <- df_long %>% mutate(disease = reorder(disease, log2FoldChange))} else {df_long <- df_long %>% mutate(disease = factor(disease, levels = disease_order))}
  if (type == "bar") {
    p <- ggplot(df_long, aes(x = disease, y = log2FoldChange, fill = log2FoldChange > 0)) + geom_col() +
      scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") + coord_flip() + theme_minimal(base_size = 10) +
      theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold")) +
      labs(title = paste("FoldChange across diseases for:", gene_of_interest), x = "Disease", y = "log2FoldChange")
  } else if (type == "trend") {
    p <- ggplot(df_long, aes(x = disease, y = log2FoldChange, group = 1)) + geom_line(color = "steelblue", linewidth = 0.3) + geom_point(aes(color = log2FoldChange > 0), size = 1) +
      scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") + theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 8), plot.title = element_text(face = "bold")) + coord_flip() +
      labs(title = paste("Trend across Diseases for:", gene_of_interest), x = "Disease", y = "log2FoldChange")
  }
  print(p); invisible(p)
}
gene_profile_plot(gene_disease_matrix, "MYCN", fc_threshold = 3, order = F, type = "bar")
# gene_profile_plot(gene_disease_matrix, "HBB", fc_threshold = 3, order = F, type = "bar")
gene_profile_plot(gene_disease_matrix, "HBB", fc_threshold = 3, order = T, type = "trend")


# gene_profile_plot(gene_disease_matrix, "GLP1R", fc_threshold = 0.1, order = T, type = "bar")
# gene_profile_plot(gene_disease_matrix, "GLP1R", fc_threshold = 0.1, order = T, type = "trend")


gene_profile_heatmap <- function(data, gene_of_interest, mode = c("pattern", "list"), fc_threshold = 0) {
  mode <- match.arg(mode)
  if (mode == "pattern") {message("Using pattern mode: matching genes with pattern '", gene_of_interest, "'")
    df <- data %>% filter(grepl(gene_of_interest, gene_name, ignore.case = TRUE))
  } else if (mode == "list") {message("Using list mode for genes: ", paste(gene_of_interest, collapse = ", "))
    df <- data %>% filter(gene_name %in% gene_of_interest)}
  if (nrow(df) == 0) {stop("No genes found matching the provided pattern or list.")}
  
  df_long <- df %>% select(gene_name, where(is.numeric)) %>%  pivot_longer(cols = -gene_name, names_to = "disease", values_to = "log2FoldChange") %>%
    filter(abs(log2FoldChange) >= fc_threshold) %>% mutate(log2FoldChange = pmax(pmin(log2FoldChange, 5), -5))
  if (nrow(df_long) == 0) stop("No values pass the fold change threshold.")
  
  mat <- df_long %>% pivot_wider(names_from = disease, values_from = log2FoldChange) %>% as.data.frame()
  rownames(mat) <- mat$gene_name
  mat <- as.matrix(mat[, -1, drop = FALSE])
  mat[is.na(mat)] <- 0
  disease_order <- colnames(mat)[hclust(dist(t(mat)))$order]
  df_long$disease <- factor(df_long$disease, levels = disease_order)
  
  p <- ggplot(df_long, aes(x = gene_name, y = disease, fill = log2FoldChange)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, limits = c(-5, 5), name = "log2FoldChange") + theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 8), panel.grid = element_blank()) +
    labs(title = paste0("Fold Change Heatmap"), x = "Gene", y = "Disease")
  print(p); invisible(p)
}


# gene_profile_heatmap(gene_disease_matrix, gene_of_interest = 'AIFM', mode = "pattern", fc_threshold = 1)
genes <- c('AIFM2', 'EPHX1', 'EPHX2', 'PLA2G6', 'FLT1')
gene_profile_heatmap(gene_disease_matrix, gene_of_interest = genes, mode = "list", fc_threshold = 1)


# genes <- c('GLP')
# gene_profile_heatmap(gene_disease_matrix, gene_of_interest = genes, mode = "pattern", fc_threshold = 0.1)


test <- gene_disease_matrix %>% filter(gene_name %in% genes)
write.csv(test, file = '../../../gene_disease_matrix.csv')




#####
gene_disease_matrix <- readRDS('output_files/signature_for_drug_discovery/gene_enrichment/gene_disease_matrix.rds')
dim(gene_disease_matrix)
head(gene_disease_matrix)[1:5]

metadata <- read.csv('output_files/signature_for_drug_discovery/sig_meta_final.csv')
length(unique(metadata$Disease))
length(unique(metadata$OrganRegion))
(unique(metadata$OrganRegion))




library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)

gene_disease_matrix <- readRDS('output_files/signature_for_drug_discovery/gene_enrichment/gene_disease_matrix.rds')

metadata <- read.csv('output_files/signature_for_drug_discovery/sig_meta_final.csv')

harmonize_organ <- function(x) {
  x_chr <- as.character(x)
  x_clean <- str_trim(x_chr)
  x_clean <- str_replace_all(x_clean, "\\s+", "_")
  x_key <- str_to_lower(x_clean)
  case_when(is.na(x_chr) | x_key %in% c("", "na", "n/a", "ua", "unknown") ~ "Unknown", 
            x_key == "adrenal_gland" ~ "Adrenal_gland", x_key %in% c("lymphatic_tissue", "lymphoid_tissue") ~ "Lymphoid_tissue", x_key == "head_and_neck" ~ "Head_and_neck",
            x_key == "lining_of_body_cavities" ~ "Lining_of_body_cavities", TRUE ~ str_to_sentence(x_key)
  )
}

metadata <- metadata %>% mutate(OrganRegion = harmonize_organ(OrganRegion))
gene_profile_heatmap <- function(data, metadata, gene_of_interest, mode = c("pattern", "list"), fc_threshold = 0, cap = 5) {
  mode <- match.arg(mode)
  
  df <- if (mode == "pattern") {data %>% filter(grepl(gene_of_interest, gene_name, ignore.case = TRUE))
  } else {data %>% filter(gene_name %in% gene_of_interest)}
  disease_annotation <- metadata %>% select(Disease, OrganRegion) %>% mutate(OrganRegion = harmonize_organ(OrganRegion)) %>% distinct()
  df_long <- df %>% select(identifier, gene_name, where(is.numeric)) %>% pivot_longer(-c(identifier, gene_name), names_to = "Disease", values_to = "log2FoldChange") %>%
    left_join(disease_annotation, by = "Disease", relationship = "many-to-many") %>% 
    mutate(OrganRegion = harmonize_organ(OrganRegion), log2FoldChange = pmax(pmin(log2FoldChange, cap), -cap), Disease_Organ_ID = paste(Disease, OrganRegion, sep = "__")) %>%
    filter(abs(log2FoldChange) >= fc_threshold)
  mat_df <- df_long %>% group_by(Disease_Organ_ID, Disease, OrganRegion, gene_name) %>% summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = gene_name, values_from = log2FoldChange, values_fill = 0) %>% arrange(OrganRegion, Disease)
  
  organ_vec <- mat_df$OrganRegion
  disease_labels <- mat_df$Disease
  mat <- mat_df %>% select(-Disease, -OrganRegion) %>% as.data.frame()
  rownames(mat) <- mat$Disease_Organ_ID
  mat <- as.matrix(mat[, -1, drop = FALSE])
  
  organ_midpoints <- mat_df %>% mutate(row_id = row_number()) %>% group_by(OrganRegion) %>% summarise(midpoint = round(mean(row_id)), .groups = "drop")
  right_anno <- rowAnnotation(Organ = anno_mark(at = organ_midpoints$midpoint, labels = organ_midpoints$OrganRegion, labels_gp = gpar(fontsize = 8, fontface = "bold")), width = unit(3.5, "cm"))
  
  ht <- Heatmap(mat, name = "log2FC", col = circlize::colorRamp2(c(-cap, 0, cap), c("steelblue", "white", "firebrick")), cluster_rows = FALSE, cluster_columns = FALSE,
                show_row_names = TRUE, row_labels = disease_labels, row_names_side = "left", row_names_gp = gpar(fontsize = 6), show_column_names = TRUE, right_annotation = right_anno)
  
  draw(ht)
  invisible(ht)
}

png("../../../GLP_heatmap.png", width = 4500, height = 9000, res = 600)
gene_profile_heatmap(data = gene_disease_matrix, metadata = metadata, gene_of_interest = "GLP", mode = "pattern", fc_threshold = 1, cap = 5)
dev.off()


