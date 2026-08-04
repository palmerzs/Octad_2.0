setwd('Desktop/Research/binchenlab/disease_signature/')

##### Required libraries #####
library(dplyr)
library(readxl)
library(tibble)
library(purrr)
#####

##### Check disease frequency which has matched controls #####
get_matched_disease_info <- function(metadata, label) {
  get_disease_info <- function(df, control_type, min_n = 2) {
    controls <- df %>% filter(Sample_type == control_type) %>% group_by(GSE_ID, OrganRegion) %>% 
      summarise(n_controls   = n(), controls_gsm = paste(GSM_ID, collapse = ","), .groups = "drop")
    cases <- df %>% filter(Sample_type == "Case") %>% group_by(GSE_ID, Disease, OrganRegion) %>% 
      summarise(n_cases   = n(), cases_gsm = paste(GSM_ID, collapse = ","), .groups = "drop")
    matched <- cases %>% inner_join(controls, by = c("GSE_ID", "OrganRegion")) %>% filter(n_cases >= min_n, n_controls >= min_n)
    disease_counts <- matched %>% group_by(Disease, OrganRegion) %>% summarise(Experiments = n_distinct(GSE_ID), .groups = "drop")
    list(disease_data = matched, disease_counts = disease_counts)
  }
  
  info_healthy    <- get_disease_info(metadata, "Healthy")
  info_normal     <- get_disease_info(metadata, "Normal")
  info_adj_normal <- get_disease_info(metadata, "AdjacentNormal")
  
  matched_disease_info <- bind_rows(info_healthy$disease_data    %>% mutate(Control_type = "Healthy"),
                                    info_normal$disease_data     %>% mutate(Control_type = "Normal"),
                                    info_adj_normal$disease_data %>% mutate(Control_type = "AdjacentNormal")) %>%
    transmute(Disease, GSE_ID, OrganRegion, Control_type, n_cases, n_controls, cases_gsm, controls_gsm) %>% arrange(Disease, GSE_ID, OrganRegion, Control_type)
  # return(matched_disease_info)
  write.csv(matched_disease_info, paste0("data/matched_disease_info_", label, ".csv"), row.names = FALSE)
}
load('data/metadata_refined.RData')
get_matched_disease_info(metadata_invitro, "invitro"); get_matched_disease_info(metadata_invivo, "invivo"); get_matched_disease_info(metadata_exvivo, "exvivo")

# Matched disease overlap frequency
matched_disease_invitro  <- read.csv('data/matched_disease_info_invitro.csv')
matched_disease_invivo   <- read.csv('data/matched_disease_info_invivo.csv')
matched_disease_exvivo   <- read.csv('data/matched_disease_info_exvivo.csv')

diseases_invitro <- unique(matched_disease_invitro$Disease); diseases_invivo  <- unique(matched_disease_invivo$Disease); diseases_exvivo  <- unique(matched_disease_exvivo$Disease)

all_unique_diseases <- unique(c(diseases_invitro, diseases_invivo, diseases_exvivo))
common_all_three    <- Reduce(intersect, list(diseases_invitro, diseases_invivo, diseases_exvivo))
disease_overlap_summary <- list(invitro = diseases_invitro, invivo = diseases_invivo, exvivo = diseases_exvivo, all_unique = all_unique_diseases, common_all = common_all_three)

# saveRDS(disease_overlap_summary, file = "data/matched_disease_overlap_info_without_blood.rds")
# saveRDS(disease_overlap_summary, file = "data/matched_disease_overlap_info_with_blood.rds")
#####

##### Extract raw reads counts from ARCHS #####
library(rhdf5)
library(dplyr)

sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
# function to extract expression
get_expression_from_archs <- function(mdm, label, expression_file_path, gsm_in_file, gene_id) {
  output_path <- file.path("data/exp_data", label)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path("data/exp_data", paste0("exp_extraction_log_", label, ".txt"))
  log_msg <- function(msg, log_file) {message(msg); cat(msg, "\n", file = log_file, append = TRUE)}
  if (file.exists(log_file)) file.remove(log_file)
  
  h5_file <- H5Fopen(expression_file_path)
  for (i in seq_len(nrow(mdm))) {
    disease      <- mdm$Disease[i]
    gse_id       <- mdm$GSE_ID[i]
    organ        <- mdm$OrganRegion[i]
    control_type <- mdm$Control_type[i]
    
    filename <- paste0(sanitize_disease(disease), "_", sanitize_simple(organ), "_", sanitize_simple(control_type), "_", sanitize_simple(gse_id), ".rds")
    output_file <- file.path(output_path, filename)
    if (file.exists(output_file)) {
      log_msg(sprintf("Skipping existing file: %s", filename), log_file)
      next
    }
    
    # collect GSMs
    case_ids    <- unlist(strsplit(mdm$cases_gsm[i], ","))
    control_ids <- unlist(strsplit(mdm$controls_gsm[i], ","))
    case_idx    <- match(case_ids, gsm_in_file)
    control_idx <- match(control_ids, gsm_in_file)
    
    valid_case_idx    <- which(!is.na(case_idx))
    valid_control_idx <- which(!is.na(control_idx))
    
    all_idx <- c(case_idx[valid_case_idx], control_idx[valid_control_idx])
    all_ids <- c(case_ids[valid_case_idx], control_ids[valid_control_idx])
    
    # handle errors
    if (length(valid_case_idx) == 0 & length(valid_control_idx) == 0) {
      log_msg(sprintf("No valid GSMs for: Disease=\"%s\" / Organ=%s / GSE_ID=%s | Missing: BOTH cases and controls", disease, organ, gse_id), log_file)
      next
    }
    if (length(valid_case_idx) == 0) {
      log_msg(sprintf("No valid CASE GSMs for: Disease=\"%s\" / Organ=%s / GSE_ID=%s", disease, organ, gse_id), log_file)
      next
    }
    if (length(valid_control_idx) == 0) {
      log_msg(sprintf("No valid CONTROL GSMs for: Disease=\"%s\" / Organ=%s / GSE_ID=%s", disease, organ, gse_id), log_file)
      next
    }
    
    # read expression
    expr <- tryCatch({
      h5read(expression_file_path, "/data/expression", index = list(all_idx, NULL))
    }, error = function(e) {
      log_msg(sprintf("Error reading: Disease=\"%s\" / Organ=%s / GSE_ID=%s — %s", disease, organ, gse_id, e$message), log_file)
      return(NULL)
    })
    if (is.null(expr)) next
    
    rownames(expr) <- all_ids; colnames(expr) <- gene_id
    
    # sample metadata
    meta <- data.frame(GSM_ID = all_ids, Group  = ifelse(all_ids %in% case_ids, "Case", "Control"), 
                       Disease = disease, OrganRegion = organ, Control_type = control_type, GSE_ID = gse_id, stringsAsFactors = FALSE)
    saveRDS(list(expr = expr, metadata = meta), output_file)
    log_msg(sprintf("Saved: %s", output_file), log_file)
  }
  H5Fclose(h5_file)
}

expression_file_path <- "data/human_gene_v2.5.h5"
gsm_in_file <- h5read(expression_file_path, "/meta/samples/geo_accession")
gene_id <- h5read(expression_file_path, "/meta/genes/ensembl_gene")

# run extraction
mdm_invitro <- read.csv("data/matched_disease_info_invitro.csv")
get_expression_from_archs(mdm_invitro, "invitro", expression_file_path, gsm_in_file, gene_id)

mdm_invivo  <- read.csv("data/matched_disease_info_invivo.csv")
get_expression_from_archs(mdm_invivo,  "invivo",  expression_file_path, gsm_in_file, gene_id)

mdm_exvivo  <- read.csv("data/matched_disease_info_exvivo.csv")
get_expression_from_archs(mdm_exvivo,  "exvivo",  expression_file_path, gsm_in_file, gene_id)
#####

##### Differential Gene Expression Analysis #####
library(edgeR)
library(RUVSeq)

# Function: remove lowly expressed genes
remLowExpr <- function(counts, sample_type) {
  x <- DGEList(counts = round(counts), group = sample_type$sample_type)
  cpm_x <- cpm(x)
  keep.exprs <- rowSums(cpm_x > 1) >= min(table(sample_type$sample_type))
  return(keep.exprs)
}
# Function: run DE and save DE results
do_DGE <- function(rds_file, out_dir = "results", use_ruv = FALSE, k = 1, n_topGenes = 2000) {
  data <- readRDS(rds_file)
  expr <- data$expr
  meta <- data$metadata
  case_counts <- expr[meta$Group == "Case", , drop = FALSE]
  control_counts <- expr[meta$Group == "Control", , drop = FALSE]
  
  if (nrow(case_counts) == 0 | nrow(control_counts) == 0) {
    message("Skipping ", rds_file, " (missing case or control samples)")
    return(NULL)
  }
  
  # edgeR
  all_counts <- rbind(case_counts, control_counts)
  sample_type <- factor(c(rep("Case", nrow(case_counts)), rep("Control", nrow(control_counts))), levels = c("Control", "Case"))
  sample_labels <- data.frame(sample_type = sample_type)
  rownames(sample_labels) <- rownames(all_counts)
  
  keep_genes <- remLowExpr(t(all_counts), sample_labels)
  filtered_counts <- all_counts[, keep_genes, drop = FALSE]

  set <- newSeqExpressionSet(as.matrix(t(filtered_counts)), phenoData = AnnotatedDataFrame(sample_labels))
  design <- model.matrix(~sample_type, data = pData(set))
  dge <- DGEList(counts = counts(set), group = set$sample_type)
  
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = 2)
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue, method = 'BH')
  res$identifier <- rownames(res)
  
  # Save
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(rds_file)), ".csv"))
  write.csv(res, out_file, row.names = FALSE)
  return(sprintf("SUCCESS: Saved DE results -> %s", out_file))
}

# # Example: Execute DE for one file
# do_DGE("data/exp_data/exvivo/alpha-1-Antitrypsin-Deficiency_Lung_Normal_GSE194313.rds", out_dir = "output_files/dge_results/exvivo")
# Execute: DE for all files in a directory
do_DGE_for_directory <- function(input_dir, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(output_dir, "00_dge_log.txt")
  if (file.exists(log_file)) file.remove(log_file)
  
  rds_files <- list.files(input_dir, full.names = TRUE, pattern = "\\.rds$")
  for (f in rds_files) {
    tryCatch({
      msg <- do_DGE(f, out_dir = output_dir); message(msg); cat(msg, "\n", file = log_file, append = TRUE)
    }, error = function(e) {
      msg <- sprintf("ERROR: %s | %s", f, e$message); message(msg); cat(msg, "\n", file = log_file, append = TRUE)
    })
  }
  message(sprintf("Log written to: %s", log_file))
}

# InVitro
do_DGE_for_directory("data/exp_data/invitro", "output_files/dge_results/invitro")
# InVivo
do_DGE_for_directory("data/exp_data/invivo",  "output_files/dge_results/invivo")
# ExVivo
do_DGE_for_directory("data/exp_data/exvivo", "output_files/dge_results/exvivo")
#####

##### QC I - Data science approach #####
library(dplyr)
summarize_dge <- function(file, exp_type) {
  res <- read.csv(file)
  n_genes <- nrow(res); sig <- res %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5); sig_genes <- nrow(sig)
  prop_pval <- mean(res$pvalue < 0.01, na.rm = TRUE); mean_abs_logFC <- mean(abs(res$log2FoldChange), na.rm = TRUE)
  snr <- if (sig_genes > 0) mean(abs(sig$log2FoldChange), na.rm = TRUE) else 0
  up <- sum(sig$log2FoldChange > 0.5); down <- sum(sig$log2FoldChange < -0.5)
  
  tibble(file = basename(file), exp_type = exp_type, n_genes = n_genes, sig_genes = sig_genes,
         prop_pval = prop_pval, mean_abs_logFC = mean_abs_logFC, snr = snr, up = up, down = down)
}

exp_dirs <- c("invitro", "invivo", "exvivo")
qc_all <- list()
for (exp in exp_dirs) {
  dge_files <- list.files(file.path("output_files/dge_results", exp), full.names = TRUE, pattern = "\\.csv$")
  qc_all[[exp]] <- bind_rows(lapply(dge_files, summarize_dge, exp_type = exp))
}
qc_summary <- bind_rows(qc_all) %>% mutate(Disease = sub("_.*", "", file))
head(qc_summary); write.csv(qc_summary, "output_files/dge_results/00_qc_summary.csv", row.names = FALSE)

# hist(qc_summary$sig_genes)
# hist(qc_summary$prop_pval)
# hist(qc_summary$mean_abs_logFC)
# hist(qc_summary$snr)

high_quality <- qc_summary %>% filter(sig_genes >= 10 & sig_genes <= 10000, prop_pval > 0.01 & prop_pval < 0.5, 
                                      mean_abs_logFC > 0.1 & mean_abs_logFC < 3, snr > 0.8 & snr < 5, up > 0 & down > 0)
head(high_quality); write.csv(high_quality, "output_files/dge_results/01_high_quality_DE.csv", row.names = FALSE)

low_quality <- anti_join(qc_summary, high_quality)
head(low_quality); write.csv(low_quality, "output_files/dge_results/02_low_quality_DE.csv", row.names = FALSE)
#####

##### Correlate high quality signature #####
library(dplyr)
library(purrr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(patchwork)

high_quality <- read.csv("output_files/dge_results/01_high_quality_DE.csv")

plot_heatmap <- function(cor_mat, disease, out_file) {
  hc <- hclust(as.dist(1 - cor_mat)); order <- hc$order; cor_mat <- cor_mat[order, order]
  cor_long <- reshape2::melt(cor_mat); names(cor_long) <- c("Dataset1", "Dataset2", "Correlation")
  cor_long$Dataset1 <- factor(cor_long$Dataset1, levels = rownames(cor_mat)); cor_long$Dataset2 <- factor(cor_long$Dataset2, levels = colnames(cor_mat))

  heat <- ggplot(cor_long, aes(Dataset1, Dataset2, fill = Correlation)) + geom_tile(color = "white") + geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred", midpoint = 0, limit = c(-1, 1), name = "Pearson\nCorrelation") +
    theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
                                          axis.text.y = element_text(size = 6), axis.title = element_blank(),
                                          panel.grid = element_blank(), plot.title = element_blank())
  dend <- ggdendrogram(hc, rotate = FALSE) + theme_void()
  final_plot <- dend / heat + plot_layout(heights = c(0.7, 5))
  final_plot <- final_plot + plot_annotation(title = paste("Pearson Correlation (log2FC):", disease), theme = theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold")))
  ggsave(out_file, final_plot, width = 8, height = 8, dpi = 600)
}
correlate_diseases <- function(experiment, de_res_dir = file.path("output_files/dge_results", experiment), 
                               out_dir = file.path("output_files/cor_results_matched", experiment)) {
  df <- high_quality %>% filter(exp_type == experiment)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (disease in unique(df$Disease)) {
    files <- df %>% filter(Disease == disease) %>% pull(file)
    if (length(files) < 2) next
    fc_matrix <- map_dfr(files, ~ read.csv(file.path(de_res_dir, .x)) %>% select(identifier, log2FoldChange) %>% mutate(file = .x)) %>%
      pivot_wider(names_from = file, values_from = log2FoldChange) %>% tibble::column_to_rownames("identifier") %>% as.matrix()
    cor_mat <- cor(fc_matrix, use = "pairwise.complete.obs", method = "pearson")
    out_file <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9]+", "-", disease), "_corr_heatmap.png"))
    plot_heatmap(cor_mat, disease, out_file)

    message("Correlation heatmap saved for: ", disease)
  }
  message("Total unique disease: ", length(list.files(out_dir)))
}

correlate_diseases(experiment = "invivo") # Total unique disease: 113
correlate_diseases(experiment = "invitro") # Total unique disease: 122
correlate_diseases(experiment = "exvivo") # Total unique disease: 18

# Disease frequency mismatch: Only 31-38% of disease has more than 2 dataset to correlate with
high_quality %>% group_by(exp_type) %>% summarise(n_unique_diseases = n_distinct(Disease))
sapply(c("exvivo", "invitro", "invivo"), function(x)
  length(list.files(file.path("output_files/cor_results_matched/", x))))
#####