setwd('Desktop/Research/binchenlab/disease_signature/')

##### Required libraries #####
library(dplyr)
library(readxl)
library(tibble)
library(purrr)
#####

##### Check disease frequency which has unmatched controls #####
load('data/metadata_refined.RData')
get_unmatched_disease_info <- function(metadata, label, min_n = 2) {
  controls <- metadata %>% filter(Sample_type %in% c("Healthy", "Normal", "AdjacentNormal")) %>% group_by(GSE_ID, OrganRegion) %>% summarise(n_controls = n(), .groups = "drop")
  cases <- metadata %>% filter(Sample_type == "Case") %>% group_by(GSE_ID, Disease, OrganRegion) %>%
    summarise(n_cases = n(), cases_gsm = paste(GSM_ID, collapse = ","), .groups = "drop") %>% filter(n_cases >= min_n)
  matched <- cases %>% inner_join(controls, by = c("GSE_ID", "OrganRegion")) %>% filter(n_controls >= min_n)
  unmatched <- anti_join(cases, matched, by = c("Disease", "GSE_ID", "OrganRegion")) %>% select(Disease, GSE_ID, OrganRegion, n_cases, cases_gsm) %>% arrange(Disease, GSE_ID, OrganRegion)
  
  out_file <- paste0("data/unmatched_disease_info_", label, ".csv"); write.csv(unmatched, out_file, row.names = FALSE)
  
  total_disease <- n_distinct(cases$Disease); matched_disease  <- n_distinct(matched$Disease); unmatched_disease <- n_distinct(unmatched$Disease)
  message(sprintf("Saved: %s | Total: %d | Matched: %d | Unmatched: %d", out_file, total_disease, matched_disease, unmatched_disease))
  return(invisible(unmatched))
}
get_unmatched_disease_info(metadata_invitro, "invitro"); get_unmatched_disease_info(metadata_invivo, "invivo"); get_unmatched_disease_info(metadata_exvivo, "exvivo")
#####

##### Extract CASES raw reads counts from ARCHS (Unmatched only) #####
library(rhdf5)
library(dplyr)

sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
# function to extract cases expression
get_cases_expression_from_archs <- function(unmatched_df, label, expression_file_path, gsm_in_file, gene_id) {
  output_path <- file.path("data/exp_data_unmatched", label)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path("data/exp_data_unmatched", paste0("exp_extraction_log_", label, ".txt"))
  log_msg <- function(msg, log_file) { message(msg); cat(msg, "\n", file = log_file, append = TRUE) }
  if (file.exists(log_file)) file.remove(log_file)
  
  h5_file <- H5Fopen(expression_file_path)
  for (i in seq_len(nrow(unmatched_df))) {
    disease <- unmatched_df$Disease[i]
    gse_id  <- unmatched_df$GSE_ID[i]
    organ   <- unmatched_df$OrganRegion[i]
    
    filename <- paste0(sanitize_disease(disease), "_", sanitize_simple(organ), "_", sanitize_simple(gse_id), ".rds")
    output_file <- file.path(output_path, filename)
    if (file.exists(output_file)) {
      log_msg(sprintf("Skipping existing file: %s", filename), log_file)
      next
    }
    
    # collect GSMs
    case_ids <- unlist(strsplit(unmatched_df$cases_gsm[i], ","))
    case_idx <- match(case_ids, gsm_in_file)
    valid_case_idx <- which(!is.na(case_idx))
    
    if (length(valid_case_idx) < 2) {
      log_msg(sprintf("No valid CASEs GSMs for: Disease=\"%s\" / Organ=%s / GSE_ID=%s", disease, organ, gse_id), log_file)
      next
    }
    
    all_idx <- case_idx[valid_case_idx]
    all_ids <- case_ids[valid_case_idx]
    
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
    meta <- data.frame(GSM_ID = all_ids, Group = "Case", Disease = disease, OrganRegion = organ, GSE_ID = gse_id, stringsAsFactors = FALSE)
    
    saveRDS(list(expr = expr, metadata = meta), output_file)
    log_msg(sprintf("Saved CASES expression: %s", output_file), log_file)
  }
  H5Fclose(h5_file)
}

expression_file_path <- "data/human_gene_v2.5.h5"
gsm_in_file <- h5read(expression_file_path, "/meta/samples/geo_accession")
gene_id <- h5read(expression_file_path, "/meta/genes/ensembl_gene")

umdm_invitro <- read.csv("data/unmatched_disease_info_invitro.csv")
length(unique(umdm_invitro$Disease))
get_cases_expression_from_archs(umdm_invitro, "invitro", expression_file_path, gsm_in_file, gene_id)

umdm_invivo <- read.csv("data/unmatched_disease_info_invivo.csv")
length(unique(umdm_invivo$Disease))
get_cases_expression_from_archs(umdm_invivo, "invivo", expression_file_path, gsm_in_file, gene_id)

umdm_exvivo <- read.csv("data/unmatched_disease_info_exvivo.csv")
length(unique(umdm_exvivo$Disease))
get_cases_expression_from_archs(umdm_exvivo, "exvivo", expression_file_path, gsm_in_file, gene_id)
#####

##### Run prediction #####
setwd('Desktop/Research/binchenlab/gpt_geo/')
source('test_autoencoder/0_prediction.R')

exp_types <- c("invitro", "invivo", "exvivo")
base_dir  <- "../disease_signature/data/exp_data_unmatched"

ae_derived_controls <- list()
for (exp in exp_types) {
  rds_files <- list.files(file.path(base_dir, exp), pattern = "\\.rds$", full.names = TRUE)
  ae_derived_controls[[exp]] <- list()

  for (file_path in rds_files) {
    message("Processing [", exp, "]: ", file_path)
    result <- tryCatch({compute_gtex_controls(file_path, method = "spearman correlation")}, error = function(e) {message("Failed to process ", file_path, ": ", e$message)
      NULL
      })
    if (!is.null(result)) {fname <- basename(file_path); ae_derived_controls[[exp]][[fname]] <- result}
  }
}
# saveRDS(ae_derived_controls, file = "../disease_signature/data/ae_derived_controls_all.rds")
ae_derived_controls <- readRDS('data/ae_derived_controls_all.rds')
filter_ae_controls_by_tissue <- function(ae_controls, min_n = 3) {
  lapply(ae_controls, function(exp_type_list) {
    filtered <- lapply(exp_type_list, function(ctrl) {
      if (length(ctrl$control_tissue_freq) == 0) return(NULL)
      dom_tissue <- names(which.max(ctrl$control_tissue_freq)); dom_count  <- max(ctrl$control_tissue_freq)
      if (dom_count < min_n) return(NULL)
      keep_ids <- head(ctrl$top_gtex_ids, dom_count)
      list(top_gtex_ids = keep_ids, control_tissue = dom_tissue, control_tissue_freq = dom_count)
    })
    filtered[!sapply(filtered, is.null)]
  })
}
ae_derived_controls_filtered <- filter_ae_controls_by_tissue(ae_derived_controls, min_n = 10)
# saveRDS(ae_derived_controls_filtered, file = "../disease_signature/data/ae_derived_controls_filtered.rds")
#####

##### Get expression for AE derived GTEx controls #####
load('../gpt_geo/data/GTEX_exp_counts_from_octad.RData')
gtex_exp_counts <- t(gtex_exp_counts)

ae_derived_controls <- readRDS('data/ae_derived_controls_filtered.rds')

sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
get_expression_from_archs_and_gtex <- function(file_name, exp_type, ae_controls, gtex_exp_counts, out_dir = "data/exp_data_unmatched_gtex") {
  dir.create(file.path(out_dir, exp_type), recursive = TRUE, showWarnings = FALSE)
  
  # case data
  case_path <- file.path("data/exp_data_unmatched", exp_type, file_name)
  case_data <- readRDS(case_path); case_expr <- case_data$expr; case_meta <- case_data$metadata
  
  # GTEx sample IDs from ae_derived_controls
  ctrl_info <- ae_controls[[exp_type]][[file_name]]
  if (is.null(ctrl_info)) {message("Skipping ", file_name, " (no valid GTEx controls)")
    return(NULL)
  }
  ctrl_ids <- ctrl_info$top_gtex_ids; ctrl_expr <- gtex_exp_counts[ctrl_ids, , drop = FALSE]
  
  common_genes <- intersect(colnames(case_expr), colnames(ctrl_expr))
  case_expr <- case_expr[, common_genes, drop = FALSE]; ctrl_expr <- ctrl_expr[, common_genes, drop = FALSE]
  combined_expr <- rbind(case_expr, ctrl_expr)
  
  # Metadata
  ctrl_meta <- data.frame(GSM_ID = rownames(ctrl_expr), Group = "Control", Disease = case_meta$Disease[1], OrganRegion = ctrl_info$control_tissue, 
                          GSE_ID = paste0("GTEx_", case_meta$GSE_ID[1]), stringsAsFactors = FALSE)
  combined_meta <- rbind(case_meta, ctrl_meta)
  
  # Save
  out_file <- file.path(out_dir, exp_type, paste0(sanitize_disease(case_meta$Disease[1]), "_", sanitize_simple(case_meta$OrganRegion[1]), 
                                                  "_GTEx_", sanitize_simple(case_meta$GSE_ID[1]), ".rds"))
  saveRDS(list(expr = combined_expr, metadata = combined_meta), out_file)
  message("Saved combined case+GTEx: ", out_file, " | Samples: ", nrow(combined_expr), " | Genes: ", ncol(combined_expr))
  
  return(out_file)
}

# get_expression_from_archs_and_gtex("Acromegaly_Adipose_tissue_GSE57803.rds", "exvivo", ae_derived_controls, gtex_exp_counts)

# Run for all disease
# wrapper for all diseases across all exp types
combine_case_gtex_for_all_disease <- function(ae_controls, gtex_exp_counts, exp_types = c("invitro", "invivo", "exvivo"), out_dir = "data/exp_data_unmatched_gtex") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "log_combine_case_gtex_for_all_disease.txt")
  if (file.exists(log_file)) file.remove(log_file)
  results <- list()
  for (etype in exp_types) {
    files <- names(ae_controls[[etype]]); results[[etype]] <- list()
    
    for (f in files) {
      msg <- tryCatch({
        output <- capture.output({out <- get_expression_from_archs_and_gtex(f, etype, ae_controls, gtex_exp_counts, out_dir)
        results[[etype]][[f]] <- out
        }, type = "message")
        output
      }, error = function(e) {
        sprintf("ERROR: %s %s | %s", etype, f, e$message)
      })
      if (length(msg) > 0) {
        cat(msg, "\n", file = log_file, append = TRUE)
      }
    }
  }
  message("Log written to: ", log_file)
}
combine_case_gtex_for_all_disease(ae_controls = ae_derived_controls, gtex_exp_counts = gtex_exp_counts)
#####

##### Differential Gene Expression Analysis with batch correction #####
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
  
  # edgeR setup
  all_counts <- rbind(case_counts, control_counts)
  sample_type <- factor(c(rep("Case", nrow(case_counts)), rep("Control", nrow(control_counts))), levels = c("Control", "Case"))
  sample_labels <- data.frame(sample_type = sample_type)
  rownames(sample_labels) <- rownames(all_counts)
  
  keep_genes <- remLowExpr(t(all_counts), sample_labels)
  filtered_counts <- all_counts[, keep_genes, drop = FALSE]
  
  set <- newSeqExpressionSet(as.matrix(t(filtered_counts)), phenoData = AnnotatedDataFrame(sample_labels))
  
  if (use_ruv) {
    ## Step 1: design without RUV
    design <- model.matrix(~sample_type, data = pData(set))
    y <- DGEList(counts = counts(set), group = set$sample_type)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)
    top <- topTags(lrt, n = nrow(set))$table
    
    ## Step 2: pick empirical genes (not strongly DE)
    top_genes <- rownames(top)[1:min(n_topGenes, nrow(top))]
    empirical <- setdiff(rownames(set), top_genes)
    if (length(empirical) < 50) {
      message("Skipping ", rds_file, " (too few empirical control genes for RUV)")
      return(NULL)
    }
    
    ## Step 3: run RUVg
    set1 <- RUVg(set, empirical, k = k)
    pheno <- pData(set1)
    ruv_terms <- paste0("W_", seq_len(k))
    formula_str <- paste("~sample_type +", paste(ruv_terms, collapse = " + "))
    design <- model.matrix(as.formula(formula_str), data = pheno)
    dge <- DGEList(counts = counts(set1), group = set1$sample_type)
    
  } else {
    design <- model.matrix(~sample_type, data = pData(set))
    dge <- DGEList(counts = counts(set), group = set$sample_type)
  }
  
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
  message("Saved DE results: ", out_file)
  return(out_file)
}

# Function: run DE for all RDS files in a directory
do_DGE_for_directory <- function(input_dir, output_dir, use_ruv = FALSE, k = 1) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(output_dir, "00_dge_log.txt")
  if (file.exists(log_file)) file.remove(log_file)
  
  rds_files <- list.files(input_dir, full.names = TRUE, pattern = "\\.rds$")
  for (f in rds_files) {
    tryCatch({
      msg <- capture.output(do_DGE(f, out_dir = output_dir, use_ruv = use_ruv, k = k), type = "message")
      cat(msg, "\n", file = log_file, append = TRUE)
    }, error = function(e) {
      msg <- sprintf("ERROR: %s | %s", f, e$message)
      message(msg)
      cat(msg, "\n", file = log_file, append = TRUE)
    })
  }
  message(sprintf("Log written to: %s", log_file))
}

# do_DGE_for_directory("data/exp_data_unmatched_gtex/invitro/", "output_files/dge_results_unmatched/invitro", use_ruv = TRUE, k = 1)
# do_DGE_for_directory("data/exp_data_unmatched_gtex/invivo/",  "output_files/dge_results_unmatched/invivo",  use_ruv = TRUE, k = 1)
# do_DGE_for_directory("data/exp_data_unmatched_gtex/exvivo/",  "output_files/dge_results_unmatched/exvivo",  use_ruv = TRUE, k = 1)
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
  dge_files <- list.files(file.path("output_files/dge_results_unmatched/", exp), full.names = TRUE, pattern = "\\.csv$")
  qc_all[[exp]] <- bind_rows(lapply(dge_files, summarize_dge, exp_type = exp))
}
qc_summary <- bind_rows(qc_all) %>% mutate(Disease = sub("_.*", "", file))
head(qc_summary); write.csv(qc_summary, "output_files/dge_results_unmatched/00_qc_summary.csv", row.names = FALSE)

# hist(qc_summary$sig_genes)
# hist(qc_summary$prop_pval)
# hist(qc_summary$mean_abs_logFC)
# hist(qc_summary$snr)

high_quality <- qc_summary %>% filter(sig_genes >= 10 & sig_genes <= 10000, prop_pval > 0.01 & prop_pval < 0.5, 
                                      mean_abs_logFC > 0.1 & mean_abs_logFC < 3, snr > 0.8 & snr < 5, up > 0 & down > 0)
head(high_quality); write.csv(high_quality, "output_files/dge_results_unmatched/01_high_quality_DE.csv", row.names = FALSE)

low_quality <- anti_join(qc_summary, high_quality)
head(low_quality); write.csv(low_quality, "output_files/dge_results_unmatched/02_low_quality_DE.csv", row.names = FALSE)
#####

##### Correlate high quality signature #####
library(dplyr)
library(purrr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(ggdendro)
library(patchwork)

high_quality <- read.csv("output_files/dge_results_unmatched/01_high_quality_DE.csv")

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
correlate_diseases <- function(experiment, de_res_dir = file.path("output_files/dge_results_unmatched/", experiment), 
                               out_dir = file.path("output_files/cor_results_unmatched", experiment)) {
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

correlate_diseases(experiment = "invivo") # Total unique disease: 41
correlate_diseases(experiment = "invitro") # Total unique disease: 27
correlate_diseases(experiment = "exvivo") # Total unique disease: 1

# Disease frequency mismatch: Only 1-33% of disease has more than 2 dataset to correlate with
high_quality %>% group_by(exp_type) %>% summarise(n_unique_diseases = n_distinct(Disease))
sapply(c("exvivo", "invitro", "invivo"), function(x)
  length(list.files(file.path("output_files/cor_results_unmatched/", x))))
#####