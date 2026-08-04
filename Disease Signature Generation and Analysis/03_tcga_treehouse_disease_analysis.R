setwd('Desktop/Research/binchenlab/disease_signature/')

##### Required libraries #####
library(dplyr)
library(rhdf5)
library(octad)
library(tidyr)
library(stringr)
library(tibble)
library(readxl)
#####

######################### TCGA #########################
##### TCGA #####
# Counts
octad_data_path <- 'data/octad.counts.and.tpm.h5'
# h5ls(octad_data_path)
samples_id <- h5read(octad_data_path, '/meta/samples')
tcga_sample_index <-  grep('TCGA-', samples_id)
tcga_sample_id <- samples_id[tcga_sample_index]
transcripts_id <- h5read(octad_data_path, '/meta/transcripts')

tcga_log2_counts <- h5read(octad_data_path, '/data/count', index = list(NULL, tcga_sample_index))
colnames(tcga_log2_counts) <- tcga_sample_id
rownames(tcga_log2_counts) <- transcripts_id

tcga_counts <- (2 ^ tcga_log2_counts - 1)
head(tcga_counts)[, 1:5]
tcga_counts <- round(tcga_counts)
tcga_counts <- t(tcga_counts)
dim(tcga_counts)
# saveRDS(tcga_counts, file = 'data/tcga/tcga_counts.rds')

# Metadata
phenoDF <- get_ExperimentHub_data('EH7274')
tcga_metadata <- phenoDF %>% filter(grepl("TCGA", sample.id, ignore.case = TRUE))
tcga_metadata <- tcga_metadata %>% select(-metastatic_site, -tumor_grade, -tumor_stage, -subtype, -age_in_year, -gender,
                                          -gain_list, -loss_list, -mutation_list, -data.source, -read.count.file)
table(tcga_metadata$sample.type)
tcga_metadata <- tcga_metadata %>% filter(sample.type %in% c('adjacent', 'primary'))

# disease_freq <- tcga_metadata %>% group_by(cancer, biopsy.site, sample.type) %>% tally() %>% spread(sample.type, n, fill = 0) %>% arrange(desc(primary))
# write.csv(disease_freq, file = 'data/tcga/disease_freq.csv')

# matched_disease_metadata <- tcga_metadata %>% filter(sample.type %in% c("primary", "adjacent")) %>% group_by(cancer, biopsy.site) %>%
#   summarise(n_cases    = sum(sample.type == "primary"), n_controls = sum(sample.type == "adjacent"), 
#             cases_id   = paste(sample.id[sample.type == "primary"], collapse = ","),
#             controls_id= paste(sample.id[sample.type == "adjacent"], collapse = ","), .groups = "drop") %>%
#   filter(n_cases >= 2 & n_controls >= 2) %>% mutate(Disease = str_to_title(cancer), OrganRegion  = str_to_title(biopsy.site), Control_type = "adjacent") %>%
#   select(Disease, OrganRegion, Control_type, n_cases, n_controls, cases_id, controls_id)
# unmatched_disease_metadata <- tcga_metadata %>% filter(sample.type == "primary") %>% group_by(cancer, biopsy.site) %>%
#   summarise(n_cases = n(), cases_id = paste(sample.id, collapse = ","), .groups = "drop") %>%
#   filter(!(cancer %in% matched_disease_metadata$Disease)) %>% mutate(Disease = str_to_title(cancer), OrganRegion = str_to_title(biopsy.site)) %>% 
#   select(Disease, OrganRegion, n_cases, cases_id)

matched_disease_metadata <- tcga_metadata %>% filter(sample.type %in% c("primary", "adjacent")) %>% group_by(cancer, biopsy.site) %>%
  summarise(n_cases = sum(sample.type == "primary"), n_controls = sum(sample.type == "adjacent"), cases_id = paste(sample.id[sample.type == "primary"], collapse = ","),
            controls_id = paste(sample.id[sample.type == "adjacent"], collapse = ","), .groups = "drop") %>% filter(n_cases >= 2 & n_controls >= 2) %>%
  mutate(Disease = str_to_title(cancer), OrganRegion  = str_to_title(biopsy.site), Control_type = "adjacent") %>%
  select(Disease, OrganRegion, Control_type, n_cases, n_controls, cases_id, controls_id)

unmatched_disease_metadata <- tcga_metadata %>% filter(sample.type == "primary") %>% group_by(cancer, biopsy.site) %>%
  summarise(n_cases = n(), cases_id = paste(sample.id, collapse = ","), .groups = "drop") %>%
  filter(!(tolower(cancer) %in% tolower(unique(matched_disease_metadata$Disease)))) %>% 
  mutate(Disease = str_to_title(cancer), OrganRegion = str_to_title(biopsy.site)) %>% select(Disease, OrganRegion, n_cases, cases_id)
intersect(matched_disease_metadata$Disease, unmatched_disease_metadata$Disease)

# write.csv(matched_disease_metadata, paste0("data/tcga/matched_disease_info.csv"), row.names = FALSE)
# write.csv(unmatched_disease_metadata, paste0("data/tcga/unmatched_disease_info.csv"), row.names = FALSE)
#####

##### Extract raw reads counts from OCTAD #####
sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

get_expression_from_tcga <- function(mdm, tcga_counts, output_dir) {
  output_path <- file.path(output_dir)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_len(nrow(mdm))) {
    # disease <- mdm$Disease[i]; organ <- mdm$OrganRegion[i]; control_type <- mdm$Control_type[i]
    disease <- mdm$Disease_Post[i]; organ <- mdm$OrganRegion_Post[i]; control_type <- mdm$Control_type[i]
    filename <- paste0(sanitize_disease(disease), "_", sanitize_simple(organ), "_", sanitize_simple(control_type), "_TCGA.rds")
    output_file <- file.path(output_path, filename)
    
    if (file.exists(output_file)) {message("Skipping existing file: ", filename)
      next
    }
    
    case_ids <- unlist(strsplit(mdm$cases_id[i], ","))
    control_ids <- unlist(strsplit(mdm$controls_id[i], ","))
    all_ids <- c(case_ids, control_ids)
    expr <- tcga_counts[rownames(tcga_counts) %in% all_ids, , drop = FALSE]
    
    if (nrow(expr) == 0) {message("No matching samples for: ", disease)
      next
    }
    
    # metadata
    meta <- data.frame(Sample_ID = rownames(expr), Group = ifelse(rownames(expr) %in% case_ids, "Case", "Control"), Disease = disease, OrganRegion = organ,
                       Control_type = control_type, Source = "TCGA", stringsAsFactors = FALSE)
    saveRDS(list(expr = expr, metadata = meta), output_file)
    message("Saved: ", output_file)
  }
}

tcga_counts <- readRDS("data/tcga/tcga_counts.rds")
# matched_disease_metadata <- read.csv("data/tcga/matched_disease_info.csv")
matched_disease_metadata <- read_excel("tcga_meta_Mapped.xlsx")
get_expression_from_tcga(matched_disease_metadata, tcga_counts, output_dir = 'data/tcga/exp_data/matched')
#####

##### Differential Gene Expression Analysis #####
library(edgeR)
library(RUVSeq)

# Function: remove lowly expressed genes
remLowExpr <- function(counts, sample_type) {
  x <- DGEList(counts = round(counts), group = sample_type$sample_type)
  cpm_x <- cpm(x)
  keep.exprs <- rowSums(cpm_x > 1) >= min(table(sample_type$sample_type))  # at least expressed in min group size
  return(keep.exprs)
}
# Function: run DE and save DE results
do_DGE <- function(rds_file, out_dir = "results", use_ruv = TRUE, k = 1, n_topGenes = 2000) {
  data <- readRDS(rds_file)
  expr <- data$expr
  meta <- data$metadata
  
  # align meta with expr rows
  expr <- expr[meta$Sample_ID, , drop = FALSE]
  
  case_counts    <- expr[meta$Group == "Case", , drop = FALSE]
  control_counts <- expr[meta$Group == "Control", , drop = FALSE]
  
  if (nrow(case_counts) == 0 | nrow(control_counts) == 0) {
    message("Skipping ", rds_file, " (missing case or control samples)")
    return(NULL)
  }
  
  sample_type <- factor(meta$Group, levels = c("Control", "Case"))
  sample_labels <- data.frame(sample_type = sample_type)
  rownames(sample_labels) <- meta$Sample_ID
  keep_genes <- remLowExpr(t(expr), sample_labels)
  filtered_counts <- expr[, keep_genes, drop = FALSE]
  set <- newSeqExpressionSet(as.matrix(t(filtered_counts)), phenoData = AnnotatedDataFrame(sample_labels))
  design <- model.matrix(~ sample_type, data = pData(set))
  
  # RUVSeq correction
  if (use_ruv) {
    set <- betweenLaneNormalization(set, which = "upper")
    mat <- counts(set)
    # empirical control genes = least variable genes
    gene_var <- apply(mat, 1, var)
    gene_var <- gene_var[is.finite(gene_var) & gene_var > 0]
    empirical <- head(names(sort(gene_var)), n = min(n_topGenes, length(gene_var)))
    raw_counts <- counts(set)
    set <- RUVg(newSeqExpressionSet(as.matrix(raw_counts), phenoData = pData(set)), cIdx = empirical, k = k, isLog = FALSE)
    design <- model.matrix(~ sample_type + W_1, data = pData(set))
  }
  dge <- DGEList(counts = counts(set), group = set$sample_type)
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = 2)
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res$identifier <- rownames(res)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(rds_file)), ".csv"))
  write.csv(res, out_file, row.names = FALSE)
  return(sprintf("SUCCESS: Saved DE results in: %s", out_file))
}

# # Example: Execute DE for one file
# do_DGE(rds_file = 'data/tcga/exp_data/Lung-Adenocarcinoma_Lung_adjacent_TCGA.rds', out_dir = 'output_files/dge_results_tcga/matched/')
# Execute: DE for all files in a directory
do_DGE_for_directory <- function(input_dir = "data/tcga/exp_data", output_dir = "output_files/dge_results_tcga/matched") {
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  lapply(rds_files, function(f) {
    message("Processing: ", f)
    do_DGE(rds_file = f, out_dir = output_dir)
  })
}
do_DGE_for_directory(input_dir = "data/tcga/exp_data/matched", output_dir = "output_files/dge_results_tcga/matched")
#####

##### Extract CASES raw reads counts from OCTAD (Unmatched only) #####
sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

get_expression_from_tcga_unmatched <- function(udm, tcga_counts, output_dir) {
  output_path <- file.path(output_dir)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  for (i in seq_len(nrow(udm))) {
    disease <- udm$Disease_mapped[i]; organ <- udm$OrganRegion[i]
    filename <- paste0(sanitize_disease(disease), "_", sanitize_simple(organ), "_TCGA.rds")
    output_file <- file.path(output_path, filename)
    
    if (file.exists(output_file)) {message("Skipping existing file: ", filename)
      next
    }
    case_ids <- unlist(strsplit(udm$cases_id[i], ","))
    expr <- tcga_counts[rownames(tcga_counts) %in% case_ids, , drop = FALSE]
    if (nrow(expr) == 0) {message("No matching samples for: ", disease)
      next
    }
    
    # Metadata
    meta <- data.frame(Sample_ID = rownames(expr), Group = ifelse(rownames(expr) %in% case_ids, "Case", "Control"), 
                       Disease = disease, OrganRegion = organ, Source = "TCGA", stringsAsFactors = FALSE)
    saveRDS(list(expr = expr, metadata = meta), output_file)
    message("Saved: ", output_file)
  }
}

tcga_counts <- readRDS("data/tcga/tcga_counts.rds")
# unmatched_disease_metadata <- read.csv("data/tcga/unmatched_disease_info.csv")
unmatched_disease_metadata <- read_excel("unmatched_disease_info_mapped.xlsx")
get_expression_from_tcga_unmatched(unmatched_disease_metadata, tcga_counts, output_dir = "data/tcga/exp_data/unmatched_cases_only_new/")
#####

##### Run prediction #####
setwd('Desktop/Research/binchenlab/gpt_geo/')

source('test_autoencoder/0_prediction.R')
base_dir <- "../disease_signature/data/tcga/exp_data/unmatched_cases_only"
ae_derived_controls <- list()
rds_files <- list.files(base_dir, pattern = "\\.rds$", full.names = TRUE)
for (file_path in rds_files) {
  message("Processing: ", file_path)
  result <- tryCatch({compute_gtex_controls(file_path, method = "spearman correlation")}, error = function(e) {message("Failed to process ", file_path, ": ", e$message)
    NULL
  })
  if (!is.null(result)) {fname <- basename(file_path); ae_derived_controls[[fname]] <- result}
}
# saveRDS(ae_derived_controls, file = "../disease_signature/data/tcga/ae_derived_controls_all.rds")

ae_derived_controls <- readRDS('data/tcga/ae_derived_controls_all.rds')
filter_ae_controls_by_tissue <- function(ae_controls, min_n = 3) {
  filtered <- lapply(ae_controls, function(ctrl) {
    if (is.null(ctrl$control_tissue_freq) || length(ctrl$control_tissue_freq) == 0) return(NULL)
    dom_tissue <- names(ctrl$control_tissue_freq)[1]; dom_count  <- ctrl$control_tissue_freq[1]
    if (dom_count < min_n) return(NULL)
    keep_ids <- head(ctrl$top_gtex_ids, dom_count)
    list(top_gtex_ids = keep_ids, control_tissue = dom_tissue, control_tissue_freq = dom_count)
  })
  filtered[!sapply(filtered, is.null)]
}
ae_derived_controls_filtered <- filter_ae_controls_by_tissue(ae_derived_controls, min_n = 10)
# saveRDS(ae_derived_controls_filtered, file = "../disease_signature/data/tcga/ae_derived_controls_filtered.rds")
#####

##### Get expression for AE derived GTEx controls #####
load('../gpt_geo/data/GTEX_exp_counts_from_octad.RData')
gtex_exp_counts <- t(gtex_exp_counts)
ae_derived_controls <- readRDS('data/tcga/ae_derived_controls_filtered.rds')

sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
get_expression_from_octad_and_gtex <- function(file_name, ae_controls, gtex_exp_counts, out_dir = "data/tcga/exp_data/unmatched") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  case_path <- file.path("data/tcga/exp_data/unmatched_cases_only", file_name)
  case_data <- readRDS(case_path)
  case_expr <- as.matrix(case_data$expr)
  case_meta <- case_data$metadata
  ctrl_info <- ae_controls[[file_name]]
  if (is.null(ctrl_info)) {message("Skipping ", file_name, " (no valid GTEx controls)")
    return(NULL)
  }
  ctrl_ids  <- ctrl_info$top_gtex_ids; ctrl_expr <- as.matrix(gtex_exp_counts[ctrl_ids, , drop = FALSE])
  common_genes <- intersect(colnames(case_expr), colnames(ctrl_expr))
  case_expr <- case_expr[, common_genes, drop = FALSE]; ctrl_expr <- ctrl_expr[, common_genes, drop = FALSE]
  ctrl_expr <- ctrl_expr[, match(colnames(case_expr), colnames(ctrl_expr)), drop = FALSE]
  stopifnot(identical(colnames(case_expr), colnames(ctrl_expr)))
  rownames(case_expr) <- make.unique(rownames(case_expr))
  rownames(ctrl_expr) <- make.unique(rownames(ctrl_expr))
  combined_expr <- rbind(case_expr, ctrl_expr)
  
  # Meta
  ctrl_meta <- data.frame(Sample_ID = rownames(ctrl_expr), Group = "Control", Disease = case_meta$Disease[1], OrganRegion = ctrl_info$control_tissue,
                          GSE_ID = paste0("GTEx_", case_meta$GSE_ID[1]), Source = "GTEx", stringsAsFactors = FALSE)
  ctrl_meta <- ctrl_meta[, colnames(case_meta), drop = FALSE]
  combined_meta <- rbind(case_meta, ctrl_meta)
  
  out_file <- file.path(out_dir, paste0(sanitize_disease(case_meta$Disease[1]), "_", sanitize_simple(case_meta$OrganRegion[1]), "_GTEx", ".rds"))
  saveRDS(list(expr = combined_expr, metadata = combined_meta), out_file)
  message("Saved combined case+GTEx: ", out_file, " | Samples: ", nrow(combined_expr), " | Genes: ", ncol(combined_expr))
  return(out_file)
}
combine_case_gtex_for_all_disease <- function(ae_controls, gtex_exp_counts, out_dir = "data/tcga/exp_data/unmatched") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "log_combine_case_gtex_for_all_disease.txt")
  if (file.exists(log_file)) file.remove(log_file)
  
  results <- list()
  files <- names(ae_controls)
  for (f in files) {
    msg <- tryCatch({
      output <- capture.output({out <- get_expression_from_octad_and_gtex(f, ae_controls, gtex_exp_counts, out_dir); results[[f]] <- out}, type = "message")
      output
    }, error = function(e) {
      sprintf("ERROR: %s | %s", f, e$message)
    })
    if (length(msg) > 0) {cat(msg, "\n", file = log_file, append = TRUE)}
  }
  message("Log written to: ", log_file); invisible(results)
}
combine_case_gtex_for_all_disease(ae_controls = ae_derived_controls, gtex_exp_counts = gtex_exp_counts, out_dir = "data/tcga/exp_data/unmatched_new/")
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
  case_counts    <- expr[meta$Group == "Case", , drop = FALSE]
  control_counts <- expr[meta$Group == "Control", , drop = FALSE]
  
  if (nrow(case_counts) == 0 | nrow(control_counts) == 0) {
    message("Skipping ", rds_file, " (missing case or control samples)")
    return(NULL)
  }
  all_counts <- rbind(case_counts, control_counts)
  sample_type <- factor(meta$Group, levels = c("Control", "Case"))
  sample_labels <- data.frame(sample_type = sample_type)
  rownames(sample_labels) <- rownames(all_counts)
  keep_genes <- remLowExpr(t(all_counts), sample_labels)
  filtered_counts <- all_counts[, keep_genes, drop = FALSE]
  set <- newSeqExpressionSet(as.matrix(t(filtered_counts)), phenoData = AnnotatedDataFrame(sample_labels))
  
  if (use_ruv) {
    design <- model.matrix(~ sample_type, data = pData(set))
    y <- DGEList(counts = counts(set), group = set$sample_type)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)
    top <- topTags(lrt, n = nrow(set))$table
    top_genes <- rownames(top)[1:min(n_topGenes, nrow(top))]
    empirical <- setdiff(rownames(set), top_genes)
    if (length(empirical) < 50) {message("Skipping ", rds_file, " (too few empirical control genes for RUV)")
      return(NULL)
    }
    set1 <- RUVg(set, empirical, k = k)
    pheno <- pData(set1)
    ruv_terms <- paste0("W_", seq_len(k))
    formula_str <- paste("~ sample_type +", paste(ruv_terms, collapse = " + "))
    design <- model.matrix(as.formula(formula_str), data = pheno)
    dge <- DGEList(counts = counts(set1), group = set1$sample_type)
    
  } else {
    design <- model.matrix(~ sample_type, data = pData(set))
    dge <- DGEList(counts = counts(set), group = set$sample_type)
  }
  
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = 2)
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res$identifier <- rownames(res)
  
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
      msg <- capture.output(do_DGE(f, out_dir = output_dir, use_ruv = use_ruv, k = k), type = "message"); cat(msg, "\n", file = log_file, append = TRUE)
    }, error = function(e) {msg <- sprintf("ERROR: %s | %s", f, e$message); message(msg); cat(msg, "\n", file = log_file, append = TRUE)})
  }
  message(sprintf("Log written to: %s", log_file))
}
do_DGE_for_directory(input_dir = "data/tcga/exp_data/unmatched", output_dir = "output_files/dge_results_tcga/unmatched", use_ruv = TRUE, k = 1)
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

exp_dirs <- c("matched", "unmatched")
qc_all <- list()
for (exp in exp_dirs) {
  dge_files <- list.files(file.path("output_files/dge_results_tcga/", exp), full.names = TRUE, pattern = "\\.csv$")
  qc_all[[exp]] <- bind_rows(lapply(dge_files, summarize_dge, exp_type = exp))
}
qc_summary <- bind_rows(qc_all) %>% mutate(Disease = sub("_.*", "", file))
head(qc_summary); write.csv(qc_summary, "output_files/dge_results_tcga/00_qc_summary.csv", row.names = FALSE)

# qc_summary <- read.csv("output_files/dge_results_tcga/00_qc_summary.csv")
# hist(qc_summary$sig_genes)
# hist(qc_summary$prop_pval)
# hist(qc_summary$mean_abs_logFC)
# hist(qc_summary$snr)

high_quality <- qc_summary %>% filter(sig_genes >= 10 & sig_genes <= 10000, prop_pval > 0.01 & prop_pval < 0.5, 
                                      mean_abs_logFC > 0.1 & mean_abs_logFC < 3, snr > 0.8 & snr < 5, up > 0 & down > 0)
head(high_quality); write.csv(high_quality, "output_files/dge_results_tcga/01_high_quality_DE.csv", row.names = FALSE)

low_quality <- anti_join(qc_summary, high_quality)
head(low_quality); write.csv(low_quality, "output_files/dge_results_tcga/02_low_quality_DE.csv", row.names = FALSE)
#####
######################### TCGA - END #########################

######################### Treehouse #########################
##### Treehouse #####
th_counts <- read.csv('data/treehouse/GSE294351_Tumor-25.01-Polya_ensembl_counts_60498genes_2025-03-03.tsv', sep = '\t', row.names = 1)
# dim(th_counts)
# head(th_counts)[, 1:10]
th_counts <- round(th_counts)
rownames(th_counts) <- sub("\\..*$", "", rownames(th_counts))
# head(th_counts)[, 1:10]
# length(unique(rownames(th_counts)))

metadata <- read.csv('data/treehouse/GSE294351_clinical_Treehouse-Tumor-Compendium-25.01-PolyA_20250131v1.tsv', sep = '\t')

head(metadata, 2)
length(intersect(colnames(th_counts), metadata$th_dataset_id))
# dim(th_counts)
# head(th_counts)[, 1:10]
# dim(metadata)
# head(metadata)[, 1:10]
# length(unique(metadata$th_dataset_id))
# sum(is.na(metadata$th_dataset_id))
# head(setdiff(colnames(th_counts), metadata$th_dataset_id))
# head(setdiff(metadata$th_dataset_id, colnames(th_counts)))
# any(grepl("^\\s|\\s$", metadata$th_dataset_id))

metadata$sample_id_norm <- gsub("-", ".", metadata$th_dataset_id)
# length(intersect(colnames(th_counts), metadata$sample_id_norm))
metadata <- metadata[metadata$sample_id_norm %in% colnames(th_counts), ]
metadata <- metadata[match(colnames(th_counts), metadata$sample_id_norm), ]
# all.equal(colnames(th_counts), metadata$sample_id_norm)

table(metadata$source_name)
# tcga_metadata <- metadata %>% filter(source_name == 'TCGA')
# length(unique(tcga_metadata$disease))
metadata <- metadata %>% filter(!source_name == 'TCGA')
# colnames(metadata)
metadata <- metadata %>% select(-age_at_dx, -organism, -study_accession, -study_dataset_id, -study_donor_id)
abbrev_table <- tribble(~source_name, ~source_full_name, ~source,
                        "CBTN", "Children’s Brain Tumor Network", "CBTN",
                        "EGA", "European Genome-phenome Archive", "EGA",
                        "ICGC", "International Cancer Genome Consortium", "ICGC",
                        "Sequence Read Archive", "Sequence Read Archive", "SRA",
                        "St. Jude Cloud", "St. Jude Cloud", "SJC",
                        "TARGET", "Therapeutically Applicable Research to Generate Effective Treatments", "TARGET",
                        "unavailable", "unavailable", "UA")
metadata <- metadata %>% left_join(abbrev_table, by = "source_name") %>% select(everything(), source_full_name, source, -source_name)
metadata <- metadata %>% mutate(disease_post = str_split_fixed(icd_disease, ":", n = 2)[,2] %>% str_trim() %>% str_remove(", NOS$") %>% str_to_title())
# length(unique(metadata$disease_post))
# length(unique(metadata$disease))
mismatched_disease <- metadata %>% distinct(disease, disease_post) %>% arrange(disease)

disease_freq <- metadata %>% distinct(source, disease_post, disease, sample_id_norm) %>% group_by(source, disease_post, disease) %>% summarise(n_samples = n(), .groups = "drop")
# length(unique(disease_freq$disease))
# length(unique(disease_freq$disease_post))
# head(disease_freq)
# disease_freq_filtered <- disease_freq %>% filter(n_samples >= 2) %>% arrange(desc(n_samples))
# length(unique(disease_freq_filtered$disease_post))
# length(unique(disease_freq_filtered$disease))

# dim(th_counts)
th_counts_filtered <- th_counts[, colnames(th_counts) %in% metadata$sample_id_norm]
# dim(th_counts_filtered)
th_counts_filtered <- t(th_counts_filtered)
# dim(th_counts_filtered)
# head(th_counts_filtered)[, 1:5]

# saveRDS(th_counts_filtered, file = 'data/treehouse/treehouse_counts.rds')
# write.csv(mismatched_disease, file = 'data/treehouse/ICD_mismatched_disease.csv')
# write.csv(disease_freq, file = 'data/treehouse/disease_freq.csv')
# write.csv(metadata, file = 'data/treehouse/refined_metadata.csv')
#####

##### Extract CASES raw reads counts from Treehouse (all samples are Unmatched/Cases) #####
sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
get_expression_from_treehouse <- function(mdm, treehouse_counts, output_dir) {
  output_path <- file.path(output_dir)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  if (nrow(treehouse_counts) == nrow(mdm) && all(mdm$sample_id_norm %in% rownames(treehouse_counts))) {treehouse_counts <- t(treehouse_counts)}
  
  groups <- mdm %>% dplyr::mutate(use_disease_post = !(is.na(disease_mapped) | disease_mapped == "")) %>% 
    dplyr::mutate(disease_final = ifelse(use_disease_post, disease_mapped, disease_post)) %>% 
    dplyr::group_by(disease_final, source, use_disease_post) %>%
    dplyr::summarise(samples = list(sample_id_norm), .groups = "drop")
  
  for (i in seq_len(nrow(groups))) {
    disease <- groups$disease_final[i]; source  <- groups$source[i]; sample_ids <- unlist(groups$samples[i])
    disease_sanitized <- sanitize_disease(disease)
    disease_sanitized <- stringr::str_to_title(disease_sanitized)
    if (!groups$use_disease_post[i]) {disease_sanitized <- paste0(disease_sanitized, "-MM-ICD")}
    filename <- paste0(disease_sanitized, "_", sanitize_simple(source), "_Treehouse.rds")
    output_file <- file.path(output_path, filename)
    if (file.exists(output_file)) {message("Skipping existing file: ", filename)
      next
    }
    expr <- treehouse_counts[, colnames(treehouse_counts) %in% sample_ids, drop = FALSE]
    if (ncol(expr) < 2) {message("No matching samples for: ", disease)
      next
    }
    meta <- data.frame(Sample_ID = colnames(expr), Group = "Case", Disease = disease_sanitized, Source = source, stringsAsFactors = FALSE)
    expr <- t(expr)
    saveRDS(list(expr = expr, metadata = meta), output_file)
    message("Saved: ", output_file)
  }
}
treehouse_counts <- readRDS("data/treehouse/treehouse_counts.rds")
# metadata <- read.csv("data/treehouse/refined_metadata.csv")
metadata <- read_excel('treehouse_meta_with_disease_mapped.xlsx')
get_expression_from_treehouse(mdm = metadata, treehouse_counts = treehouse_counts, output_dir = 'data/treehouse/exp_data')
#####

##### Run prediction #####
setwd('Desktop/Research/binchenlab/gpt_geo/')

source('test_autoencoder/0_prediction.R')
base_dir <- "../disease_signature/data/treehouse/exp_data"
ae_derived_controls <- list()
rds_files <- list.files(base_dir, pattern = "\\.rds$", full.names = TRUE)
for (file_path in rds_files) {
  message("Processing: ", file_path)
  result <- tryCatch({compute_gtex_controls(file_path, method = "spearman correlation")}, error = function(e) {message("Failed to process ", file_path, ": ", e$message)
    NULL
  })
  if (!is.null(result)) {fname <- basename(file_path); ae_derived_controls[[fname]] <- result}
}
saveRDS(ae_derived_controls, file = "../disease_signature/data/treehouse/ae_derived_controls_all.rds")

ae_derived_controls <- readRDS('data/treehouse/ae_derived_controls_all.rds')
filter_ae_controls_by_tissue <- function(ae_controls, min_n = 3) {
  filtered <- lapply(ae_controls, function(ctrl) {
    if (is.null(ctrl$control_tissue_freq) || length(ctrl$control_tissue_freq) == 0) return(NULL)
    dom_tissue <- names(ctrl$control_tissue_freq)[1]; dom_count  <- ctrl$control_tissue_freq[1]
    if (dom_count < min_n) return(NULL)
    keep_ids <- head(ctrl$top_gtex_ids, dom_count)
    list(top_gtex_ids = keep_ids, control_tissue = dom_tissue, control_tissue_freq = dom_count)
  })
  filtered[!sapply(filtered, is.null)]
}
ae_derived_controls_filtered <- filter_ae_controls_by_tissue(ae_derived_controls, min_n = 10)
saveRDS(ae_derived_controls_filtered, file = "../disease_signature/data/treehouse/ae_derived_controls_filtered.rds")
#####

##### Get expression for AE derived GTEx controls #####
load('../gpt_geo/data/GTEX_exp_counts_from_octad.RData')
gtex_exp_counts <- t(gtex_exp_counts)
ae_derived_controls <- readRDS('data/treehouse/ae_derived_controls_filtered.rds')

sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
sanitize_simple  <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
get_expression_from_treehouse_and_gtex <- function(file_name, ae_controls, gtex_exp_counts, out_dir = "data/treehouse/exp_data_gtex") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  case_path <- file.path("data/treehouse/exp_data", file_name)
  case_data <- readRDS(case_path)
  case_expr <- as.matrix(case_data$expr)
  case_meta <- case_data$metadata
  ctrl_info <- ae_controls[[file_name]]
  if (is.null(ctrl_info)) {message("Skipping ", file_name, " (no valid GTEx controls)")
    return(NULL)
  }
  ctrl_ids  <- ctrl_info$top_gtex_ids; ctrl_expr <- as.matrix(gtex_exp_counts[ctrl_ids, , drop = FALSE])
  common_genes <- intersect(colnames(case_expr), colnames(ctrl_expr))
  case_expr <- case_expr[, common_genes, drop = FALSE]; ctrl_expr <- ctrl_expr[, common_genes, drop = FALSE]
  ctrl_expr <- ctrl_expr[, match(colnames(case_expr), colnames(ctrl_expr)), drop = FALSE]
  stopifnot(identical(colnames(case_expr), colnames(ctrl_expr)))
  rownames(case_expr) <- make.unique(rownames(case_expr))
  rownames(ctrl_expr) <- make.unique(rownames(ctrl_expr))
  combined_expr <- rbind(case_expr, ctrl_expr)
  
  # Meta
  ctrl_meta <- data.frame(Sample_ID = rownames(ctrl_expr), Group = "Control", Disease = case_meta$Disease[1], OrganRegion = ctrl_info$control_tissue,
                          GSE_ID = paste0("GTEx_", case_meta$GSE_ID[1]), Source = "GTEx", stringsAsFactors = FALSE)
  ctrl_meta <- ctrl_meta[, colnames(case_meta), drop = FALSE]
  combined_meta <- rbind(case_meta, ctrl_meta)
  
  out_file <- file.path(out_dir, paste0(sanitize_disease(case_meta$Disease[1]), "_", sanitize_simple(case_meta$Source[1]), "_GTEx", ".rds"))
  saveRDS(list(expr = combined_expr, metadata = combined_meta), out_file)
  message("Saved combined case+GTEx: ", out_file, " | Samples: ", nrow(combined_expr), " | Genes: ", ncol(combined_expr))
  return(out_file)
}
combine_case_gtex_for_all_disease <- function(ae_controls, gtex_exp_counts, out_dir = "data/treehouse/exp_data_gtex") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "log_combine_case_gtex_for_all_disease.txt")
  if (file.exists(log_file)) file.remove(log_file)
  
  results <- list()
  files <- names(ae_controls)
  for (f in files) {
    msg <- tryCatch({
      output <- capture.output({out <- get_expression_from_treehouse_and_gtex(f, ae_controls, gtex_exp_counts, out_dir); results[[f]] <- out}, type = "message")
      output
    }, error = function(e) {
      sprintf("ERROR: %s | %s", f, e$message)
    })
    if (length(msg) > 0) {cat(msg, "\n", file = log_file, append = TRUE)}
  }
  message("Log written to: ", log_file); invisible(results)
}
combine_case_gtex_for_all_disease(ae_controls = ae_derived_controls, gtex_exp_counts = gtex_exp_counts, out_dir = "data/treehouse/exp_data_gtex_new")
#####

# geo <- readRDS('data/exp_data_unmatched_gtex/invivo/Abortion-Habitual_Endometrium_GTEx_GSE65099.rds')
# tcga <- readRDS('data/tcga/exp_data/unmatched/Adrenal-Cortex-Cancer_Adrenal_Gland_GTEx.rds')
# tree <- readRDS('data/treehouse/exp_data_gtex/Acute-Megakaryoblastic-Leukemia_SJC_GTEx.rds')
# e_geo <- geo$expr; e_tcga <- tcga$expr; e_tree <- tree$expr
# m_geo <- geo$metadata; m_tcga <- tcga$metadata; m_tree <- tree$metadata

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
  case_counts    <- expr[meta$Group == "Case", , drop = FALSE]
  control_counts <- expr[meta$Group == "Control", , drop = FALSE]
  
  if (nrow(case_counts) == 0 | nrow(control_counts) == 0) {
    message("Skipping ", rds_file, " (missing case or control samples)")
    return(NULL)
  }
  all_counts <- rbind(case_counts, control_counts)
  sample_type <- factor(meta$Group, levels = c("Control", "Case"))
  sample_labels <- data.frame(sample_type = sample_type)
  rownames(sample_labels) <- rownames(all_counts)
  keep_genes <- remLowExpr(t(all_counts), sample_labels)
  filtered_counts <- all_counts[, keep_genes, drop = FALSE]
  set <- newSeqExpressionSet(as.matrix(t(filtered_counts)), phenoData = AnnotatedDataFrame(sample_labels))
  
  if (use_ruv) {
    design <- model.matrix(~ sample_type, data = pData(set))
    y <- DGEList(counts = counts(set), group = set$sample_type)
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = 2)
    top <- topTags(lrt, n = nrow(set))$table
    top_genes <- rownames(top)[1:min(n_topGenes, nrow(top))]
    empirical <- setdiff(rownames(set), top_genes)
    if (length(empirical) < 50) {message("Skipping ", rds_file, " (too few empirical control genes for RUV)")
      return(NULL)
    }
    set1 <- RUVg(set, empirical, k = k)
    pheno <- pData(set1)
    ruv_terms <- paste0("W_", seq_len(k))
    formula_str <- paste("~ sample_type +", paste(ruv_terms, collapse = " + "))
    design <- model.matrix(as.formula(formula_str), data = pheno)
    dge <- DGEList(counts = counts(set1), group = set1$sample_type)
    
  } else {
    design <- model.matrix(~ sample_type, data = pData(set))
    dge <- DGEList(counts = counts(set), group = set$sample_type)
  }
  
  dge <- calcNormFactors(dge, method = "TMM")
  dge <- estimateDisp(dge, design)
  fit <- glmFit(dge, design)
  lrt <- glmLRT(fit, coef = 2)
  res <- lrt$table
  colnames(res) <- c("log2FoldChange", "logCPM", "LR", "pvalue")
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res$identifier <- rownames(res)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(rds_file)), ".csv"))
  write.csv(res, out_file, row.names = FALSE)
  message("Saved DE results: ", out_file)
  return(out_file)
}
# do_DGE(rds_file = 'data/treehouse/exp_data_gtex/Acute-Megakaryoblastic-Leukemia_SJC_GTEx.rds', out_dir = 'output_files/dge_treehouse', use_ruv = T, k = 1, n_topGenes = 2000)
# Function: run DE for all RDS files in a directory
do_DGE_for_directory <- function(input_dir, output_dir, use_ruv = FALSE, k = 1) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(output_dir, "00_dge_log.txt")
  if (file.exists(log_file)) file.remove(log_file)
  
  rds_files <- list.files(input_dir, full.names = TRUE, pattern = "\\.rds$")
  for (f in rds_files) {
    tryCatch({
      msg <- capture.output(do_DGE(f, out_dir = output_dir, use_ruv = use_ruv, k = k), type = "message"); cat(msg, "\n", file = log_file, append = TRUE)
    }, error = function(e) {msg <- sprintf("ERROR: %s | %s", f, e$message); message(msg); cat(msg, "\n", file = log_file, append = TRUE)})
  }
  message(sprintf("Log written to: %s", log_file))
}
do_DGE_for_directory(input_dir = "data/treehouse/exp_data_gtex", output_dir = "output_files/dge_results_treehouse/unmatched", use_ruv = TRUE, k = 1)
#####

##### QC I - Data science approach #####
summarize_dge <- function(file) {
  res <- read.csv(file)
  n_genes <- nrow(res); sig <- res %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5); sig_genes <- nrow(sig)
  prop_pval <- mean(res$pvalue < 0.01, na.rm = TRUE); mean_abs_logFC <- mean(abs(res$log2FoldChange), na.rm = TRUE)
  snr <- if (sig_genes > 0) mean(abs(sig$log2FoldChange), na.rm = TRUE) else 0
  up <- sum(sig$log2FoldChange > 0.5); down <- sum(sig$log2FoldChange < -0.5)
  
  tibble(file = basename(file), n_genes = n_genes, sig_genes = sig_genes,
         prop_pval = prop_pval, mean_abs_logFC = mean_abs_logFC, snr = snr, up = up, down = down)
}
de_res_dir <- "output_files/dge_results_treehouse/unmatched"
dge_files <- list.files(de_res_dir, full.names = TRUE, pattern = "\\.csv$")
qc_summary <- bind_rows(lapply(dge_files, summarize_dge)) %>% mutate(Disease = sub("_.*", "", file))
qc_summary <- qc_summary %>% mutate(Source = str_extract(file, "_(.*?)_GTEx\\.csv$") %>% str_remove_all("^_|_GTEx\\.csv$")); table(qc_summary$Source)
head(qc_summary); write.csv(qc_summary, file.path(de_res_dir,'../', "00_qc_summary.csv"), row.names = FALSE)

# hist(qc_summary$sig_genes)
# hist(qc_summary$prop_pval)
# hist(qc_summary$mean_abs_logFC)
# hist(qc_summary$snr)

high_quality <- qc_summary %>% filter(sig_genes >= 10 & sig_genes <= 10000, prop_pval > 0.01 & prop_pval < 0.5, 
                                      mean_abs_logFC > 0.1 & mean_abs_logFC < 3, snr > 0.8 & snr < 5, up > 0 & down > 0)
head(high_quality); write.csv(high_quality, "output_files/dge_results_treehouse/01_high_quality_DE.csv", row.names = FALSE)

low_quality <- anti_join(qc_summary, high_quality)
head(low_quality); write.csv(low_quality, "output_files/dge_results_treehouse/02_low_quality_DE.csv", row.names = FALSE)
#####
######################### Treehouse - END #########################