setwd('Desktop/Research/binchenlab/disease_signature/')

##### Required libraries #####
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)
#####

##### EnrichR internal evaluation #####
# parse_row <- function(row) {
#   tokens <- strsplit(paste(row, collapse = "\t"), "\t")[[1]]
#   tokens <- tokens[tokens != ""]
#   disease_gse <- tokens[1]
#   gse <- str_extract(disease_gse, "GSE[0-9]+")
#   disease <- str_trim(sub(gse, "", disease_gse))
#   rest <- tokens[-1]
#   is_num <- !is.na(suppressWarnings(as.numeric(rest[1])))
#   if (is_num) {
#     values <- as.numeric(rest[c(TRUE, FALSE)])
#     genes  <- rest[c(FALSE, TRUE)]
#   } else {
#     genes  <- rest[c(TRUE, FALSE)]
#     values <- as.numeric(rest[c(FALSE, TRUE)])
#   }
#   if (length(genes) > length(values)) values <- c(values, NA_real_)
#   tibble(Disease = disease, GSE = gse, Gene = genes, Value = values)
# }
#
# up_data <- read.csv("output_files/signature_for_drug_discovery/validation/enrichr/Disease_Signatures_from_GEO_up_2014.txt", stringsAsFactors = FALSE)
# up_df <- map_dfr(1:nrow(up_data), ~ parse_row(unlist(up_data[.x, ])))
# up_list <- up_df %>% split(.$Disease) %>% map(~ split(.x, .$GSE)) %>% map(~ map(.x, ~ setNames(.x$Value, .x$Gene)))
#
# dn_data <- read.csv("output_files/signature_for_drug_discovery/validation/enrichr/Disease_Signatures_from_GEO_down_2014.txt", stringsAsFactors = FALSE)
# dn_df <- map_dfr(1:nrow(dn_data), ~ parse_row(unlist(dn_data[.x, ])))
# dn_list <- dn_df %>% split(.$Disease) %>% map(~ split(.x, .$GSE)) %>% map(~ map(.x, ~ setNames(.x$Value, .x$Gene)))
#
# sig_list <- union(names(up_list), names(dn_list)) %>% set_names() %>%
#   map(function(disease) {
#     gses <- union(names(up_list[[disease]]), names(dn_list[[disease]]))
#     set_names(gses) %>% map(function(gse) {
#       list(up   = if (!is.null(up_list[[disease]][[gse]])) up_list[[disease]][[gse]] else NULL, down = if (!is.null(dn_list[[disease]][[gse]])) dn_list[[disease]][[gse]] else NULL)
#       })
#   })
# str(sig_list, max.level = 3)
# saveRDS(sig_list, file = 'output_files/signature_for_drug_discovery/validation/enrichr/sig_enrichr.rds')
sig_enrichr <- readRDS('output_files/signature_for_drug_discovery/validation/enrichr/sig_enrichr.rds')
# length(intersect(names(sig_enrichr$`Acute Lung Injury`$GSE10474$up), names(sig_enrichr$`Acute Lung Injury`$GSE1871$up)))
# length(intersect(names(sig_enrichr$`Acute Lung Injury`$GSE10474$down), names(sig_enrichr$`Acute Lung Injury`$GSE1871$down)))
#
# length(intersect(names(sig_enrichr$`Breast Cancer`$GSE2429$up), names(sig_enrichr$`Breast Cancer`$GSE3744$up)))
# length(intersect(names(sig_enrichr$`Breast Cancer`$GSE2429$down), names(sig_enrichr$`Breast Cancer`$GSE3744$down)))
correlation_for_direction <- function(disease_list, type = "up") {
  gses <- names(disease_list)
  combn(gses, 2, simplify = FALSE) %>% map_dfr(function(pair) {
    g1 <- disease_list[[pair[1]]][[type]]
    g2 <- disease_list[[pair[2]]][[type]]
    common_genes <- intersect(names(g1), names(g2))
    if (length(common_genes) > 1) {
      vals1 <- g1[common_genes]; vals2 <- g2[common_genes]; r <- cor(vals1, vals2, use = "complete.obs")
    } else {r <- NA_real_}
    tibble(GSE1 = pair[1], GSE2 = pair[2], Direction = type, n_common = length(common_genes), Correlation = r, CommonGenes = paste(common_genes, collapse = ","))
  })
}
sig_enrichr_corr_table <- map_dfr(names(sig_enrichr), function(disease) {
  disease_list <- sig_enrichr[[disease]]
  if (length(disease_list) >= 2) {
    bind_rows(correlation_for_direction(disease_list, "up"), correlation_for_direction(disease_list, "down")) %>% mutate(Disease = disease, .before = 1)
  }
})
# write.csv(sig_enrichr_corr_table, file = 'output_files/signature_for_drug_discovery/validation/enrichr/sig_enrichr_corr_table.csv')
#####


read_enrichr_perturbations <- function(up_path, down_path) {
  parse_file <- function(data, direction) {
    lines <- apply(data, 1, paste, collapse = " ")
    lines <- str_squish(lines)
    m <- str_match(lines, "^(.*?)\\s+(human|mouse|rat)\\s+(GSE\\d+)\\s+sample\\s*(\\d+)\\s+(.*)$")
    tibble(Disease = m[, 2], Species = m[, 3], GSE = m[, 4], Genes = str_replace_all(m[, 6], "\\s+", ","), Direction = direction) %>%
      filter(!is.na(GSE)) %>% mutate(Disease_ID = str_extract(Disease, "(DOID-[0-9]+|C[0-9]{4,}|[0-9]{3,6})"), Disease = str_trim(str_remove(Disease, "(DOID-[0-9]+|C[0-9]{4,}|[0-9]{3,6})")))
  }
  up_data <- read.delim(up_path, header = FALSE, stringsAsFactors = FALSE); dn_data <- read.delim(down_path, header = FALSE, stringsAsFactors = FALSE)
  up_df <- parse_file(up_data, "Up"); dn_df <- parse_file(dn_data, "Down")
  all_df <- bind_rows(up_df, dn_df) %>% dplyr::select(Disease, Disease_ID, Species, GSE, Direction, Genes)
  all_df <- all_df %>% group_by(Disease, Disease_ID, Species, GSE) %>%
    summarise(Direction = ifelse(n_distinct(Direction) > 1, "Both", first(Direction)), Genes = paste(unique(unlist(str_split(paste(Genes, collapse = ","), ","))), collapse = ","), .groups = "drop")
  return(all_df)
}
sig_all <- read_enrichr_perturbations("output_files/signature_for_drug_discovery/validation/enrichr/Disease_Perturbations_from_GEO_up.txt",
                                      "output_files/signature_for_drug_discovery/validation/enrichr/Disease_Perturbations_from_GEO_down.txt")
head(sig_all, 2)
# saveRDS(sig_all, file = 'output_files/signature_for_drug_discovery/validation/enrichr/sig_enrichr_big.rds')
sig_human <- sig_all %>% filter(Species == "human")
calc_jaccard <- function(g1, g2) {
  genes1 <- unique(trimws(unlist(str_split(g1, ",")))); genes2 <- unique(trimws(unlist(str_split(g2, ","))))
  inter <- length(intersect(genes1, genes2)); union <- length(union(genes1, genes2))
  if (union == 0) return(NA_real_)
  inter / union
}

jaccard_results <- sig_human %>% group_by(Disease) %>% group_modify(~{
  df <- .x
  gses <- unique(df$GSE)
  if (length(gses) < 2) return(tibble())
  pairs <- expand.grid(GSE1 = gses, GSE2 = gses, stringsAsFactors = FALSE) %>% filter(GSE1 < GSE2)
  pairs %>% mutate(Disease = df$Disease[1], Jaccard = map2_dbl(GSE1, GSE2, ~{
    genes1 <- df$Genes[df$GSE == .x]; genes2 <- df$Genes[df$GSE == .y]; calc_jaccard(genes1, genes2)
    }))
  }) %>% ungroup() %>% dplyr::select(Disease, GSE1, GSE2, Jaccard)
head(jaccard_results, 2)
range(jaccard_results$Jaccard)
length(unique(jaccard_results$Disease))
jaccard_summary <- jaccard_results %>% group_by(Disease) %>% summarise(mean_jaccard = mean(Jaccard, na.rm = TRUE)) %>% arrange(desc(mean_jaccard))
ggplot(jaccard_summary, aes(x = reorder(Disease, mean_jaccard), y = mean_jaccard)) + geom_col(fill = "steelblue") + coord_flip() +
  labs(title = "Overlap per Disease", x = "Disease", y = "Jaccard Index") + theme_minimal(base_size = 10)


sig_geo_meta <- read.csv('output_files/signature_for_drug_discovery/sig_metadata_matched.csv', row.names = 1)
head(sig_geo_meta, 2)

exp_meta <- read.csv('data/matched_disease_info_invivo.csv')
sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
exp_meta$Disease <- sanitize_disease(exp_meta$Disease)
head(exp_meta, 2)

sig_geo_meta$GSE_ID <- sub(".*_(GSE[0-9]+)\\.csv$", "\\1", sig_geo_meta$file)
sig_geo_meta$Control_type <- str_extract(sig_geo_meta$file, "(?<=_)[^_]+(?=_GSE[0-9]+)")
sig_geo_meta$OrganRegion <- sig_geo_meta$file %>% str_replace(paste0("_", sig_geo_meta$Control_type, "_GSE[0-9]+.*"), "") %>% sub("^[^_]+_", "", .)

sanitize_region <- function(x) {x %>% str_trim() %>% gsub("[^A-Za-z0-9]+", "-", .) %>% gsub("^-|-$", "", .)}
sig_geo_meta$OrganRegion <- sanitize_region(sig_geo_meta$OrganRegion)
exp_meta$OrganRegion     <- sanitize_region(exp_meta$OrganRegion)

merged_meta <- merge(sig_geo_meta, exp_meta, by = c('GSE_ID', 'Disease', 'OrganRegion', 'Control_type'))#, all.x = TRUE)
length(unique(merged_meta$Disease))
merged_meta <- merged_meta %>% filter(!is.na(n_cases) & !is.na(n_controls))
merged_meta <- merged_meta %>% na.omit()
length(unique(merged_meta$Disease))
head(merged_meta, 2)
colnames(merged_meta)
merged_meta <- merged_meta %>% dplyr::select(-exp_type, -n_genes, -sig_genes, -prop_pval, -mean_abs_logFC, -snr, -up, -down)
geo_meta <- merged_meta
head(geo_meta, 2)

geo_sig_path <- 'output_files/signature_for_drug_discovery/signatures_all/'
list.files(geo_sig_path)[1:5] # one directory might have multiple csv files for that disease
example <- read.csv('output_files/signature_for_drug_discovery/signatures_all/Abortion-Habitual/Abortion-Habitual_Endometrium_GTEx_GSE65099.csv')
head(example)
mapping <- read.csv('../gene_info_biomart.csv')
head(mapping, 2)

sig_enrichr <- readRDS('output_files/signature_for_drug_discovery/validation/enrichr/sig_enrichr_big.rds')
head(sig_enrichr)
sig_enrichr <- sig_enrichr %>% filter(Species == 'human')
head(sig_enrichr, 2)

dsa_meta <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/Disease_information_Datasets.csv')
dsa_meta <- dsa_meta %>% filter(library_strategy == 'Microarray') %>% separate(control_case_sample_count, into = c("control", "case"), sep = "\\|", convert = TRUE) %>% filter(control > 3, case > 3)
head(dsa_meta, 2)

length(intersect(dsa_meta$accession, sig_enrichr$GSE))
common_gse <- intersect(dsa_meta$accession, sig_enrichr$GSE)

dsa_common <- dsa_meta %>% filter(accession %in% common_gse)
enrichr_common <- sig_enrichr %>% filter(GSE %in% common_gse)

dsaid_to_gse <- setNames(dsa_common$accession, dsa_common$dsaid)
dsa_sig <- readRDS('output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_signature.rds')

head(names(dsa_sig))
head(names(dsaid_to_gse))
intersect(names(dsa_sig), names(dsaid_to_gse))

names(dsa_sig) <- dsaid_to_gse[names(dsa_sig)]
dsa_sig <- dsa_sig[!is.na(names(dsa_sig))]
dsa_sig





library(dplyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)

# --- Helper: Jaccard ---
calc_jaccard <- function(g1, g2) {
  inter <- intersect(g1, g2)
  union_len <- length(union(g1, g2))
  if (union_len == 0) return(list(jaccard = NA, inter = character()))
  list(jaccard = length(inter) / union_len, inter = inter)
}

# --- Load mapping ---
mapping <- read.csv("../gene_info_biomart.csv") %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()

# --- Prepare Enrichr (human only) ---
enrichr_data <- sig_enrichr %>%
  dplyr::select(Disease, Genes) %>%
  mutate(Genes = str_split(Genes, ",")) %>%
  group_by(Disease) %>%
  summarise(Genes = list(unique(unlist(Genes)))) %>%
  ungroup()

# --- Prepare local signatures ---
geo_sig_path <- "output_files/signature_for_drug_discovery/signatures_all/"
disease_dirs <- list.dirs(geo_sig_path, recursive = FALSE, full.names = FALSE)

local_data <- map_dfr(disease_dirs, function(disease_name) {
  csv_files <- list.files(file.path(geo_sig_path, disease_name),
                          pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) return(NULL)
  
  local_genes <- map_dfr(csv_files, ~ read.csv(.x) %>% dplyr::select(identifier)) %>%
    left_join(mapping, by = c("identifier" = "gene_id")) %>%
    filter(!is.na(gene_name)) %>%
    pull(gene_name) %>%
    unique()
  
  tibble(Disease = disease_name, Genes = list(local_genes))
})

# --- Match diseases by best string match ---
# normalize names (lowercase, remove punctuation)
normalize <- function(x) {
  x %>%
    tolower() %>%
    gsub("[-_/]", " ", .) %>%        # replace -, _, / with space
    gsub("[[:punct:]]", "", .) %>%   # remove other punctuation
    gsub("\\s+", " ", .) %>%         # collapse multiple spaces
    trimws()                         # trim leading/trailing spaces
}

local_data <- local_data %>% mutate(Disease_clean = normalize(Disease))
enrichr_data <- enrichr_data %>% mutate(Disease_clean = normalize(Disease))

common_diseases <- intersect(local_data$Disease_clean, enrichr_data$Disease_clean)

# --- Compute overlaps ---
disease_overlap <- map_dfr(common_diseases, function(disease_key) {
  local_genes <- local_data$Genes[local_data$Disease_clean == disease_key][[1]]
  enrichr_genes <- enrichr_data$Genes[enrichr_data$Disease_clean == disease_key][[1]]
  
  j <- calc_jaccard(local_genes, enrichr_genes)
  
  tibble(
    Disease = local_data$Disease[local_data$Disease_clean == disease_key][1],
    Overlap_Gene_Count = length(j$inter),
    Jaccard_Index = round(j$jaccard, 4),
    Overlap_Genes = paste(j$inter, collapse = ", ")
  )
})

hist(disease_overlap$Jaccard_Index)
hist(disease_overlap$Overlap_Gene_Count)



length(intersect(sig_all$GSE, dsa_meta$accession))






# ##### DiSignAtlas internal evaluation #####
# meta_1 has dataset info, meta_2 has sample level info (GSM ids)
meta_1 <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/Disease_information_Datasets.csv', sep = ',')
# head(meta_1, 5)[, 1:5]
# colnames(meta_1)
# table(meta_1$organism)
# table(meta_1$library_strategy)
# table(meta_1$data_source)
# head(meta_1$control_case_sample_count)
# meta_1 <- meta_1 %>% filter(deg_count > 10, organism == "Homo sapiens", library_strategy == 'RNA-Seq') %>%
#   separate(control_case_sample_count, into = c("control", "case"), sep = "\\|", convert = TRUE) %>% filter(control > 3, case > 3)
meta_1 <- meta_1 %>% separate(control_case_sample_count, into = c("control", "case"), sep = "\\|", convert = TRUE) %>% filter(control > 3, case > 3)
# meta_1 <- meta_1 %>% filter(!is.na(tissue), tissue != "") # Data from TCGA doesn't have tissue labelled so don't remove blank tissue values
meta_1 <- meta_1 %>% dplyr::rename(tissue_1 = tissue)

parse_gmt_meta <- function(gmt_file) {
  gmt <- read.delim(gmt_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(gmt)[1:2] <- c("dataset_id", "info")

  rows <- map_dfr(seq_len(nrow(gmt)), function(i) {
    tokens <- strsplit(gmt$info[i], "\\|")[[1]]
    deg_index <- which(!is.na(suppressWarnings(as.integer(tokens))))[1]
    tibble(dataset_id = gmt$dataset_id[i], accession = tokens[1], control_accession = gsub(";", ",", tokens[2]), case_accession = gsub(";", ",", tokens[3]),
           deg_count = as.integer(tokens[deg_index]), disease = tokens[deg_index + 1], tissue = tokens[deg_index + 3], source = tokens[deg_index + 4],
           strategy = tokens[deg_index + 5], organism = tokens[deg_index + 6])})
  rows
}
meta_2 <- parse_gmt_meta("output_files/signature_for_drug_discovery/validation/disignatlas/Disease_information_DEGs.gmt")
meta_2 <- meta_2 %>% dplyr::rename(tissue_2 = tissue)

meta <- meta_1 %>% inner_join(meta_2, by = c("dsaid" = "dataset_id", "accession" = "accession", "deg_count" = "deg_count", "disease" = "disease",
                                             "data_source" = "source", "library_strategy" = "strategy", "organism" = "organism"))
colnames(meta)
length(unique(meta$dsaid)); length(unique(meta$accession)); length(unique(meta$platform)); length(unique(meta$disease))
range(meta$deg_count); table(meta$organism); table(meta$library_strategy); table(meta$data_source)
# write.csv(meta, file = 'output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_metadata.csv', row.names = F)
# 
# 
dsa_meta <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_metadata.csv')
# 
dsa_sig_path <- 'output_files/signature_for_drug_discovery/validation/disignatlas/dsa_diff_download'
dsa_sig_files <- list.files(dsa_sig_path, full.names = T)
head(dsa_sig_files)

file_ids <- str_extract(basename(dsa_sig_files), "DSA\\d+")
files_df <- tibble(dataset_id = file_ids, file_path = dsa_sig_files)
matched_files <- dsa_meta %>% dplyr::rename(dataset_id = dsaid) %>% inner_join(files_df, by = "dataset_id")
# test <- read.csv(matched_files$file_path[1], sep = '\t')
# head(test)

id_to_anno <- function(ids) {
  tibble(GeneID = ids, Symbol = mapIds(org.Hs.eg.db, keys = as.character(ids), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"),
         Ensembl = mapIds(org.Hs.eg.db, keys = as.character(ids), keytype = "ENTREZID", column = "ENSEMBL", multiVals = "first"))
}
read_dsa_df <- function(file, dataset_id) {
  df <- readr::read_delim(file, delim = "\t", show_col_types = FALSE) %>% as.data.frame()
  df <- dplyr::select(df, GeneID, Log2FC, PValue, AdjPValue)
  anno <- id_to_anno(df$GeneID)
  df <- df %>% left_join(anno, by = "GeneID") %>% dplyr::select(GeneID, Symbol, Ensembl, Log2FC, PValue, AdjPValue)
  return(df)
}
dataset_dfs <- purrr::map2(matched_files$file_path, matched_files$dataset_id, read_dsa_df)
names(dataset_dfs) <- matched_files$dataset_id
head(dataset_dfs[[1]])
saveRDS(dataset_dfs, file = 'output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_signature.rds')

dsa_sig <- readRDS('output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_signature.rds')
head(dsa_sig[[1]])
head(names(dsa_sig))

dsa_meta <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_metadata.csv')
head(dsa_meta, 2)

multi_disease <- dsa_meta %>% group_by(disease) %>% filter(n() > 2) %>% pull(disease) %>% unique()
head(multi_disease, 5)

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

corr_analysis_disease <- function(disease, dsa_sig, dsa_meta) {
  dsas <- dsa_meta %>% filter(disease == !!disease) %>% pull(dsaid)
  if (length(dsas) < 2) {
    message("Not enough datasets for: ", disease)
    return(NULL)
  }
  pairs <- combn(dsas, 2, simplify = FALSE)
  results <- lapply(pairs, function(pair) {
    d1 <- dsa_sig[[pair[1]]] %>% dplyr::select(Symbol, Log2FC) %>% filter(!is.na(Symbol))
    d2 <- dsa_sig[[pair[2]]] %>% dplyr::select(Symbol, Log2FC) %>% filter(!is.na(Symbol))
    merged <- inner_join(d1, d2, by = "Symbol", suffix = c("_1", "_2"))
    if (nrow(merged) < 5) return(NULL)
    cor_val <- cor(merged$Log2FC_1, merged$Log2FC_2, use = "pairwise.complete.obs")
    list(pair = pair, cor = cor_val, n_common = nrow(merged))
  })

  results <- Filter(Negate(is.null), results)
  if (length(results) == 0) {
    message("No valid comparisons for: ", disease)
    return(NULL)
  }
  corr_table <- data.frame(Dataset1 = sapply(results, function(x) x$pair[1]), Dataset2 = sapply(results, function(x) x$pair[2]),
                           N_common = sapply(results, function(x) x$n_common), Correlation = sapply(results, function(x) x$cor))
  dsa2meta <- dsa_meta %>% dplyr::select(dsaid, accession, control_accession, case_accession)
  corr_table <- corr_table %>% left_join(dsa2meta, by = c("Dataset1" = "dsaid")) %>%
    dplyr::rename(GSE1 = accession, Control1 = control_accession, Case1 = case_accession) %>%
    left_join(dsa2meta, by = c("Dataset2" = "dsaid")) %>% dplyr::rename(GSE2 = accession, Control2 = control_accession, Case2 = case_accession) %>%
    dplyr::select(Dataset1, GSE1, Control1, Case1, Dataset2, GSE2, Control2, Case2, N_common, Correlation)
  corr_matrix <- corr_table %>% dplyr::select(Dataset1, Dataset2, Correlation) %>% pivot_wider(names_from = Dataset2, values_from = Correlation) %>% as.data.frame()
  rownames(corr_matrix) <- corr_matrix$Dataset1
  corr_matrix <- corr_matrix[, -1, drop = FALSE]
  corr_matrix[is.na(corr_matrix)] <- 0
  diag(corr_matrix) <- 1
  colnames(corr_matrix) <- setdiff(colnames(corr_matrix), "Dataset1")

  dsa_labels <- dsa_meta %>% dplyr::mutate(label = paste0(dsaid, "_", accession)) %>% dplyr::select(dsaid, label)
  rownames(corr_matrix) <- dsa_labels$label[match(rownames(corr_matrix), dsa_labels$dsaid)]
  colnames(corr_matrix) <- dsa_labels$label[match(colnames(corr_matrix), dsa_labels$dsaid)]

  heatmap_obj <- pheatmap(corr_matrix,
                          cluster_rows = if (ncol(corr_matrix) > 1) TRUE else FALSE, cluster_cols = if (ncol(corr_matrix) > 1) TRUE else FALSE,
                          display_numbers = TRUE, main = paste(disease, "corr heatmap"), color = colorRampPalette(c("steelblue", "white", "darkred"))(50),
                          breaks = seq(-1, 1, length.out = 51), angle_col = 45)
  return(list(corr_table = corr_table, corr_matrix = corr_matrix, heatmap = heatmap_obj))
}

# alz_output <- corr_analysis_disease("Alzheimer's Disease", dsa_sig, dsa_meta)
# head(alz_output$corr_table, 2)
# alz_output$corr_matrix[1:5, 1:5]

multi_output <- map(head(multi_disease, 20), ~ corr_analysis_disease(.x, dsa_sig, dsa_meta))
names(multi_output) <- head(multi_disease, 20)
multi_output <- compact(multi_output)
pdf("output_files/signature_for_drug_discovery/validation/disignatlas/top20_disignatlas_sig_correlation_heatmaps.pdf", width = 20, height = 16)
for (disease in names(multi_output)) {
  grid::grid.newpage()
  print(multi_output[[disease]]$heatmap)
}
dev.off()
# saveRDS(multi_output, file = "output_files/signature_for_drug_discovery/validation/disignatlas/top20_disignatlas_sig_correlation_results.rds")
# multi_output <- readRDS("output_files/signature_for_drug_discovery/validation/disignatlas/top20_disignatlas_sig_correlation_results.rds")
#####

##### DrugRepo #####
# https://repo-hub.broadinstitute.org/repurposing#download-data
drug <- read.csv("output_files/signature_for_drug_discovery/repurposing_drugs_20200324.txt", sep = "\t", skip = 9, header = TRUE)
head(drug, 5)

# dim(drug)
# colnames(drug)
# table(drug$clinical_phase)
# length(unique(drug$disease_area))
# length(unique(drug$indication))

drug <- drug %>% mutate(clinical_phase = trimws(as.character(clinical_phase)))
disease_indication_summary <- drug %>% filter(disease_area != "", clinical_phase != "") %>% group_by(disease_area, clinical_phase) %>%
  summarise(n = n(), .groups = "drop") %>% pivot_wider(names_from  = clinical_phase, values_from = n, names_prefix = "n_", values_fill = 0) %>%
  left_join(drug %>% filter(disease_area != "", clinical_phase != "", indication != "") %>% group_by(disease_area) %>%
              summarise(unique_indication = n_distinct(indication), .groups = "drop"), by = "disease_area") %>%
  mutate(total   = rowSums(across(starts_with("n_"))), percent = round(100 * total / sum(total), 2)) %>%
  relocate(total, .after = disease_area) %>% relocate(c(unique_indication, percent), .after = last_col()) %>% arrange(desc(total))

write.csv(drug, file = "output_files/signature_for_drug_discovery/validation/drug_info.csv", row.names = FALSE)
write.csv(disease_indication_summary, file = "output_files/signature_for_drug_discovery/validation/disease_indication_summary.csv", row.names = FALSE)

# neuro <- drug %>% filter(disease_area == 'neurology/psychiatry')
# neuro_pairs <- neuro %>% select(pert_iname, indication, clinical_phase) %>% filter(indication != "") %>% arrange(indication, pert_iname)
# write.csv(neuro_pairs, file = "output_files/signature_for_drug_discovery/validation/neuro_pairs.csv", row.names = FALSE)
#####

##### External evaluation #####
sig_geo_meta <- read.csv('output_files/signature_for_drug_discovery/sig_metadata_matched.csv', row.names = 1)
head(sig_geo_meta, 2)

exp_meta <- read.csv('data/matched_disease_info_invivo.csv')
sanitize_disease <- function(x) gsub("[^A-Za-z0-9]+", "-", x)
exp_meta$Disease <- sanitize_disease(exp_meta$Disease)
head(exp_meta, 2)

sig_geo_meta$GSE_ID <- sub(".*_(GSE[0-9]+)\\.csv$", "\\1", sig_geo_meta$file)
sig_geo_meta$Control_type <- str_extract(sig_geo_meta$file, "(?<=_)[^_]+(?=_GSE[0-9]+)")
sig_geo_meta$OrganRegion <- sig_geo_meta$file %>% str_replace(paste0("_", sig_geo_meta$Control_type, "_GSE[0-9]+.*"), "") %>% sub("^[^_]+_", "", .)

sanitize_region <- function(x) {x %>% str_trim() %>% gsub("[^A-Za-z0-9]+", "-", .) %>% gsub("^-|-$", "", .)}
sig_geo_meta$OrganRegion <- sanitize_region(sig_geo_meta$OrganRegion)
exp_meta$OrganRegion     <- sanitize_region(exp_meta$OrganRegion)

merged_meta <- merge(sig_geo_meta, exp_meta, by = c('GSE_ID', 'Disease', 'OrganRegion', 'Control_type'))#, all.x = TRUE)
length(unique(merged_meta$Disease))
merged_meta <- merged_meta %>% filter(!is.na(n_cases) & !is.na(n_controls))
merged_meta <- merged_meta %>% na.omit()
length(unique(merged_meta$Disease))
head(merged_meta, 2)
colnames(merged_meta)
merged_meta <- merged_meta %>% dplyr::select(-exp_type, -n_genes, -sig_genes, -prop_pval, -mean_abs_logFC, -snr, -up, -down)
geo_meta <- merged_meta
head(geo_meta, 2)

sig_enrichr <- readRDS('output_files/signature_for_drug_discovery/validation/enrichr/sig_enrichr.rds')
# length(names(sig_enrichr))
# head(sig_enrichr$`Acute Lung Injury`$GSE10474$up)
# head(sig_enrichr$`Acute Lung Injury`$GSE10474$down)
# head(sig_enrichr$`Acute Lung Injury`$GSE1871$up)
# head(sig_enrichr$`Acute Lung Injury`$GSE1871$down)
metadata_list <- list()
for (disease in names(sig_enrichr)) {
  disease_data <- sig_enrichr[[disease]]
  for (gse in names(disease_data)) {
    gse_data <- disease_data[[gse]]
    metadata_list[[length(metadata_list) + 1]] <- data.frame(Disease = disease, GSE = gse, stringsAsFactors = FALSE)
  }
}
enrichr_meta <- do.call(rbind, metadata_list)
head(enrichr_meta)
enrichr_meta$Disease <- sanitize_disease(enrichr_meta$Disease)
length(unique(enrichr_meta$GSE))
length(unique(enrichr_meta$Disease))
write.csv(enrichr_meta, file = 'output_files/signature_for_drug_discovery/validation/enrichr/sig_enrichr_metadata.csv')

length(intersect(unique(enrichr_meta$GSE), unique(geo_meta$GSE_ID)))
length(intersect(unique(enrichr_meta$Disease), unique(geo_meta$Disease)))
merged_meta <- merge(enrichr_meta, geo_meta, by = 'Disease')
length(unique(merged_meta$Disease))

# geo_sig_example <- read.csv('output_files/signature_for_drug_discovery/signatures_all/Acute-Lung-Injury/Acute-Lung-Injury_Nasal_Normal_GSE192364.csv')
# head(geo_sig_example)
# enrichr_up_sig_example <- sig_enrichr$`Acute Lung Injury`$GSE10474$up
# enrichr_dn_sig_example <- sig_enrichr$`Acute Lung Injury`$GSE10474$dn
# head(enrichr_up_sig_example)
# head(enrichr_dn_sig_example)

library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
map_ids <- function(ensembl_ids) {
  mapping <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filters="ensembl_gene_id", values=ensembl_ids, mart=mart)
  mapping
}
compute_signature_corr <- function(disease, geo_file, enr_up, enr_dn) {
  geo_sig <- read.csv(geo_file)
  mapping <- map_ids(geo_sig$identifier)
  geo_sig <- merge(geo_sig, mapping, by.x="identifier", by.y="ensembl_gene_id")
  geo_vec <- setNames(geo_sig$log2FoldChange, geo_sig$hgnc_symbol)
  enr_up <- enr_up; enr_dn <- -1 * enr_dn; enr_vec <- c(enr_up, enr_dn)
  common_genes <- intersect(names(geo_vec), names(enr_vec))
  if (length(common_genes) < 5) {return(data.frame(Disease=disease, Correlation=NA, N=length(common_genes)))}
  geo_aligned <- geo_vec[common_genes]
  enr_aligned <- enr_vec[common_genes]
  cor_val <- cor(geo_aligned, enr_aligned, method="pearson")
  data.frame(Disease=disease, Correlation=cor_val, N=length(common_genes))
}

merged_meta <- merged_meta %>% mutate(geo_path = file.path("output_files/signature_for_drug_discovery/signatures_all", Disease, file))
head(merged_meta, 2)

results <- data.frame()
results <- do.call(rbind, lapply(unique(merged_meta$Disease), function(d) {
  geo_files <- unique(merged_meta$geo_path[merged_meta$Disease == d])
  enr_up <- sig_enrichr[[d]][[1]]$up; enr_dn <- sig_enrichr[[d]][[1]]$dn
  do.call(rbind, lapply(geo_files, function(f) {compute_signature_corr(d, f, enr_up, enr_dn)}))}))
head(results)
write.csv(results, file = 'output_files/signature_for_drug_discovery/validation/enrichr/corr_results.csv')

dsa_meta <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_metadata.csv')
head(dsa_meta, 2)
colnames(dsa_meta)
dsa_meta <- dsa_meta %>% dplyr::select(-platform, -deg_count, -diseaseid, -library_strategy, -data_source, -tissue_1, -tissue_2, -organism, -definition)
dsa_sig <- readRDS('output_files/signature_for_drug_discovery/validation/disignatlas/disignatlas_signature.rds')

normalize_gsm <- function(x) {x %>% str_split(",") %>% unlist() %>% str_trim() %>% sort() %>% paste(collapse = ",")}
merged_meta$cases_gsm    <- sapply(merged_meta$cases_gsm, normalize_gsm)
merged_meta$controls_gsm <- sapply(merged_meta$controls_gsm, normalize_gsm)
dsa_meta$case_accession    <- sapply(dsa_meta$case_accession, normalize_gsm)
dsa_meta$control_accession <- sapply(dsa_meta$control_accession, normalize_gsm)


final_merged <- merged_meta %>% inner_join(dsa_meta, by = c("GSE_ID" = "accession"), relationship = 'many-to-many')
# final_merged <- merged_meta %>% inner_join(dsa_meta, by = c("GSE_ID" = "accession", "cases_gsm" = "case_accession"), relationship = 'many-to-many')
# final_merged <- merged_meta %>% inner_join(dsa_meta, by = c("GSE_ID" = "accession", "controls_gsm" = "control_accession"), relationship = 'many-to-many')
# final_merged <- merged_meta %>% inner_join(dsa_meta, by = c("GSE_ID" = "accession", "cases_gsm" = "case_accession", "controls_gsm" = "control_accession"))
length(unique(final_merged$Disease))
colnames(final_merged)
final_merged <- final_merged %>% dplyr::select(GSE_ID, Disease, OrganRegion, Control_type, disease, dsaid, n_controls, n_cases, control, case,
                                               file, controls_gsm, cases_gsm, control_accession, case_accession)

# final_merged <- final_merged %>% dplyr::select(GSE_ID, Disease, OrganRegion, Control_type, disease, dsaid, n_controls, n_cases, control, case,
#                                                file, controls_gsm, cases_gsm)

unmatched_from_geo <- merged_meta %>% anti_join(dsa_meta, by = c("GSE_ID" = "accession", "cases_gsm" = "case_accession", "controls_gsm" = "control_accession"))
unmatched_from_dsa <- dsa_meta %>% anti_join(merged_meta, by = c("accession" = "GSE_ID", "case_accession" = "cases_gsm", "control_accession" = "controls_gsm"))

final_merged <- final_merged %>% mutate(geo_path = file.path("output_files/signature_for_drug_discovery/signatures_all", Disease, file))
# write.csv(final_merged, file = 'output_files/signature_for_drug_discovery/validation/disignatlas/merged_metadata_for_corr.csv')

final_merged <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/merged_metadata_for_corr.csv', row.names = 1)
head(final_merged, 2)
colnames(final_merged)
names(final_merged)[names(final_merged) == "Disease"]           <- "Disease_GEOMeta"
names(final_merged)[names(final_merged) == "disease"]           <- "Disease_DSA"
names(final_merged)[names(final_merged) == "dsaid"]             <- "DSA_ID"
names(final_merged)[names(final_merged) == "n_controls"]        <- "n_controls_GEOMeta"
names(final_merged)[names(final_merged) == "n_cases"]           <- "n_cases_GEOMeta"
names(final_merged)[names(final_merged) == "control"]           <- "n_controls_DSA"
names(final_merged)[names(final_merged) == "case"]              <- "n_cases_DSA"
names(final_merged)[names(final_merged) == "controls_gsm"]      <- "controls_gsm_GEOMeta"
names(final_merged)[names(final_merged) == "cases_gsm"]         <- "cases_gsm_GEOMeta"
names(final_merged)[names(final_merged) == "control_accession"] <- "controls_gsm_DSA"
names(final_merged)[names(final_merged) == "case_accession"]    <- "cases_gsm_DSA"

results <- list()
for (i in seq_len(nrow(final_merged))) {
  geo_file <- final_merged$geo_path[i]
  DSA_ID   <- final_merged$DSA_ID[i]
  if (!file.exists(geo_file)) {
    message("File not found: ", geo_file, " (Disease: ", final_merged$Disease_GEOMeta[i], ", GSE_ID: ", final_merged$GSE_ID[i], ", DSA_ID: ", DSA_ID, ")")
    next
  }
  if (file.exists(geo_file) && !is.null(dsa_sig[[DSA_ID]])) {
    geo_sig    <- read.csv(geo_file); dsa_sig_df <- dsa_sig[[DSA_ID]]
    merged_sig <- merge(geo_sig, dsa_sig_df, by.x = "identifier", by.y = "Ensembl", all = FALSE)
    if (nrow(merged_sig) > 0) {
      pearson  <- cor(merged_sig$log2FoldChange, merged_sig$Log2FC, method = "pearson",  use = "pairwise.complete.obs")
      spearman <- cor(merged_sig$log2FoldChange, merged_sig$Log2FC, method = "spearman", use = "pairwise.complete.obs")
    } else {pearson <- spearman <- NA}
    meta_row <- final_merged[i, , drop = FALSE]; meta_row$Cor_Pearson <- pearson; meta_row$Cor_Spearman <- spearman; meta_row$n_genes <- nrow(merged_sig)
    results[[i]] <- meta_row
  }
}
cor_results <- do.call(rbind, results)
head(cor_results, 5)
length(unique(cor_results$Disease_GEOMeta))
length(unique(cor_results$Disease_DSA))
length(unique(cor_results$GSE_ID))
write.csv(cor_results, file = 'output_files/signature_for_drug_discovery/validation/disignatlas/cor_results.csv')
cor_results <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/cor_results.csv', row.names = 1)
hist(cor_results$Cor_Pearson)
low_corr <- subset(cor_results, Cor_Pearson < 0.5)
length(unique(low_corr$Disease_GEOMeta))
length(unique(low_corr$Disease_DSA))
length(unique(low_corr$GSE_ID))
write.csv(low_corr, file = 'output_files/signature_for_drug_discovery/validation/disignatlas/low_cor_results.csv')

# Jaccard Overall
jaccard_results <- list()
for (i in seq_len(nrow(final_merged))) {
  geo_file <- final_merged$geo_path[i]
  dsa_id   <- final_merged$dsaid[i]
  if (!file.exists(geo_file) || is.null(dsa_sig[[dsa_id]])) next
  # Load GEO and DSA data
  geo_sig    <- read.csv(geo_file)
  dsa_sig_df <- dsa_sig[[dsa_id]]
  # Extract gene IDs (strip version numbers if present)
  geo_ids <- unique(gsub("\\..*", "", geo_sig$identifier))
  dsa_ids <- unique(gsub("\\..*", "", dsa_sig_df$Ensembl))
  # Overlap and Jaccard
  overlap_genes <- intersect(geo_ids, dsa_ids)
  n_overlap <- length(overlap_genes)
  union_set <- length(union(geo_ids, dsa_ids))
  jaccard <- ifelse(union_set > 0, n_overlap / union_set, NA)
  overlap_str <- if (n_overlap > 0) paste(overlap_genes, collapse = ",") else ""
  
  jaccard_results[[i]] <- data.frame(Disease = final_merged$Disease[i], GSE_ID = final_merged$GSE_ID[i], DSA_ID = dsa_id,
                                     n_geo = length(geo_ids), n_dsa = length(dsa_ids), n_overlap = n_overlap,
                                     Jaccard = jaccard, OverlapGenes = overlap_str, stringsAsFactors = FALSE)
}
jaccard_table <- do.call(rbind, jaccard_results)
head(jaccard_table)
range(jaccard_table$Jaccard)
hist(jaccard_table$Jaccard)
write.csv(jaccard_table, file = 'output_files/signature_for_drug_discovery/validation/disignatlas/jaccard_table.csv')

mean_jaccard   <- mean(jaccard_table$Jaccard, na.rm = TRUE)
median_jaccard <- median(jaccard_table$Jaccard, na.rm = TRUE)
p_jaccard <- ggplot(jaccard_table, aes(x = Jaccard)) + geom_density(fill = "steelblue", alpha = 0.5) + theme_bw(base_size = 10) +
  geom_vline(aes(xintercept = mean_jaccard), color = "darkred", linetype = "dashed", size = 0.7) +
  geom_vline(aes(xintercept = median_jaccard), color = "darkgreen", linetype = "dashed", size = 0.7) +
  labs(title = "Overall Gene Set Overlap (Jaccard Index)", x = "Jaccard Index", y = "Density", 
       caption = paste("Mean (Dark Red) = ", round(mean_jaccard, 3), ", Median (Green) = ", round(median_jaccard, 3), sep = ""))
p_jaccard; ggsave("output_files/signature_for_drug_discovery/validation/disignatlas/jaccard_density.png", p_jaccard, width = 7, height = 5, dpi = 300)

# Jaccard based on cutoff
padj_cutoff <- 0.05
lfc_cutoff  <- 0.5
# head(geo_sig, 2)
# head(dsa_sig$DSA00004, 2)
jaccard_sig_results <- list()
for (i in seq_len(nrow(final_merged))) {
  geo_file <- final_merged$geo_path[i]
  dsa_id   <- final_merged$dsaid[i]
  if (!file.exists(geo_file) || is.null(dsa_sig[[dsa_id]])) next
  # GEO
  geo_sig <- read.csv(geo_file)
  geo_sig <- geo_sig[geo_sig$padj < padj_cutoff & abs(geo_sig$log2FoldChange) > lfc_cutoff, ]
  geo_ids <- unique(gsub("\\..*", "", geo_sig$identifier))
  # DSA
  dsa_sig_df <- dsa_sig[[dsa_id]]
  dsa_sig_df <- dsa_sig_df[!is.na(dsa_sig_df$Ensembl), ]
  if ("AdjPValue" %in% colnames(dsa_sig_df)) {
    dsa_sig_df <- dsa_sig_df[(is.na(dsa_sig_df$AdjPValue) | dsa_sig_df$AdjPValue < padj_cutoff) & abs(dsa_sig_df$Log2FC) > lfc_cutoff, ]
  } else {
    dsa_sig_df <- dsa_sig_df[abs(dsa_sig_df$Log2FC) > lfc_cutoff, ]
  }
  dsa_ids <- unique(gsub("\\..*", "", dsa_sig_df$Ensembl))
  
  # Overlap & Jaccard
  overlap     <- intersect(geo_ids, dsa_ids)
  n_overlap   <- length(overlap)
  union_set   <- length(union(geo_ids, dsa_ids))
  jaccard     <- ifelse(union_set > 0, n_overlap / union_set, NA)
  overlap_str <- if (length(overlap) > 0) paste(overlap, collapse = ",") else ""
  jaccard_sig_results[[i]] <- data.frame(Disease = final_merged$Disease[i], GSE_ID = final_merged$GSE_ID[i], DSA_ID = dsa_id,
                                         n_geo_sig = length(geo_ids), n_dsa_sig = length(dsa_ids), n_overlap = n_overlap, Jaccard = jaccard,
                                         OverlapGenes = overlap_str, stringsAsFactors = FALSE)
}
jaccard_sig_table <- do.call(rbind, jaccard_sig_results)
head(jaccard_sig_table)
range(jaccard_sig_table$Jaccard)
hist(jaccard_sig_table$Jaccard)
write.csv(jaccard_sig_results, file = 'output_files/signature_for_drug_discovery/validation/disignatlas/jaccard_table_sig.csv')

mean_jaccard   <- mean(jaccard_sig_table$Jaccard, na.rm = TRUE)
median_jaccard <- median(jaccard_sig_table$Jaccard, na.rm = TRUE)
p_jaccard <- ggplot(jaccard_sig_table, aes(x = Jaccard)) + geom_density(fill = "steelblue", alpha = 0.5) + theme_bw(base_size = 10) +
  geom_vline(aes(xintercept = mean_jaccard), color = "darkred", linetype = "dashed", size = 0.7) +
  geom_vline(aes(xintercept = median_jaccard), color = "darkgreen", linetype = "dashed", size = 0.7) +
  labs(title = "Significant Gene Set Overlap (Jaccard Index)", x = "Jaccard Index", y = "Density",
       caption = paste("Mean (Dark Red) = ", round(mean_jaccard, 3), ", Median (Green) = ", round(median_jaccard, 3), sep = ""))
p_jaccard; ggsave("output_files/signature_for_drug_discovery/validation/disignatlas/jaccard_density_sig.png", p_jaccard, width = 7, height = 5, dpi = 300)
#####