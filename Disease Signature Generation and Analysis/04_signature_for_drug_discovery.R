setwd('Desktop/Research/binchenlab/disease_signature/')

##### Required libraries #####
library(dplyr)
library(readxl)
library(tibble)
library(purrr)
#####

meta <- read.csv('output_files/dge_results/01_high_quality_DE.csv')

##### Select INVIVO and remove BLOOD signature #####
# meta <- read.csv('output_files/dge_results/01_high_quality_DE.csv')
# meta <- read.csv('output_files/dge_results_unmatched/01_high_quality_DE.csv')
invivo <- meta %>% filter(exp_type == "invivo")
length(unique(meta$Disease))
length(unique(meta$exp_type))


tapply(meta$Disease, meta$exp_type, function(x) length(unique(x)))
dlist <- split(meta$Disease, meta$exp_type)
dlist <- lapply(dlist, unique)
invitro_invivo <- intersect(dlist$invitro, dlist$invivo)
exvivo_invitro <- intersect(dlist$exvivo, dlist$invitro)
exvivo_invivo  <- intersect(dlist$exvivo, dlist$invivo)
common_all <- Reduce(intersect, dlist)
sapply(list(exvivo_invitro = exvivo_invitro, exvivo_invivo = exvivo_invivo, invitro_invivo = invitro_invivo,all_three = common_all), length)

blood_all <- invivo %>% filter(grepl("Blood", file, ignore.case = TRUE))
length(unique(blood_all$Disease))
blood_related <- c("Anemia-Sickle-Cell","Antiphospholipid-Syndrome","Arthritis-Juvenile","Arthritis-Psoriatic","Arthritis-Rheumatoid",
                   "Asthma","Autoimmune-Lymphoproliferative-Syndrome","Cryopyrin-Associated-Periodic-Syndromes","Dermatomyositis",
                   "Familial-Mediterranean-Fever","Giant-Cell-Arteritis","Graft-vs-Host-Disease","Hemorrhagic-Fever-Ebola","Hereditary-Autoinflammatory-Diseases",
                   "HIV-Infections","Hodgkin-Disease","Leukemia-Lymphocytic-Chronic-B-Cell","Leukemia-Myeloid-Acute","Lupus-Erythematosus-Systemic",
                   "Macrophage-Activation-Syndrome","Malaria","Monoclonal-Gammopathy-of-Undetermined-Significance","Multiple-Myeloma","Neoplasms-Plasma-Cell",
                   "Precursor-Cell-Lymphoblastic-Leukemia-Lymphoma","Psoriasis","Rheumatic-Fever","Rheumatoid-Arthritis-Systemic-Juvenile","Sarcoidosis",
                   "Sepsis","Sj-gren-Mikulicz-syndrome","Sjogren-s-Syndrome","Still-s-Disease-Adult-Onset","Systemic-Inflammatory-Response-Syndrome",
                   "Tuberculosis","Tuberculosis-Pulmonary","VEXAS-syndrome",
                   'Atherosclerosis', 'Chickenpox', 'Chronic-recurrent-multifocal-osteomyelitis', 'COVID-19', 'Dengue', 'fetal-inflammatory-response-syndrome',
                   'Virus-Diseases', 'Pneumonia-Bacterial', 'Hepatitis-Chronic', 'Hepatitis-B', 'Hepatitis-C', 'Leishmaniasis', 'Metapneumovirus')
infection_related <- grep("Infection", unique(blood_all$Disease), value = TRUE, ignore.case = TRUE)
blood_related <- unique(c(blood_related, infection_related))
non_blood <- setdiff(unique(blood_all$Disease), blood_related)
# all_diseases <- unique(blood_all$Disease)
# disease_df <- data.frame(Disease = all_diseases,
#                          blood_related = ifelse(all_diseases %in% blood_related, all_diseases, ""),
#                          non_blood = ifelse(all_diseases %in% non_blood, all_diseases, ""))
# write.csv(disease_df, "output_files/signature_for_drug_discovery/disease_classification.csv", row.names = FALSE)

invivo_clean <- invivo %>% filter(!(grepl("Blood", file, ignore.case = TRUE) & !Disease %in% blood_related))
cat('Non-Blood related signature removed: ', (length((invivo$file)) - length((invivo_clean$file))), '\n')
cat('Total Disease: ', length(unique(invivo_clean$Disease)), '\n')
cat('Total Signature: ', length(invivo_clean$file), '\n')

# sig_meta_matched <- invivo_clean
# write.csv(sig_meta_matched, file = 'output_files/signature_for_drug_discovery/sig_metadata_matched.csv')

# sig_meta_unmatched <- invivo_clean
# write.csv(sig_meta_unmatched, file = 'output_files/signature_for_drug_discovery/sig_metadata_unmatched.csv')
#####

##### GEO X GTEX #####
sig_geo_meta <- read.csv('output_files/signature_for_drug_discovery/sig_metadata_matched.csv', row.names = 1)
sig_geo_meta$control_source  <- "GEO"

sig_gtex_meta <- read.csv('output_files/signature_for_drug_discovery/sig_metadata_unmatched.csv', row.names = 1)
sig_gtex_meta$control_source <- "GTEX"
length(sig_gtex_meta$file)
length(unique(sig_gtex_meta$Disease))

sig_merged_meta <- rbind(sig_geo_meta, sig_gtex_meta)
head(sig_merged_meta)
colnames(sig_merged_meta)
# write.csv(sig_merged_meta, file = 'output_files/signature_for_drug_discovery/sig_metadata_merged.csv')

disease_freq <- sig_merged_meta %>% group_by(control_source, Disease) %>% summarise(count = n(), files = paste(file, collapse = ", "), .groups = "drop") %>%
  pivot_wider(names_from = control_source, values_from = c(count, files), values_fill = list(count = 0, files = "")) %>%
  mutate(total = count_GEO + count_GTEX) %>% select(Disease, count_GEO, count_GTEX, total, files_GEO, files_GTEX) %>% arrange(desc(total))
write.csv(disease_freq, file = 'output_files/signature_for_drug_discovery/sig_freq.csv')

aggregate(cbind(n_genes, sig_genes, prop_pval, mean_abs_logFC, snr) ~ control_source, data = sig_merged_meta, mean)
metrics <- c("n_genes", "sig_genes", "prop_pval", "mean_abs_logFC", "snr")
sig_long <- melt(sig_merged_meta[, c("control_source", metrics)], id.vars = "control_source")
ggplot(sig_long, aes(x = control_source, y = value, color = control_source)) + geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) + facet_wrap(~ variable, scales = "free_y") + theme_bw() + labs(title = "Signature metrics", x = "Source", y = "Value")
#####

##### Freq plots #####
disease_freq <- read.csv('output_files/signature_for_drug_discovery/sig_freq.csv')

library(ggplot2); library(dplyr); library(tidyr)
# data <- data.frame(Category = factor(c("GEOMeta", "Matched", "ARCHS4", "QC_I", "QC_II"), levels = c("GEOMeta", "Matched", "ARCHS4", "QC_I", "QC_II")),
#                    Disease = c(623, 405, 375, 323, 265), Dataset   = c(2513, 1141, 984, 772, 638))
# data <- data.frame(Category = factor(c("GEOMeta", "Matched", "ARCHS4_GTEx", "QC_I", "QC_II"), levels = c("GEOMeta", "Matched", "ARCHS4_GTEx", "QC_I", "QC_II")),
#                    Disease = c(623, 323, 175, 126, 107), Dataset   = c(2513, 1529, 399, 251, 227))
df_long <- data %>% pivot_longer(cols = c(Disease, Dataset), names_to = "Source", values_to = "Count")
ggplot(df_long, aes(x = Category, y = Count, group = Source, color = Source)) + geom_line(size = 1) + geom_point(size = 2) +
  geom_text(aes(label = Count), vjust = -0.7, size = 3.5) + scale_color_manual(values = c("Disease" = "darkred", "Dataset" = "steelblue")) +
  labs(title = "Freq. of Disease and Dataset", x = "", y = "") + theme_minimal(base_size = 14) + 
  theme(plot.title = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "top")

# phenoDF <- get_ExperimentHub_data('EH7274')
# phenoDF <- phenoDF %>% filter(data.source == 'GTEX')
# table(phenoDF$biopsy.site)
# ae_controls <- readRDS('data/ae_derived_controls_filtered.rds')
# ae_controls <- ae_controls$invivo

sig_meta <- read.csv('output_files/signature_for_drug_discovery/sig_metadata_merged.csv')
length(sig_meta$file)
table(sig_meta$control_source)
length(unique(sig_meta$Disease))

geo <- sig_meta %>% filter(control_source == 'GEO')
length(unique(geo$Disease))

gtex <- sig_meta %>% filter(control_source == 'GTEX')
length(unique(gtex$Disease))

data <- data.frame(Category = c("Diseases", "Datasets"), GEO = c(265, 638), GTEx = c(107, 227), Total = c(312, 865))
data_long <- tidyr::pivot_longer(data, cols = c("GEO", "GTEx"), names_to = "Source", values_to = "Count")
ggplot(data_long, aes(x = Category, y = Count, fill = Source, label = Count)) + geom_bar(stat = "identity") + 
  geom_text(position = position_stack(vjust = 0.5), color = "white", size = 4) + 
  geom_text(data = data, aes(x = Category, y = Total, label = paste0("", Total)), inherit.aes = FALSE, vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("GEO" = "steelblue", "GTEx" = "firebrick")) + labs(title = "HQ signatures", y = "Count", x = "") +
  theme_minimal(base_size = 10)
#####

##### Metadata #####
sig_geo_meta <- read.csv('output_files/signature_for_drug_discovery/sig_metadata_matched.csv', row.names = 1)
sig_geo_meta$signature_source  <- "GEO"
sig_geo_meta$case_source  <- "GEO"
sig_geo_meta$control_source  <- "GEO"

sig_gtex_meta <- read.csv('output_files/signature_for_drug_discovery/sig_metadata_unmatched.csv', row.names = 1)
sig_gtex_meta$signature_source  <- "GEO"
sig_gtex_meta$case_source  <- "GEO"
sig_gtex_meta$control_source <- "GTEX"

sig_tcga_meta <- read.csv('output_files/dge_results_tcga/00_qc_summary.csv')
sig_tcga_meta$signature_source  <- "TCGA"
sig_tcga_meta$case_source  <- "TCGA"
sig_tcga_meta <- sig_tcga_meta %>% mutate(control_source = case_when(exp_type == "matched" ~ "TCGA", exp_type == "unmatched" ~ "GTEX", TRUE ~ NA_character_))

sig_treehouse_meta <- read.csv('output_files/dge_results_treehouse/00_qc_summary.csv')
sig_treehouse_meta$signature_source  <- "TreeHouse"
sig_treehouse_meta$case_source  <- "TreeHouse"
sig_treehouse_meta$control_source <- "GTEX"
sig_treehouse_meta <- sig_treehouse_meta[!grepl("MM-ICD", sig_treehouse_meta$Disease, fixed = TRUE), ]
any(grepl("MM-ICD", sig_treehouse_meta$Disease, fixed = TRUE))

colnames(sig_geo_meta)
colnames(sig_gtex_meta)
colnames(sig_tcga_meta)
colnames(sig_treehouse_meta)
cols <- unique(c(colnames(sig_geo_meta), colnames(sig_gtex_meta), colnames(sig_tcga_meta), colnames(sig_treehouse_meta)))
cols[!(cols %in% Reduce(intersect, list(colnames(sig_geo_meta), colnames(sig_gtex_meta), colnames(sig_tcga_meta), colnames(sig_treehouse_meta))))]
sig_treehouse_meta <- sig_treehouse_meta %>% rename(exp_type = Source)

sig_meta <- bind_rows(sig_geo_meta %>% mutate(meta_source = "GEO"), sig_gtex_meta %>% mutate(meta_source = "GTEX"), 
                      sig_tcga_meta %>% mutate(meta_source = "TCGA"), sig_treehouse_meta %>% mutate(meta_source = "Treehouse"))
table(sig_meta$meta_source)
table(sig_meta$signature_source)
table(sig_meta$signature_source, sig_meta$control_source)
# table(sig_meta$case_source)
# table(sig_meta$control_source)
# table(sig_meta$exp_type)

# Impute organ
head(sig_meta, 2)
sig_meta$OrganRegion <- NA_character_
is_geo <- sig_meta$signature_source == "GEO"; is_tcga <- sig_meta$signature_source == "TCGA"
is_treehouse<- sig_meta$signature_source == "TreeHouse"
sig_meta$OrganRegion[is_geo] <- sub("^[^_]+_([^_]+(?:_[^_]+)*)_(?i:Normal|AdjacentNormal|Healthy|GTEx).*", "\\1", sig_meta$file[is_geo], perl = TRUE)
sig_meta$OrganRegion[is_tcga] <- sub("^[^_]+_([^_]+(?:_[^_]+)*)_(?i:Adjacent|GTEx)(?:_|\\.|$).*", "\\1", sig_meta$file[is_tcga], perl = TRUE)
sig_meta$OrganRegion[is_treehouse] <- "UA"
unique(sig_meta$OrganRegion)

length(unique(sig_meta$Disease))
freq_table <- sig_meta %>% group_by(Disease) %>% 
  summarise(total_signatures = n(), GEO = sum(signature_source == "GEO"), TCGA = sum(signature_source == "TCGA"), TreeHouse = sum(signature_source == "TreeHouse"),
            control_GEO = sum(control_source == "GEO"), control_TCGA = sum(control_source == "TCGA"), control_GTEX = sum(control_source == "GTEX")) %>% 
  arrange(desc(total_signatures))

# write.csv(sig_meta, file = 'output_files/signature_for_drug_discovery/sig_meta_final.csv')
# write.csv(freq_table, file = 'output_files/signature_for_drug_discovery/disease_freq_sig_meta_final.csv')

sig_meta <- read.csv('output_files/signature_for_drug_discovery/sig_meta_final.csv', row.names = 1)
head(sig_meta, 2)
dim(sig_meta)
sig_meta$Disease <- tools::toTitleCase(tolower(sig_meta$Disease))
length(unique(sig_meta$Disease))
# sig_meta$Disease[grep("Fibrolamellar", sig_meta$Disease)]

output_dir <- "output_files/signature_for_drug_discovery/signature_new"
dir.create(output_dir, showWarnings = FALSE)
all_files <- list.files("output_files", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
file_names <- basename(all_files)
for (disease in unique(sig_meta$Disease)) {
  sub_dir <- file.path(output_dir, disease)
  dir.create(sub_dir, showWarnings = FALSE)
  files_to_copy <- sig_meta$file[sig_meta$Disease == disease]
  matched <- all_files[file_names %in% files_to_copy]
  
  if (length(matched) > 0)
    file.copy(matched, sub_dir, overwrite = TRUE)
}
cat("Done! Created", length(unique(sig_meta$Disease)), "folders and copied matching signatures CSVs.\n")
#####


# ##### For manual annotation #####
# library(dplyr)
# library(stringr)
# library(purrr)
# library(writexl)
# 
# meta <- read.csv('output_files/signature_for_drug_discovery/sig_meta_final.csv', row.names = 1)
# # colnames(meta)
# meta <- meta %>% filter(signature_source == 'GEO')
# # table(meta$signature_source)
# # table(meta$case_source)
# # table(meta$control_source)
# # meta <- meta %>% filter(control_source == 'GEO')
# # table(meta$meta_source)
# # table(meta$OrganRegion)
# # table(meta$exp_type)
# # length(unique(meta$Disease))
# # head(meta, 2)
# # table(meta$control_source)
# 
# # length(list.files('data/exp_data/invivo/'))
# # example <- readRDS('data/exp_data/invivo/Acute-Cellular-Rejection_Kidney_Normal_GSE131179.rds')
# # names(example)
# # head(example$metadata, 2)
# # head(example$expr)[, 1:5]
# 
# meta <- meta %>% mutate(sig_file = str_replace(file, "\\.csv$", ".rds"), sig_type = if_else(control_source == "GEO", "matched", "unmatched"))
# # table(meta$control_source)
# # table(meta$sig_type)
# 
# matched_exp_dir <- "data/exp_data/invivo"
# unmatched_exp_dir <- "data/exp_data_unmatched/invivo"
# matched_files <- list.files(matched_exp_dir,   pattern = "\\.rds$")
# unmatched_files <- list.files(unmatched_exp_dir, pattern = "\\.rds$")
# 
# meta_matched <- meta %>% filter(sig_type == "matched", sig_file %in% matched_files)
# length(unique(meta_matched$sig_file))
# read_gsm_meta_by_filevec <- function(file_vec, dir) {
#   map_dfr(file_vec, function(f) {
#     obj <- readRDS(file.path(dir, f))
#     md  <- obj$metadata
#     md %>% mutate(sig_file = f)
#   })
# }
# 
# gsm_meta_matched <- read_gsm_meta_by_filevec(meta_matched$sig_file, matched_exp_dir) %>% mutate(sig_type = "matched")
# meta_unmatched <- meta %>% filter(sig_type == "unmatched") %>% mutate(rds_file = str_replace(sig_file, "_GTEx_", "_"))
# # missing_in_unmatched <- setdiff(meta_unmatched$rds_file, unmatched_files)
# # length(missing_in_unmatched); head(missing_in_unmatched)
# meta_unmatched <- meta_unmatched %>% filter(rds_file %in% unmatched_files)
# # length(unique(meta_unmatched$sig_file))
# read_gsm_meta_unmatched <- function(df, dir) {
#   map_dfr(seq_len(nrow(df)), function(i) {
#     f_sig <- df$sig_file[i] # signature ID used in `meta`
#     f_rds <- df$rds_file[i] # actual file name on disk
#     obj <- readRDS(file.path(dir, f_rds))
#     md  <- obj$metadata
#     md %>% mutate(sig_file = f_sig)
#   })
# }
# gsm_meta_unmatched <- read_gsm_meta_unmatched(meta_unmatched, unmatched_exp_dir) %>% mutate(sig_type = "unmatched", Control_type = "GTEX")
# 
# # colnames(gsm_meta_matched)
# # colnames(gsm_meta_unmatched)
# gsm_meta <- bind_rows(gsm_meta_matched, gsm_meta_unmatched)
# length(unique(gsm_meta$sig_file))
# colnames(gsm_meta)
# colnames(meta)
# 
# 
# cols_to_drop <- c("file", "Disease", "OrganRegion", "signature_source", "case_source", "control_source", "meta_source")
# meta_keep <- meta %>% select(-all_of(cols_to_drop))
# gsm_sig_meta_merged <- gsm_meta %>% left_join(meta_keep, by = c("sig_file", "sig_type"))
# # colnames(gsm_sig_meta_merged)
# # colSums(is.na(gsm_sig_meta_merged))
# # length(unique(gsm_sig_meta_merged$sig_file))
# # length(unique(gsm_sig_meta_merged$GSM_ID))
# # table(gsm_sig_meta_merged$sig_type)
# 
# val_data <- read.csv('output_files/signature_for_drug_discovery/validation/disignatlas/cor_results.csv', row.names = 1)
# # colnames(val_data)
# # head(val_data, 2)
# val_map <- val_data %>% mutate(sig_file = str_replace(file, "\\.csv$", ".rds")) %>%
#   select(sig_file, Cor_Pearson, Cor_Spearman, Cor_genes_n = n_genes) %>% distinct(sig_file, .keep_all = TRUE)
# # any(duplicated(val_map$sig_file))
# gsm_sig_meta_val_merged <- gsm_sig_meta_merged %>% left_join(val_map, by = "sig_file") %>% mutate(Validation = if_else(is.na(Cor_Pearson), "No", "Yes")) %>%
#   relocate(Validation, Cor_Pearson, Cor_Spearman, Cor_genes_n, .after = sig_file)
# head(gsm_sig_meta_val_merged, 2)
# 
# # table(gsm_sig_meta_val_merged$Validation)
# # y <- gsm_sig_meta_val_merged %>% filter(Validation == 'Yes')
# # n <- gsm_sig_meta_val_merged %>% filter(Validation == 'No')
# # length(unique(y$Disease))
# # length(unique(n$Disease))
# write_xlsx(gsm_sig_meta_val_merged, "output_files/gsm_used_for_sig.xlsx")
# #####