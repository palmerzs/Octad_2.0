################################################################################################
library(dplyr)
library(stringr)

## ---- 1. Category -> pathway keyword regex ------------------------
category_keywords <- c(
  VIR              = "INTERFERON|VIRAL|ANTIVIRAL|IFN|TOLL_LIKE_RECEPTOR|RIG_I|PATTERN_RECOGNITION|INNATE_IMMUNE|DEFENSE_RESPONSE_TO_VIRUS",
  BAC              = "TOLL_LIKE_RECEPTOR|INNATE_IMMUNE|INFLAMMATORY_RESPONSE|NEUTROPHIL|COMPLEMENT|NFKB|IL6|TNF|DEFENSE_RESPONSE_TO_BACTERIUM",
  PAR              = "INNATE_IMMUNE|EOSINOPHIL|TH2|IL4|IL5|IL13|COMPLEMENT|DEFENSE_RESPONSE",
  ACUTE_INJURY     = "INFLAMMATORY_RESPONSE|APOPTOSIS|HYPOXIA|COMPLEMENT|COAGULATION|TNF|IL6|NFKB|REACTIVE_OXYGEN",
  TRANSPLANT       = "ALLOGRAFT_REJECTION|T_CELL|INTERFERON_GAMMA|CYTOTOXIC|ANTIGEN_PROCESSING|MHC",
  AUTOIMMUNE       = "T_CELL|B_CELL|INTERFERON|CYTOKINE|COMPLEMENT|TNF|IL6|JAK_STAT|AUTOIMMUNE|ADAPTIVE_IMMUNE",
  ALLERGY          = "IL4|IL5|IL13|TH2|MAST_CELL|EOSINOPHIL|ALLERGIC|HISTAMINE",
  SKIN_BARRIER     = "KERATINOCYTE|EPIDERMIS|CORNIFIED_ENVELOPE|SKIN_DEVELOPMENT|KERATINIZATION|KERATIN",
  SKIN_INFLAM      = "IL17|IL23|TH17|KERATINOCYTE|INFLAMMATORY_RESPONSE|INTERFERON|EPIDERMIS",
  GI_IBD           = "TNF|IL6|IL23|TH17|NFKB|INFLAMMATORY_BOWEL|MUCOSAL|EPITHELIAL_BARRIER|DEFENSE_RESPONSE",
  GI_MOTILITY      = "MUSCLE_CONTRACTION|SMOOTH_MUSCLE|NEUROMUSCULAR|GASTRIC|ENTERIC_NERVOUS",
  LIVER            = "BILE_ACID|CHOLESTEROL|FATTY_ACID|LIPID_METABOLISM|XENOBIOTIC|FIBROSIS|STELLATE|BILIARY|COMPLEMENT",
  KIDNEY           = "KIDNEY|NEPHRON|APOPTOSIS|COMPLEMENT|INFLAMMATORY|FIBROSIS|TGF_BETA|OXIDATIVE_PHOSPHORYLATION|HYPOXIA",
  CARDIO           = "MUSCLE_CONTRACTION|CARDIAC|OXIDATIVE_PHOSPHORYLATION|HYPERTROPHY|EXTRACELLULAR_MATRIX|COAGULATION|COMPLEMENT|ANGIOGENESIS|VASCULAR",
  METABOLIC        = "OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|ADIPOGENESIS|INSULIN|GLYCOLYSIS|MTORC1|PPAR|LIPID_METABOLISM|INFLAMMATORY_RESPONSE",
  NEURODEGEN       = "NEURON|SYNAPT|OXIDATIVE_PHOSPHORYLATION|PROTEASOME|UBIQUITIN|APOPTOSIS|NEUROINFLAMMATION|MITOCHONDRI",
  PSYCH            = "SYNAPT|NEURON|GABA|GLUTAMATE|NEUROTRANSMITTER|AXON|DOPAMINE",
  NEURODEV         = "SYNAPT|NEURON_DIFFERENTIATION|AXON|CHROMATIN|TRANSCRIPTION|NERVOUS_SYSTEM_DEVELOPMENT",
  CANCER           = "CELL_CYCLE|E2F_TARGETS|MYC_TARGETS|G2M_CHECKPOINT|DNA_REPAIR|MITOTIC|APOPTOSIS|P53|PROLIFERATION|EPITHELIAL_MESENCHYMAL",
  REPRO            = "ESTROGEN|PROGESTERONE|HORMONE|ANGIOGENESIS|EXTRACELLULAR_MATRIX|EPITHELIAL",
  RESP_FIBROTIC    = "EXTRACELLULAR_MATRIX|COLLAGEN|TGF_BETA|EPITHELIAL_MESENCHYMAL_TRANSITION|WOUND_HEALING|FIBROSIS",
  SKELETAL         = "MUSCLE|MYOGENESIS|EXTRACELLULAR_MATRIX|COLLAGEN|OSSIFICATION|SARCOMERE",
  HEME             = "HEME|ERYTHROCYTE|OXIDATIVE_PHOSPHORYLATION|HYPOXIA|COMPLEMENT|HEMOGLOBIN",
  ORAL             = "INFLAMMATORY_RESPONSE|OSTEOCLAST|NFKB|NEUTROPHIL|BONE_RESORPTION",
  FIBROSIS_GENERAL = "EXTRACELLULAR_MATRIX|COLLAGEN|TGF_BETA|EPITHELIAL_MESENCHYMAL_TRANSITION|WOUND_HEALING"
)

## ---- 2. Disease -> category assignment ----------------------------

VIR <- c("Adenoviridae-Infections","Caliciviridae-Infections","Chickenpox","Common-Cold",
         "Coronavirus-Infections","Covid-19","Cytomegalovirus-Infections","Dengue","Dna-Virus-Infections",
         "Enterovirus-Infections","Hantavirus-Infections","Hemorrhagic-Fever-Ebola","Hepatitis-b","Hepatitis-c",
         "Hiv-Infections","Influenza-Human","Metapneumovirus","Respiratory-Syncytial-Virus-Infections",
         "Rna-Virus-Infections","Virus-Diseases","Zika-Virus-Infection")

BAC <- c("Bacterial-Infections","Campylobacter-Infections","Chlamydia-Infections","Gonorrhea",
         "Haemophilus-Infections","Pneumococcal-Infections","Pneumonia-Bacterial","Streptococcal-Infections",
         "Tuberculosis","Tuberculosis-Pulmonary","Vaginosis-Bacterial")

PAR <- c("Hookworm-Infections","Leishmaniasis","Malaria","Schistosomiasis","Scabies")

ACUTE_INJURY <- c("Acute-Kidney-Injury","Acute-Lung-Injury","Burns-Chemical","Chorioamnionitis",
                  "Fetal-Inflammatory-Response-Syndrome","Hypoxia-Ischemia-Brain","Macrophage-Activation-Syndrome",
                  "Reperfusion-Injury","Respiratory-Distress-Syndrome","Sepsis","Systemic-Inflammatory-Response-Syndrome")

TRANSPLANT <- c("Acute-Cellular-Rejection","Acute-t-Cell-Mediated-Rejection",
                "Bronchiolitis-Obliterans-Syndrome","Graft-vs-Host-Disease")

AUTOIMMUNE <- c("Antiphospholipid-Syndrome","Arthritis-Juvenile","Arthritis-Psoriatic",
                "Arthritis-Rheumatoid","Autoimmune-Lymphoproliferative-Syndrome","Castleman-Disease","Celiac-Disease",
                "Common-Variable-Immunodeficiency","Cryopyrin-Associated-Periodic-Syndromes","Dermatomyositis",
                "Familial-Mediterranean-Fever","Giant-Cell-Arteritis","Hepatitis-Autoimmune",
                "Hereditary-Autoinflammatory-Diseases","Lupus-Erythematosus-Systemic","Lupus-Nephritis","Polymyositis",
                "Rheumatic-Fever","Rheumatoid-Arthritis-Systemic-Juvenile","Sarcoidosis","Scleroderma-Localized",
                "Sj-Gren-Mikulicz-Syndrome","Sjogren-s-Syndrome","Spondylitis-Ankylosing","Still-s-Disease-Adult-Onset",
                "Vexas-Syndrome")

ALLERGY <- c("Asthma","Dermatitis-Allergic-Contact","Dermatitis-Atopic","Drug-Eruptions",
             "Eosinophilic-Esophagitis","Rhinitis-Allergic","Rhinitis-Allergic-Seasonal")

SKIN_BARRIER <- c("Erythroderma-Congenital-with-Palmoplantar-Keratoderma-Hypotrichosis-and-Hyper-Ige",
                  "Hyperkeratosis-Epidermolytic","Ichthyosiform-Erythroderma-Congenital","Ichthyosis-Lamellar",
                  "Netherton-Syndrome","Keratosis-Actinic")

SKIN_INFLAM <- c("Alopecia-Areata","Hidradenitis-Suppurativa","Keloid","Lichen-Planus",
                 "Necrobiosis-Lipoidica","Pityriasis-Rubra-Pilaris","Psoriasis","Pruritus")

GI_IBD <- c("Colitis-Lymphocytic","Colitis-Ulcerative","Crohn-Disease","Inflammatory-Bowel-Diseases",
            "Irritable-Bowel-Syndrome")

GI_MOTILITY <- c("Esophageal-Achalasia","Esophagitis","Gastroesophageal-Reflux","Gastroparesis",
                 "Heartburn","Non-Erosive-Reflux-Disease","Barrett-Esophagus","Hirschsprung-Disease")

LIVER <- c("Alagille-Syndrome","Biliary-Atresia","Budd-Chiari-Syndrome","Cholangitis-Sclerosing",
           "Choledochal-Cyst","Hepatitis","Hepatitis-Chronic","Idiopathic-Noncirrhotic-Portal-Hypertension",
           "Liver-Cirrhosis","Liver-Diseases","Non-Alcoholic-Fatty-Liver-Disease")

KIDNEY <- c("Diabetic-Nephropathies","Glomerulonephritis-Iga","Glomerulonephritis-Membranoproliferative",
            "Glomerulonephritis-Membranous","Glomerulosclerosis-Focal-Segmental","Kidney-Calculi","Nephrosis-Lipoid",
            "Uremia")

CARDIO <- c("Arrhythmogenic-Right-Ventricular-Dysplasia","Atherosclerosis","Cardiomyopathies",
            "Cardiomyopathy-Dilated","Cardiomyopathy-Hypertrophic","Cardiomyopathy-Restrictive",
            "Chronic-Limb-Threatening-Ischemia","Heart-Failure","Heart-Failure-Systolic",
            "Heart-Septal-Defects-Ventricular","Hypertension","Hypertension-Pulmonary",
            "Hypoplastic-Left-Heart-Syndrome","Intracranial-Aneurysm","Myocardial-Ischemia",
            "Peripheral-Arterial-Disease","Pulmonary-Arterial-Hypertension","Ventricular-Dysfunction-Right")

METABOLIC <- c("Acromegaly","Cachexia","Diabetes-Mellitus","Diabetes-Mellitus-Type-1",
               "Diabetes-Mellitus-Type-2","Glucose-Intolerance","Hypogonadism-and-Testicular-Atrophy",
               "Insulin-Resistance","Metabolically-Unhealthy-Obesity","Obesity","Obesity-Morbid",
               "Pituitary-Acth-Hypersecretion","Polycystic-Ovary-Syndrome","Prediabetic-State","Sarcopenia")

NEURODEGEN <- c("Alzheimer-Disease","Amyotrophic-Lateral-Sclerosis","Epilepsy-Temporal-Lobe",
                "Focal-Cortical-Dysplasia","Huntington-Disease","Lewy-Body-Disease","Multiple-Sclerosis",
                "Multiple-System-Atrophy","Myositis-Inclusion-Body","Parkinson-Disease","Pseudotumor-Cerebri")

PSYCH <- c("Bipolar-Disorder","Depressive-Disorder-Major","Opioid-Related-Disorders","Schizophrenia",
           "Substance-Related-Disorders")

NEURODEV <- c("Down-Syndrome","Fragile-x-Syndrome","Rett-Syndrome","Williams-Syndrome",
              "Myotonic-Dystrophy","Tuberous-Sclerosis")

CANCER <- c("Adenocarcinoma-of-Esophagus","Adenoma","Adenoma-Liver-Cell","Adenoma-Pleomorphic",
            "Adenomatous-Polyposis-Coli","Adenomatous-Polyps","Adrenal-Cortex-Cancer","Adrenocortical-Carcinoma",
            "Astrocytoma","Brain-Lower-Grade-Glioma","Brain-Neoplasms","Breast-Neoplasms","Buschke-Lowenstein-Tumor",
            "Carcinoid-Tumor","Carcinoma-Basal-Cell","Carcinoma-Ductal-Breast","Carcinoma-Hepatocellular",
            "Carcinoma-Large-Cell","Carcinoma-Neuroendocrine","Carcinoma-Non-Small-Cell-Lung",
            "Carcinoma-Pancreatic-Ductal","Carcinoma-Papillary","Carcinoma-Renal-Cell","Carcinoma-Signet-Ring-Cell",
            "Carcinoma-Squamous-Cell","Cholangiocarcinoma","Chondroblastoma","Choroid-Plexus-Carcinoma",
            "Colonic-Neoplasms","Colorectal-Neoplasms","Colorectal-Neoplasms-Hereditary-Nonpolyposis",
            "Conjunctival-Neoplasms","Cystadenocarcinoma-Serous","Desmoplastic-Small-Round-Cell-Tumor",
            "Diffuse-Intrinsic-Pontine-Glioma","Diffuse-Large-b-Cell-Lymphoma","Endometrial-Neoplasms",
            "Endometrial-Stromal-Tumors","Ependymoma","Esophageal-Neoplasms","Esophageal-Squamous-Cell-Carcinoma",
            "Fibrolamellar-Carcinoma","Fibrolamellar-Hepatocellular-Carcinoma","Fibroma","Fibrosarcoma",
            "Follicular-Neoplasm","Gallbladder-Neoplasms","Ganglioglioma","Ganglioneuroma",
            "Gastro-Enteropancreatic-Neuroendocrine-Tumor","Gastrointestinal-Stromal-Tumors","Glioblastoma","Glioma",
            "Head-and-Neck-Neoplasms","Hemangiosarcoma","Hepatoblastoma","Hodgkin-Disease","Hyperplasia",
            "Kidney-Chromophobe","Kidney-Neoplasms","Laryngeal-Neoplasms","Leiomyoma","Leiomyosarcoma",
            "Leukemia-Lymphocytic-Chronic-b-Cell","Leukemia-Megakaryoblastic-Acute","Leukemia-Myeloid-Acute",
            "Leukemia-Myelomonocytic-Chronic","Leukemia-Myelomonocytic-Juvenile","Lipoma","Liposarcoma",
            "Lung-Neoplasms","Lymphoma","Lymphoma-b-Cell-Marginal-Zone","Lymphoma-Follicular","Lymphoma-t-Cell",
            "Lymphoma-t-Cell-Peripheral","Medulloblastoma","Melanocytic-Nevus-Syndrome-Congenital","Melanoma",
            "Meningioma","Mesothelioma","Metaplasia","Monoclonal-Gammopathy-of-Undetermined-Significance",
            "Mouth-Neoplasms","Multiple-Myeloma","Mycosis-Fungoides","Myelodysplastic-Syndromes",
            "Myeloproliferative-Disorders","Myoepithelioma","Myofibromatosis","Nasopharyngeal-Carcinoma",
            "Neoplasms-Neuroepithelial","Neoplasms-Plasma-Cell","Neoplasms-Squamous-Cell","Neuroblastoma",
            "Neuroectodermal-Tumor-Melanotic","Neurofibroma","Neurofibrosarcoma","Nevus-Pigmented",
            "Oligodendroastrocytoma","Oligodendroglioma","Osteosarcoma","Ovarian-Neoplasms","Pancreatic-Adenoma",
            "Pancreatic-Neoplasms","Perivascular-Epithelioid-Cell-Neoplasms","Pheochromocytoma-Paraganglioma",
            "Pituitary-Neoplasms","Polyps","Intestinal-Polyposis","Precursor-Cell-Lymphoblastic-Leukemia-Lymphoma",
            "Prostatic-Neoplasms","Pseudolymphoma","Rectal-Neoplasms","Retinoblastoma","Rhabdoid-Tumor",
            "Rhabdomyosarcoma","Rhabdomyosarcoma-Alveolar","Rhabdomyosarcoma-Embryonal","Sarcoma",
            "Sarcoma-Alveolar-Soft-Part","Sarcoma-Ewing","Sarcoma-Kaposi","Sarcoma-Synovial",
            "Small-Cell-Lung-Carcinoma","Squamous-Cell-Carcinoma-of-Head-and-Neck","Stomach-Neoplasms","Thymoma",
            "Thyroid-Cancer-Papillary","Thyroid-Carcinoma-Anaplastic","Thyroid-Neoplasms","Tongue-Neoplasms",
            "Urinary-Bladder-Neoplasms","Uterine-Cervical-Dysplasia","Uterine-Cervical-Neoplasms","Wilms-Tumor")

REPRO <- c("Abortion-Habitual","Adenomyosis","Endometriosis","Epididymitis","Infertility",
           "Pelvic-Organ-Prolapse","Pre-Eclampsia","Pregnancy-Loss-Recurrent-Susceptibility-to-1",
           "Pregnancy-Tubal","Premature-Birth","Fetal-Growth-Retardation","Uterine-Cervicitis",
           "Prostatic-Hyperplasia")

RESP_FIBROTIC <- c("Idiopathic-Pulmonary-Fibrosis","Lung-Diseases-Interstitial",
                   "Pulmonary-Disease-Chronic-Obstructive")

SKELETAL <- c("Collagen-Vi-Related-Muscular-Dystrophy","Chronic-Recurrent-Multifocal-Osteomyelitis",
              "Fibrous-Dysplasia-of-Bone","Mitochondrial-Myopathies","Muscular-Dystrophy-Facioscapulohumeral")

HEME <- c("Anemia-Sickle-Cell")

ORAL <- c("Gingivitis","Periodontitis")

FIBROSIS_GENERAL <- c("Fibrosis","Implant-Capsular-Contracture","Pterygium")

## ---- 3. Build base category map ------------------------------------
category_map <- c(
  setNames(rep("VIR", length(VIR)), VIR),
  setNames(rep("BAC", length(BAC)), BAC),
  setNames(rep("PAR", length(PAR)), PAR),
  setNames(rep("ACUTE_INJURY", length(ACUTE_INJURY)), ACUTE_INJURY),
  setNames(rep("TRANSPLANT", length(TRANSPLANT)), TRANSPLANT),
  setNames(rep("AUTOIMMUNE", length(AUTOIMMUNE)), AUTOIMMUNE),
  setNames(rep("ALLERGY", length(ALLERGY)), ALLERGY),
  setNames(rep("SKIN_BARRIER", length(SKIN_BARRIER)), SKIN_BARRIER),
  setNames(rep("SKIN_INFLAM", length(SKIN_INFLAM)), SKIN_INFLAM),
  setNames(rep("GI_IBD", length(GI_IBD)), GI_IBD),
  setNames(rep("GI_MOTILITY", length(GI_MOTILITY)), GI_MOTILITY),
  setNames(rep("LIVER", length(LIVER)), LIVER),
  setNames(rep("KIDNEY", length(KIDNEY)), KIDNEY),
  setNames(rep("CARDIO", length(CARDIO)), CARDIO),
  setNames(rep("METABOLIC", length(METABOLIC)), METABOLIC),
  setNames(rep("NEURODEGEN", length(NEURODEGEN)), NEURODEGEN),
  setNames(rep("PSYCH", length(PSYCH)), PSYCH),
  setNames(rep("NEURODEV", length(NEURODEV)), NEURODEV),
  setNames(rep("CANCER", length(CANCER)), CANCER),
  setNames(rep("REPRO", length(REPRO)), REPRO),
  setNames(rep("RESP_FIBROTIC", length(RESP_FIBROTIC)), RESP_FIBROTIC),
  setNames(rep("SKELETAL", length(SKELETAL)), SKELETAL),
  setNames(rep("HEME", length(HEME)), HEME),
  setNames(rep("ORAL", length(ORAL)), ORAL),
  setNames(rep("FIBROSIS_GENERAL", length(FIBROSIS_GENERAL)), FIBROSIS_GENERAL)
)

## ---- 4. Two-category overrides (diseases spanning two mechanisms) --
overrides <- c(
  "Acute-Kidney-Injury"     = "ACUTE_INJURY+KIDNEY",
  "Lupus-Nephritis"         = "AUTOIMMUNE+KIDNEY",
  "Diabetic-Nephropathies"  = "METABOLIC+KIDNEY",
  "Diabetic-Neuropathies"   = "METABOLIC+NEURODEGEN",
  "Barrett-Esophagus"       = "GI_MOTILITY+CANCER",
  "Keratosis-Actinic"       = "SKIN_BARRIER+CANCER",
  "Scleroderma-Localized"   = "AUTOIMMUNE+FIBROSIS_GENERAL",
  "Hypoxia-Ischemia-Brain"  = "ACUTE_INJURY+NEURODEGEN",
  "Myositis-Inclusion-Body" = "NEURODEGEN+SKELETAL",
  "Prostatic-Hyperplasia"   = "REPRO+CANCER"
)
# Diabetic-Neuropathies wasn't in a base list above - add it via override only
category_map[names(overrides)] <- overrides

## ---- 5. Sanity check against your actual signature folder ----------
# Run this after setting signature_dir, to catch any disease present in
# your data but missing from category_map (or vice versa) so you can
# patch gaps before running Section 7.
#
signature_dir <- 'output_files/signature_for_drug_discovery/signature_may20/'
actual_diseases <- list.files(signature_dir)
missing_from_map   <- setdiff(actual_diseases, names(category_map))
extra_in_map       <- setdiff(names(category_map), actual_diseases)
cat(length(missing_from_map), "diseases in your data have no category assigned:\n")
print(missing_from_map)
cat(length(extra_in_map), "diseases in the map are not in your data (harmless, just unused):\n")
print(extra_in_map)

## ---- 6. Expand to per-disease keyword regex and build final table --
resolve_regex <- function(codes) {
  parts <- str_split(codes, "\\+")[[1]]
  paste(unique(category_keywords[parts]), collapse = "|")
}

expected_biology <- tibble(
  disease = names(category_map),
  category = unname(category_map)
) %>%
  rowwise() %>%
  mutate(expected_keywords = resolve_regex(category)) %>%
  ungroup()

cat("expected_biology built for", nrow(expected_biology), "diseases\n")

################################################################################################


################################################################################################
## ---- 0. Setup --------------------------------------------------
required_pkgs <- c("msigdbr", "fgsea", "data.table", "dplyr",
                   "purrr", "stringr", "tibble")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new_pkgs, update = FALSE, ask = FALSE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

## ---- 1. Paths ---------------------------------------------------
meta_signature_dir <- "../output_files/signature_for_drug_discovery/meta_signatures/"
out_dir <- "../output_files/pathway_enrichment/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta_files <- list.files(meta_signature_dir, pattern = "_meta\\.csv$", full.names = TRUE)
disease_names <- str_remove(basename(meta_files), "_meta\\.csv$")
names(meta_files) <- disease_names

stopifnot(length(meta_files) > 0)
cat("Found", length(meta_files), "disease meta-signature files\n")

## ---- 2. MSigDB gene sets ----------------------------------------
get_pathways <- function() {
  h        <- msigdbr(species = "Homo sapiens", category = "H")
  reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  gobp     <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
  
  msig_all <- bind_rows(
    mutate(h,        collection = "HALLMARK"),
    mutate(reactome, collection = "REACTOME"),
    mutate(gobp,     collection = "GOBP")
  )
  
  pathway_list   <- split(msig_all$gene_symbol, msig_all$gs_name)
  collection_map <- msig_all %>% distinct(gs_name, collection) %>% tibble::deframe()
  
  list(pathways = pathway_list, collection_map = collection_map)
}

msig <- get_pathways()
cat("Loaded", length(msig$pathways), "gene sets (Hallmark + Reactome + GO:BP)\n")

## ---- 3. Ranking function -----------------------------------------
# Signed rank statistic per gene: sign(log2FC) * -log10(p).
# LR itself is unsigned, so direction comes from log2FoldChange_meta.
# Collapse duplicate gene symbols by keeping the most extreme stat.

make_ranks <- function(df, min_p = 1e-300) {
  df <- df %>%
    filter(!is.na(gene_symbol), gene_symbol != "",
           !is.na(pval_meta), !is.na(log2FoldChange_meta)) %>%
    mutate(
      p_floor   = pmax(pval_meta, min_p),
      rank_stat = sign(log2FoldChange_meta) * -log10(p_floor)
    ) %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = abs(rank_stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  ranks <- setNames(df$rank_stat, df$gene_symbol)
  ranks <- ranks[is.finite(ranks)]
  sort(ranks, decreasing = TRUE)
}

## ---- 4. Run fgsea for one disease --------------------------------
run_one_disease <- function(disease, file, pathways, collection_map, min_genes = 5) {
  message("Running: ", disease)
  df <- tryCatch(read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) < 50) {
    warning(disease, ": file missing or too few rows, skipped")
    return(NULL)
  }
  
  ranks <- make_ranks(df)
  if (length(ranks) < 500) {
    warning(disease, ": too few ranked genes (", length(ranks), "), skipped")
    return(NULL)
  }
  
  res <- tryCatch(
    fgsea(pathways = pathways, stats = ranks, minSize = min_genes,
          maxSize = 500, eps = 0),
    error = function(e) { warning(disease, ": fgsea failed - ", e$message); NULL }
  )
  if (is.null(res)) return(NULL)
  
  res <- as.data.table(res)
  res[, disease := disease]
  res[, collection := collection_map[pathway]]
  res[, leadingEdge := sapply(leadingEdge, paste, collapse = "/")]
  setcolorder(res, c("disease", "collection", "pathway"))
  res[order(padj)]
}

## ---- 5. Loop over all diseases -----------------------------------
all_results <- purrr::imap(meta_files, function(file, disease) {
  run_one_disease(disease, file, msig$pathways, msig$collection_map)
})

all_results <- all_results[!sapply(all_results, is.null)]
cat(length(all_results), "of", length(meta_files), "diseases completed successfully\n")

results_dt <- rbindlist(all_results, use.names = TRUE, fill = TRUE)

saveRDS(all_results, file.path(out_dir, "fgsea_results_by_disease.rds"))
fwrite(results_dt, file.path(out_dir, "fgsea_results_all_diseases.csv"))

## ---- 6. Significant-pathway summary table -------------------------
sig_summary <- results_dt[padj < 0.05, .(
  n_sig_pathways = .N,
  n_sig_hallmark = sum(collection == "HALLMARK"),
  n_sig_reactome = sum(collection == "REACTOME"),
  n_sig_gobp     = sum(collection == "GOBP"),
  top_pathway    = pathway[which.min(padj)],
  top_padj       = min(padj)
), by = disease]

fwrite(sig_summary, file.path(out_dir, "significant_pathway_summary.csv"))
cat(nrow(sig_summary), "diseases have at least one pathway at padj < 0.05\n")

## ---- 7. Manual concordance check across all diseases ----------------
# Comprehensive category-based lookup (all 350 diseases), built by
# grouping diseases into shared-mechanism categories rather than
# curating 350 individual keyword sets. See expected_biology_lookup.R
# for the full category assignments and a sanity-check block that flags
# any disease in your signature_dir with no category assigned.
# HEURISTIC SCREENING TOOL - review assignments for diseases you plan
# to feature prominently before reporting concordance numbers.

source("scripts/validation_expected_biology_lookup.R")  # defines `expected_biology` (disease, category, expected_keywords)

check_concordance <- function(disease_name, results_dt, expected_biology, top_n = 30) {
  pattern <- expected_biology$expected_keywords[expected_biology$disease == disease_name]
  if (length(pattern) == 0) return(NULL)
  
  top_hits <- results_dt[disease == disease_name][order(padj)][seq_len(min(top_n, .N))]
  top_hits[, concordant := str_detect(pathway, regex(pattern, ignore_case = TRUE))]
  
  list(
    disease              = disease_name,
    n_concordant          = sum(top_hits$concordant, na.rm = TRUE),
    n_checked              = nrow(top_hits),
    concordant_pathways   = top_hits[concordant == TRUE, pathway],
    top_hits              = top_hits
  )
}

# Runs across every disease present in both your results and the lookup.
representative_diseases <- intersect(expected_biology$disease, unique(results_dt$disease))
cat(length(representative_diseases), "diseases have both fgsea results and a category assignment\n")

concordance_reports <- lapply(representative_diseases, check_concordance,
                              results_dt = results_dt,
                              expected_biology = expected_biology)
names(concordance_reports) <- representative_diseases

# Build one summary table instead of printing 350 blocks to console
concordance_summary <- rbindlist(lapply(representative_diseases, function(d) {
  r <- concordance_reports[[d]]
  if (is.null(r)) return(NULL)
  data.table(
    disease = d,
    category = expected_biology$category[expected_biology$disease == d],
    n_concordant = r$n_concordant,
    n_checked = r$n_checked,
    concordance_rate = round(r$n_concordant / r$n_checked, 2),
    example_concordant_pathway = if (length(r$concordant_pathways) > 0) r$concordant_pathways[1] else NA_character_
  )
}), fill = TRUE)

fwrite(concordance_summary, file.path(out_dir, "concordance_summary_all_diseases.csv"))
cat("Median concordance rate:", median(concordance_summary$concordance_rate, na.rm = TRUE), "\n")

# Spot-check a handful in the console
for (d in head(representative_diseases, 5)) {
  r <- concordance_reports[[d]]
  if (is.null(r)) next
  cat("\n==", d, "==\n")
  cat(r$n_concordant, "/", r$n_checked, "top pathways matched expected biology\n")
  print(head(r$concordant_pathways, 5))
}

# ============================================================
# Step 8: Visualize top enriched pathways by disease category
# ============================================================
# Run this after Sections 0-7 have produced:
#   - results_dt              (from fgsea_results_all_diseases.csv)
#   - expected_biology        (disease -> category, from expected_biology_lookup.R)
#   - sig_summary             (per-disease significant pathway counts)
#
# If starting a fresh R session, uncomment the reload block below.

library(ggplot2)
library(dplyr)
library(data.table)
library(forcats)
library(stringr)

# ---- reload if starting fresh ----------------------------------------
out_dir <- "../output_files/pathway_enrichment/"
results_dt <- fread(file.path(out_dir, "fgsea_results_all_diseases.csv"))
sig_summary <- fread(file.path(out_dir, "significant_pathway_summary.csv"))
source("validation_expected_biology_lookup.R")  # defines expected_biology

fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 8.1 Attach disease category to every pathway result -------------
results_cat <- merge(results_dt, expected_biology[, c("disease", "category")],
                     by = "disease", all.x = TRUE)

## ---- 8.2 Per-category, per-pathway summary ----------------------------
# For each category, how many diseases (within that category) is a given
# pathway significant in (padj < 0.05), and what's its average direction
# (NES) and magnitude?
n_diseases_per_category <- expected_biology %>% count(category, name = "n_diseases_total")

category_pathway_summary <- results_cat[padj < 0.05, .(
  n_diseases_sig = uniqueN(disease),
  mean_NES = mean(NES),
  mean_abs_NES = mean(abs(NES))
), by = .(category, pathway, collection)] %>%
  left_join(n_diseases_per_category, by = "category") %>%
  mutate(pct_diseases_sig = n_diseases_sig / n_diseases_total)

fwrite(category_pathway_summary, file.path(out_dir, "category_pathway_summary.csv"))

## ---- 8.3 Top N pathways per category (by breadth across diseases) ----
top_n_per_category <- 8

top_pathways <- category_pathway_summary %>%
  group_by(category) %>%
  slice_max(order_by = pct_diseases_sig, n = top_n_per_category, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pathway_short = str_trunc(str_replace_all(pathway, "_", " "), 45))

## ---- 8.4 Faceted bar chart: top pathways per category -----------------
# Split into chunks of ~6 categories per figure so facets stay legible
categories_all <- unique(top_pathways$category)
chunks <- split(categories_all, ceiling(seq_along(categories_all) / 6))

for (i in seq_along(chunks)) {
  plot_data <- top_pathways %>% filter(category %in% chunks[[i]])
  
  p <- ggplot(plot_data, aes(x = fct_reorder(pathway_short, pct_diseases_sig),
                             y = pct_diseases_sig, fill = mean_NES)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ category, scales = "free_y", ncol = 2) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
                         name = "Mean NES") +
    labs(x = NULL, y = "% of diseases in category significant (padj < 0.05)",
         title = paste0("Top enriched pathways by disease category (", i, "/", length(chunks), ")")) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave(file.path(fig_dir, paste0("top_pathways_by_category_", i, ".png")),
         p, width = 12, height = 8, dpi = 300)
}

## ---- 8.5 Category x collection heatmap: enrichment "burden" ----------
heat_data <- results_cat[padj < 0.05, .(n_sig = .N), by = .(category, collection)] %>%
  left_join(n_diseases_per_category, by = "category") %>%
  mutate(sig_per_disease = n_sig / n_diseases_total)

p_heat <- ggplot(heat_data, aes(x = collection, y = fct_rev(category), fill = sig_per_disease)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Sig. pathways\nper disease") +
  labs(x = NULL, y = NULL,
       title = "Enrichment burden by disease category and gene set collection") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(fig_dir, "category_collection_heatmap.png"), p_heat, width = 8, height = 10, dpi = 300)

## ---- 8.6 GSEA-style dot plot per category ------------------------------
# Classic "pathway x disease" bubble plot: color = NES (direction),
# size = significance. Caps the number of diseases shown per category so
# it stays readable (prioritizes diseases with the most signal).
plot_category_dotplot <- function(cat_name, top_pathways_n = 10, max_diseases = 15) {
  cat_diseases <- expected_biology$disease[expected_biology$category == cat_name]
  cat_diseases <- intersect(cat_diseases, unique(results_cat$disease))
  
  if (length(cat_diseases) > max_diseases) {
    keep <- sig_summary[disease %in% cat_diseases][order(-n_sig_pathways)][1:max_diseases, disease]
    cat_diseases <- keep
  }
  
  top_paths <- top_pathways %>% filter(category == cat_name) %>% pull(pathway)
  if (length(top_paths) == 0 || length(cat_diseases) == 0) {
    warning(cat_name, ": nothing to plot, skipped")
    return(invisible(NULL))
  }
  
  dot_data <- results_cat %>%
    filter(disease %in% cat_diseases, pathway %in% top_paths) %>%
    mutate(pathway_short = str_trunc(str_replace_all(pathway, "_", " "), 40))
  
  p <- ggplot(dot_data, aes(x = disease, y = fct_reorder(pathway_short, NES),
                            size = -log10(pmax(padj, 1e-10)), color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0) +
    scale_size_continuous(name = "-log10(padj)", range = c(1, 6)) +
    labs(x = NULL, y = NULL, title = paste0(cat_name, ": top pathways across diseases")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(fig_dir, paste0("dotplot_", cat_name, ".png")), p, width = 10, height = 7, dpi = 300)
  p
}

# Generate for a starter set of categories - add/remove as you like
for (cat_name in c("CANCER", "AUTOIMMUNE", "NEURODEGEN", "METABOLIC", "VIR", "KIDNEY")) {
  plot_category_dotplot(cat_name)
}

## ---- 8.8 Biologically-driven theme summary (replaces the raw-count heatmap) --
# The category x collection heatmap in 8.5 is dominated by two artifacts that
# have nothing to do with disease biology: (1) GO:BP has ~7,300 terms vs.
# Hallmark's 50, so raw counts always favor GO:BP regardless of biology, and
# (2) GO:BP is highly redundant - thousands of near-duplicate terms can
# represent a handful of distinct programs. Composite categories with only
# 1-2 diseases (e.g. SKIN_BARRIER+CANCER = Keratosis-Actinic alone) also
# produce dramatic-looking but statistically meaningless single-disease
# outliers. This section replaces "how many pathways" with "which canonical
# biological programs show up, how broadly, and in which direction" - a
# claim you can defend to a biologist.

theme_keywords <- c(
  "Interferon/antiviral"         = "INTERFERON|ANTIVIRAL|DEFENSE_RESPONSE_TO_VIRUS",
  "Innate immune/inflammation"   = "\\bTNF\\b|NFKB|INFLAMMATORY_RESPONSE|\\bIL6\\b|COMPLEMENT",
  "Adaptive immune (T/B cell)"   = "T_CELL|B_CELL|ADAPTIVE_IMMUNE|ANTIGEN_PROCESSING|ALLOGRAFT_REJECTION",
  "Th2/allergic"                 = "\\bIL4\\b|\\bIL5\\b|\\bIL13\\b|MAST_CELL|EOSINOPHIL",
  "Th17/IL23 axis"                = "\\bIL17\\b|\\bIL23\\b|TH17",
  "Cell cycle/proliferation"     = "CELL_CYCLE|E2F_TARGETS|G2M_CHECKPOINT|MITOTIC",
  "Growth signaling (MYC/mTOR)"  = "MYC_TARGETS|MTORC1|PI3K_AKT",
  "Apoptosis/cell death"         = "APOPTOSIS|CASPASE|PROGRAMMED_CELL_DEATH",
  "DNA damage/repair"            = "DNA_REPAIR|DNA_DAMAGE|\\bP53\\b",
  "EMT/ECM/fibrosis"             = "EPITHELIAL_MESENCHYMAL|EXTRACELLULAR_MATRIX|COLLAGEN|TGF_BETA",
  "Angiogenesis/vascular"        = "ANGIOGENESIS|VEGF|VASCULAR_DEVELOPMENT",
  "Coagulation/complement"       = "COAGULATION|COMPLEMENT|PLATELET",
  "Oxidative phosphorylation"    = "OXIDATIVE_PHOSPHORYLATION|ELECTRON_TRANSPORT|MITOCHONDRI",
  "Lipid/cholesterol metabolism" = "FATTY_ACID|LIPID_METABOLISM|CHOLESTEROL|ADIPOGENESIS",
  "Glycolysis/glucose handling"  = "GLYCOLYSIS|GLUCOSE|\\bINSULIN\\b",
  "Hypoxia response"             = "HYPOXIA",
  "Neuronal/synaptic"            = "\\bNEURON\\b|SYNAPT|\\bAXON\\b|NEUROTRANSMITTER",
  "Muscle contraction"           = "MUSCLE_CONTRACTION|MYOGENESIS|SARCOMERE",
  "Epithelial/keratinization"    = "KERATINOCYTE|EPIDERMIS|CORNIFIED|KERATINIZATION",
  "Hormone signaling"            = "ESTROGEN|ANDROGEN|PROGESTERONE|HORMONE",
  "Protein homeostasis"          = "UNFOLDED_PROTEIN|PROTEASOME|UBIQUITIN",
  "Xenobiotic/bile metabolism"   = "XENOBIOTIC|BILE_ACID"
)

# Collapse composite "A+B" categories to their primary label for this rollup,
# so a single disease (n=1) doesn't get displayed as if it were a category
# pattern. Per-disease results (concordance_summary etc.) still use the full
# composite category - this collapsing is only for the aggregate visualization.
expected_biology <- expected_biology %>%
  mutate(category_primary = str_split(category, "\\+", simplify = TRUE)[, 1])

n_diseases_per_primary <- expected_biology %>% count(category_primary, name = "n_total")

results_primary <- merge(results_dt,
                         expected_biology[, c("disease", "category_primary")],
                         by = "disease", all.x = TRUE)

theme_summary <- rbindlist(lapply(names(theme_keywords), function(theme_name) {
  pattern <- theme_keywords[[theme_name]]
  hits <- results_primary[padj < 0.05 & str_detect(pathway, regex(pattern, ignore_case = TRUE))]
  if (nrow(hits) == 0) return(NULL)
  
  out <- hits[, .(n_diseases_hit = uniqueN(disease), mean_NES = mean(NES)),
              by = category_primary]
  out[, theme := theme_name]
  out
}))

theme_summary <- theme_summary %>%
  left_join(n_diseases_per_primary, by = "category_primary") %>%
  mutate(pct_diseases_hit = n_diseases_hit / n_total)

fwrite(theme_summary, file.path(out_dir, "biological_theme_summary.csv"))

## ---- 8.9 Bubble matrix: category x biological theme --------------------
# size = breadth (% of diseases in category showing this theme)
# color = direction (mean NES: red = up-regulated, blue = down-regulated)
# Rows/columns ordered by clustering, so categories with similar biological
# programs land next to each other - itself a form of validation: if
# categories that are known to share mechanisms (e.g. AUTOIMMUNE and
# TRANSPLANT, both T-cell driven) cluster together, that's independent
# support that the signatures are capturing real, expected biology.
# theme_summary <- read.csv('../output_files/pathway_enrichment/biological_theme_summary.csv')
wide_mat <- theme_summary %>%
  select(category_primary, theme, pct_diseases_hit) %>%
  tidyr::pivot_wider(names_from = theme, values_from = pct_diseases_hit, values_fill = 0) %>%
  tibble::column_to_rownames("category_primary") %>%
  as.matrix()

row_order <- rownames(wide_mat)[hclust(dist(wide_mat))$order]
col_order <- colnames(wide_mat)[hclust(dist(t(wide_mat)))$order]

theme_summary <- theme_summary %>%
  mutate(category_primary = factor(category_primary, levels = row_order),
         theme = factor(theme, levels = col_order))

p_theme <- ggplot(theme_summary, aes(x = theme, y = category_primary, size = pct_diseases_hit, color = mean_NES)) + geom_point() +
  scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0, name = "Mean NES\n(direction)") +
  scale_size_continuous(name = "% diseases\nin category", range = c(0.8, 6), labels = scales::percent) +
  labs(x = NULL, y = NULL, title = "Biological programs enriched by disease category", subtitle = "Bubble size = breadth across diseases in category; color = direction of enrichment; rows/columns clustered by similarity") +
  theme_minimal(base_size = 18) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15), axis.title = element_text(size = 17),
                                        plot.title = element_text(size = 20, face = "bold"), plot.subtitle = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 13),
                                        panel.grid.major = element_line(color = "grey92"))
p_theme
ggsave(file.path(fig_dir, "biological_theme_bubble_matrix.png"), p_theme, width = 14, height = 10, dpi = 300)


# Using the `treemap` package instead of treemapify/ggplot2 - it has a much
# smaller dependency footprint (no svglite/systemfonts version chain to fight)
if (!requireNamespace("treemap", quietly = TRUE)) install.packages("treemap")
library(treemap)

category_totals <- results_cat[padj < 0.05, .(n_sig_pathways = .N), by = category] %>%
  left_join(n_diseases_per_category, by = "category") %>%
  mutate(sig_per_disease = n_sig_pathways / n_diseases_total)

png(file.path(fig_dir, "category_treemap.png"), width = 10, height = 7, units = "in", res = 300)
treemap(category_totals,
        index = "category",
        vSize = "n_sig_pathways",
        vColor = "sig_per_disease",
        type = "value",
        palette = "RdYlBu",
        title = "Total significant-pathway volume by disease category",
        fontsize.labels = 10,
        border.col = "white")
dev.off()

cat("Figures saved to:", fig_dir, "\n")
cat("Generated: top_pathways_by_category_*.png,",
    "category_collection_heatmap.png,",
    "dotplot_<CATEGORY>.png,",
    "category_treemap.png\n")

fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 8.1 Attach disease category to every pathway result -------------
results_cat <- merge(results_dt, expected_biology[, c("disease", "category")],
                     by = "disease", all.x = TRUE)

## ---- 8.2 Per-category, per-pathway summary ----------------------------
# For each category, how many diseases (within that category) is a given
# pathway significant in (padj < 0.05), and what's its average direction
# (NES) and magnitude?
n_diseases_per_category <- expected_biology %>% count(category, name = "n_diseases_total")

category_pathway_summary <- results_cat[padj < 0.05, .(
  n_diseases_sig = uniqueN(disease),
  mean_NES = mean(NES),
  mean_abs_NES = mean(abs(NES))
), by = .(category, pathway, collection)] %>%
  left_join(n_diseases_per_category, by = "category") %>%
  mutate(pct_diseases_sig = n_diseases_sig / n_diseases_total)

fwrite(category_pathway_summary, file.path(out_dir, "category_pathway_summary.csv"))

## ---- 8.3 Top N pathways per category (by breadth across diseases) ----
top_n_per_category <- 8

top_pathways <- category_pathway_summary %>%
  group_by(category) %>%
  slice_max(order_by = pct_diseases_sig, n = top_n_per_category, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pathway_short = str_trunc(str_replace_all(pathway, "_", " "), 45))

## ---- 8.4 Faceted bar chart: top pathways per category -----------------
# Split into chunks of ~6 categories per figure so facets stay legible
categories_all <- unique(top_pathways$category)
chunks <- split(categories_all, ceiling(seq_along(categories_all) / 6))

for (i in seq_along(chunks)) {
  plot_data <- top_pathways %>% filter(category %in% chunks[[i]])
  
  p <- ggplot(plot_data, aes(x = fct_reorder(pathway_short, pct_diseases_sig),
                             y = pct_diseases_sig, fill = mean_NES)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ category, scales = "free_y", ncol = 2) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0,
                         name = "Mean NES") +
    labs(x = NULL, y = "% of diseases in category significant (padj < 0.05)",
         title = paste0("Top enriched pathways by disease category (", i, "/", length(chunks), ")")) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"))
  
  ggsave(file.path(fig_dir, paste0("top_pathways_by_category_", i, ".png")),
         p, width = 12, height = 8, dpi = 300)
}

## ---- 8.5 Category x collection heatmap: enrichment "burden" ----------
heat_data <- results_cat[padj < 0.05, .(n_sig = .N), by = .(category, collection)] %>%
  left_join(n_diseases_per_category, by = "category") %>%
  mutate(sig_per_disease = n_sig / n_diseases_total)

p_heat <- ggplot(heat_data, aes(x = collection, y = fct_rev(category), fill = sig_per_disease)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Sig. pathways\nper disease") +
  labs(x = NULL, y = NULL,
       title = "Enrichment burden by disease category and gene set collection") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(fig_dir, "category_collection_heatmap.png"), p_heat, width = 8, height = 10, dpi = 300)

## ---- 8.6 GSEA-style dot plot per category ------------------------------
# Classic "pathway x disease" bubble plot: color = NES (direction),
# size = significance. Caps the number of diseases shown per category so
# it stays readable (prioritizes diseases with the most signal).
plot_category_dotplot <- function(cat_name, top_pathways_n = 10, max_diseases = 15) {
  cat_diseases <- expected_biology$disease[expected_biology$category == cat_name]
  cat_diseases <- intersect(cat_diseases, unique(results_cat$disease))
  
  if (length(cat_diseases) > max_diseases) {
    keep <- sig_summary[disease %in% cat_diseases][order(-n_sig_pathways)][1:max_diseases, disease]
    cat_diseases <- keep
  }
  
  top_paths <- top_pathways %>% filter(category == cat_name) %>% pull(pathway)
  if (length(top_paths) == 0 || length(cat_diseases) == 0) {
    warning(cat_name, ": nothing to plot, skipped")
    return(invisible(NULL))
  }
  
  dot_data <- results_cat %>%
    filter(disease %in% cat_diseases, pathway %in% top_paths) %>%
    mutate(pathway_short = str_trunc(str_replace_all(pathway, "_", " "), 40))
  
  p <- ggplot(dot_data, aes(x = disease, y = fct_reorder(pathway_short, NES),
                            size = -log10(pmax(padj, 1e-10)), color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B", midpoint = 0) +
    scale_size_continuous(name = "-log10(padj)", range = c(1, 6)) +
    labs(x = NULL, y = NULL, title = paste0(cat_name, ": top pathways across diseases")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(fig_dir, paste0("dotplot_", cat_name, ".png")), p, width = 10, height = 7, dpi = 300)
  p
}

# Generate for a starter set of categories - add/remove as you like
for (cat_name in c("CANCER", "AUTOIMMUNE", "NEURODEGEN", "METABOLIC", "VIR", "KIDNEY")) {
  plot_category_dotplot(cat_name)
}

## ---- 8.7 Treemap: how much signal does each category carry? -----------
if (!requireNamespace("treemapify", quietly = TRUE)) install.packages("treemapify")
library(treemapify)
install.packages("svglite")
library(treemapify)

category_totals <- results_cat[padj < 0.05, .(n_sig_pathways = .N), by = category] %>%
  left_join(n_diseases_per_category, by = "category") %>%
  mutate(sig_per_disease = n_sig_pathways / n_diseases_total)

p_treemap <- ggplot(category_totals, aes(area = n_sig_pathways, fill = sig_per_disease,
                                         label = category)) +
  geom_treemap() +
  geom_treemap_text(color = "white", place = "centre", grow = FALSE, reflow = TRUE, size = 11) +
  scale_fill_viridis_c(name = "Sig. pathways\nper disease") +
  labs(title = "Total significant-pathway volume by disease category")

ggsave(file.path(fig_dir, "category_treemap.png"), p_treemap, width = 10, height = 7, dpi = 300)

cat("Figures saved to:", fig_dir, "\n")
cat("Generated: top_pathways_by_category_*.png,",
    "category_collection_heatmap.png,",
    "dotplot_<CATEGORY>.png,",
    "category_treemap.png\n")
################################################################################################
