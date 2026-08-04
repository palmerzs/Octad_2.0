# ==============================================================================
# DiseaseSig Portal — improved app.R
# Changes vs original:
#   1. Lazy memoized signature reads (no startup pre-loading) — fast startup
#   2. Gene Search: genes on X-axis, diseases ranked on Y-axis (by signed log2FC)
#   3. Disease Similarity: own disease + signature selectors (no dependency on Explorer)
#   4. Pathway Enrichment: own disease + signature selectors + padj/logFC controls
#   5. GLP Explorer — forest plot of a target gene's log2FC across diseases
#   6. Volcano plot: standard red/blue colours, top-N gene labels
#   7. Disease name prettification helper
#   8. URL bookmarking & state sharing (enableBookmarking)
#   9. Drug repurposing reversal logic surfaced in tables/banners
# ==============================================================================

.libPaths(c("/home/ubuntu/shubham_data/octad2_data/R_libs_4.3", .libPaths()))

library(shiny)
library(dplyr)
library(data.table)
library(DT)
library(ggplot2)
library(plotly)
library(stringr)
library(tidyr)
library(bslib)
library(bsicons)
library(enrichR)

# octad + octad.db are required for Drug Repurposing tab
# Both must be installed in the lib path set at line 1
suppressPackageStartupMessages({
  .octad_available <- tryCatch({
    library(octad, quietly = TRUE)
    library(octad.db, quietly = TRUE)
    TRUE
  }, error = function(e) {
    message("octad not available: ", conditionMessage(e))
    FALSE
  })
})

# Unique compound count from octad.db's LINCS library (EH7270), for the Home
# page value box. Computed once at startup; degrades to "N/A" if octad.db or
# ExperimentHub isn't reachable (e.g. local dev without the full DB). Never
# prompts interactively — ExperimentHub's "create cache dir?" prompt is
# suppressed explicitly, since that would hang a non-interactive R session.
n_lincs_compounds <- tryCatch({
  if (!.octad_available) stop("octad not loaded")
  # Preemptively create the cache dir so ExperimentHub's interactive
  # "create directory?" prompt never triggers (the exact path shown in the
  # prompt matches tools::R_user_dir("ExperimentHub", "cache")).
  eh_cache <- tryCatch(tools::R_user_dir("ExperimentHub", "cache"), error = function(e) NULL)
  if (!is.null(eh_cache) && !dir.exists(eh_cache)) dir.create(eh_cache, recursive = TRUE, showWarnings = FALSE)
  suppressMessages(suppressWarnings(
    try(ExperimentHub::setExperimentHubOption("ASK", FALSE), silent = TRUE)
  ))
  obj <- tryCatch(octad.db::get_ExperimentHub_data("EH7270"), error = function(e) NULL)
  if (is.null(obj)) {
    obj <- tryCatch({
      e <- new.env(); utils::data(list = "lincs_sig_info", package = "octad.db", envir = e)
      get(ls(e)[1], envir = e)
    }, error = function(e) NULL)
  }
  name_col <- intersect(c("pert_iname", "Drug", "drug", "pert_name"), names(obj))[1]
  if (is.null(obj) || is.na(name_col)) "N/A" else format(length(unique(obj[[name_col]])), big.mark = ",")
}, error = function(e) "N/A")

# ── paths ──────────────────────────────────────────────────────────────────────
signature_dir <- "/home/ubuntu/shubham_data/octad2_data/signature_may20/"
meta_signature_dir <- "/home/ubuntu/shubham_data/octad2_data/meta_signatures/"
example_drug_dir <- "/home/ubuntu/shubham_data/octad2_data/example_drug_signatures/"

# Example drug/perturbation signatures (Semaglutide/GLP-1, Chen Lab), listed
# dynamically so any file placed in example_drug_dir is automatically offered
# — no hardcoded filename or tissue.
example_drug_files <- tryCatch({
  f <- list.files(example_drug_dir, pattern = "\\.csv$", full.names = TRUE)
  setNames(f, tools::file_path_sans_ext(basename(f)))
}, error = function(e) character(0))

# Loads a precomputed <disease>_meta.csv if it exists; NULL if not found yet.
# This lets the app work before the meta-analysis has been generated/uploaded
# for every disease, degrading gracefully to single-dataset mode.
read_meta_signature <- function(disease_val) {
  path <- file.path(meta_signature_dir, paste0(disease_val, "_meta.csv"))
  if (!file.exists(path)) return(NULL)
  fread(path)
}

# ── helpers ────────────────────────────────────────────────────────────────────
pretty_disease <- function(x) {
  # "Carcinoma-Non-Small-Cell-Lung" → "Carcinoma Non Small Cell Lung"
  gsub("-", " ", x)
}

# Best-effort tissue/organ extraction from filename convention:
# <Disease>_<Tissue...>_<Source>_<GSExxxxx>.csv
# e.g. "Alzheimer-Disease_Brain_Basal_ganglia_Normal_GSE193438.csv" -> "Brain_Basal_ganglia"
#      "Abortion-Habitual_Endometrium_GTEx_GSE65099.csv"            -> "Endometrium"
# Falls back to "Unknown" if the convention doesn't match (heuristic, not guaranteed for every file).
parse_tissue <- function(disease, file) {
  base <- tools::file_path_sans_ext(file)
  esc  <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", disease)
  esc  <- gsub("[-_]", "[-_]", esc)  # tolerate hyphen/underscore swaps between folder & filename
  rest <- sub(paste0("^", esc, "_?"), "", base, ignore.case = TRUE)
  parts <- strsplit(rest, "_")[[1]]
  parts <- parts[parts != ""]
  if (length(parts) >= 3 && grepl("^GSE[0-9]+$", parts[length(parts)], ignore.case = TRUE)) {
    tissue_parts <- parts[seq_len(length(parts) - 2)]
    tissue <- paste(tissue_parts, collapse = "_")
  } else if (length(parts) == 2 && grepl("^GSE[0-9]+$", parts[length(parts)], ignore.case = TRUE)) {
    tissue <- parts[1]
  } else if (length(parts) >= 1) {
    tissue <- parts[1]
  } else {
    tissue <- "Unknown"
  }
  if (tissue == "") tissue <- "Unknown"
  tissue
}

build_signature_index <- function(dir) {
  files <- list.files(dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  disease_v <- basename(dirname(files))
  file_v    <- basename(files)
  data.frame(
    path         = files,
    disease      = disease_v,
    file         = file_v,
    signature_id = paste0(disease_v, " | ", file_v),
    tissue       = mapply(parse_tissue, disease_v, file_v, USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}

# ── startup: build index ───────────────────────────────────────────────────────
message("[startup] Building signature index …")
sig_index <- build_signature_index(signature_dir)
disease_choices  <- sort(unique(sig_index$disease))
disease_choices_named <- setNames(disease_choices, pretty_disease(disease_choices))
# Some datasets come from non-GEO sources (repository acronyms, not tissue
# names) that the filename-based tissue parser can't distinguish from a real
# organ. Excluded here so they don't show as misleading "tissue" options —
# those datasets remain included when "All tissues" is selected.
non_tissue_labels <- c("CBTN", "EGA", "ICGC", "SJC", "SRA", "TARGET", "UA")
tissue_choices <- sort(setdiff(unique(sig_index$tissue), non_tissue_labels))

# ── lazy, memoized signature reader ────────────────────────────────────────────
# No startup pre-loading: startup stays instant. Each CSV is read from disk only
# when first requested, then kept in a small in-session cache so repeated access
# (e.g. re-running a tab) is fast without holding all 979 files in memory upfront.
.sig_cache <- new.env(parent = emptyenv())

read_signature <- function(path) {
  if (!is.null(.sig_cache[[path]])) return(.sig_cache[[path]])
  dt <- fread(path)
  dt[, `:=`(
    gene_symbol      = fifelse(is.na(gene_symbol) | gene_symbol == "", identifier, gene_symbol),
    gene_description = fifelse(is.na(gene_description), "", gene_description)
  )]
  df <- as.data.frame(dt)
  .sig_cache[[path]] <- df
  df
}

# Pull one signature by disease + file (returns data.frame with disease/file cols)
get_sig <- function(disease_val, file_val) {
  path <- sig_index$path[sig_index$disease == disease_val & sig_index$file == file_val][1]
  read_signature(path)
}

# ── UI ─────────────────────────────────────────────────────────────────────────
ui <- function(request) navbarPage(
  title = div(
    style = "display:flex; align-items:center; gap:10px;",
    bs_icon("activity"), "OCTAD 2.0"
  ),
  windowTitle = "OCTAD 2.0",
  theme = bs_theme(
    version   = 5,
    bootswatch = "flatly",
    primary   = "#2C6EAB",
    secondary = "#5DADE2"
  ),
  id = "main_nav",
  
  # ── Home ────────────────────────────────────────────────────────────────────
  tabPanel("Home",
           div(class = "container-fluid p-0",
               div(style = "background: #2C6EAB; color: #fff; padding: 3rem 1rem;",
                   div(class = "container",
                       h1(style = "font-weight: 600; margin-bottom: 0.5rem;", bs_icon("activity"), " OCTAD 2.0"),
                       p(style = "font-size: 1.15rem; max-width: 900px; opacity: 0.95; margin-bottom: 0;",
                         "Transcriptomics-based drug discovery: disease signature exploration, ",
                         strong("disease-centric drug repurposing"), ", and ",
                         strong("drug-centric indication expansion"), "."
                       )
                   )
               ),
               div(class = "container mt-4",
                   layout_columns(
                     fill = FALSE,
                     value_box("Signatures", nrow(sig_index),                   showcase = bs_icon("files"),           theme = "info"),
                     value_box("Diseases",   length(unique(sig_index$disease)), showcase = bs_icon("folder2-open"),   theme = "primary"),
                     value_box("Drugs",      n_lincs_compounds,                 showcase = bs_icon("capsule"),         theme = "success")
                   ),
                   br(),
                   h4(style = "font-weight: 600; color: #2C6EAB;", bs_icon("heart-pulse"), " Disease Explorer"),
                   layout_columns(
                     fill = FALSE,
                     col_widths = c(3, 3, 3, 3),
                     card(style = "border-left: 4px solid #2C6EAB;",
                          card_body(
                            h6(bs_icon("bar-chart-line"), " Signature Explorer", style = "font-weight: 600;"),
                            p(class = "text-muted small mb-0", "Volcano plot, top DEGs, and signature summary for any selected disease.")
                          )),
                     card(style = "border-left: 4px solid #2C6EAB;",
                          card_body(
                            h6(bs_icon("search"), " Gene Explorer", style = "font-weight: 600;"),
                            p(class = "text-muted small mb-0", "Search genes across diseases, or trace a target gene's log2FC across every disease.")
                          )),
                     card(style = "border-left: 4px solid #2C6EAB;",
                          card_body(
                            h6(bs_icon("diagram-2"), " Disease Similarity", style = "font-weight: 600;"),
                            p(class = "text-muted small mb-0", "Correlate any disease signature against all others. Identify reversal candidates.")
                          )),
                     card(style = "border-left: 4px solid #2C6EAB;",
                          card_body(
                            h6(bs_icon("signpost-split"), " Pathway Enrichment", style = "font-weight: 600;"),
                            p(class = "text-muted small mb-0", "GO, KEGG, Reactome, WikiPathways, and Hallmark enrichment for any disease/direction.")
                          ))
                   ),
                   br(),
                   h4(style = "font-weight: 600; color: #196F3D;", bs_icon("capsule"), " Drug Explorer"),
                   layout_columns(
                     fill = FALSE,
                     col_widths = c(6, 6),
                     card(style = "border-left: 4px solid #196F3D;",
                          card_body(
                            h6(bs_icon("capsule"), " Drug Repurposing", style = "font-weight: 600;"),
                            p(class = "text-muted small mb-0", "Provide a disease signature, identify existing drugs that reverse its expression.")
                          )),
                     card(style = "border-left: 4px solid #196F3D;",
                          card_body(
                            h6(bs_icon("bullseye"), " Indication Expansion", style = "font-weight: 600;"),
                            p(class = "text-muted small mb-0", "Provide a drug signature, identify all potential indications.")
                          ))
                   ),
                   br(),
                   layout_columns(
                     col_widths = c(8, 4),
                     card(
                       card_header(style = "font-weight: 600;", bs_icon("info-circle"), " Drug repurposing logic"),
                       card_body(
                         tags$ul(class = "mb-0",
                                 tags$li("Disease up-regulated genes = activated programs."),
                                 tags$li("Disease down-regulated genes = suppressed protective programs."),
                                 tags$li("A drug that ", strong("reverses"), " the disease signature (negative correlation) is a candidate repurposing compound."),
                                 tags$li("Use Disease Similarity to find diseases with strongly negative correlations as proxy drug contexts.")
                         )
                       )
                     ),
                     card(
                       card_header(style = "font-weight: 600;", bs_icon("envelope"), " Contact"),
                       card_body(
                         p(class = "mb-0", tags$a(href = "mailto:contact@octad.org", "contact@octad.org"))
                       )
                     )
                   ),
                   br()
               )
           )
  ),
  
  # ── Disease Explorer menu ───────────────────────────────────────────────────
  navbarMenu("Disease Explorer", icon = bs_icon("heart-pulse"),
             
             # ── Disease Explorer ────────────────────────────────────────────────────────
             tabPanel("Signature Explorer",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     h5(bs_icon("funnel"), " Select Disease"),
                                     selectizeInput("disease", "Disease",
                                                    choices = c("Choose a disease…" = "", disease_choices_named),
                                                    selected = "",
                                                    options = list(placeholder = "Search disease…")),
                                     uiOutput("signature_selector"),
                                     hr(),
                                     h5(bs_icon("sliders"), " Thresholds"),
                                     numericInput("padj_cutoff", "Adjusted p-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
                                     numericInput("logfc_cutoff", "|log2FC| cutoff", 1, min = 0, step = 0.1),
                                     numericInput("top_n", "Top genes to show", 50, min = 5, step = 5),
                                     numericInput("label_n", "Label top N genes on volcano", 10, min = 0, step = 5),
                                     hr(),
                                     h5(bs_icon("search"), " Filter table"),
                                     textInput("gene_search", "Search gene / Ensembl ID / description", ""),
                                     downloadButton("download_filtered", "Download filtered genes", class = "btn-sm btn-outline-primary w-100")
                        ),
                        mainPanel(width = 9,
                                  br(),
                                  layout_columns(
                                    fill = FALSE,
                                    value_box("Total genes",       textOutput("total_genes"), showcase = bs_icon("diagram-3"),       theme = "secondary"),
                                    value_box("Significant genes", textOutput("sig_genes"),   showcase = bs_icon("activity"),        theme = "primary"),
                                    value_box("Up-regulated",      textOutput("up_genes"),    showcase = bs_icon("arrow-up-circle"), theme = "danger"),
                                    value_box("Down-regulated",    textOutput("down_genes"),  showcase = bs_icon("arrow-down-circle"), theme = "info")
                                  ),
                                  br(),
                                  navset_tab(
                                    nav_panel("Volcano Plot",
                                              br(),
                                              plotlyOutput("volcano_plot", height = "650px")
                                    ),
                                    nav_panel("Top DEGs",
                                              br(),
                                              h5(class = "text-danger", bs_icon("arrow-up-circle"), " Top up-regulated genes"),
                                              DTOutput("top_up"),
                                              br(),
                                              h5(class = "text-primary", bs_icon("arrow-down-circle"), " Top down-regulated genes"),
                                              DTOutput("top_down")
                                    ),
                                    nav_panel("Summary",
                                              br(),
                                              DTOutput("summary_table"),
                                              br(),
                                              plotlyOutput("direction_bar", height = "300px")
                                    )
                                  )
                        )
                      )
             ),
             
             # ── Gene Search ─────────────────────────────────────────────────────────────
             tabPanel("Gene Explorer",
                      navset_tab(
                        nav_panel("Search across diseases",
                                  sidebarLayout(
                                    sidebarPanel(width = 3,
                                                 h5(bs_icon("search"), " Gene Query"),
                                                 textInput("global_gene_search", "Gene symbol(s) or Ensembl ID(s)", ""),
                                                 helpText("Partial match supported. Comma-separated, e.g. GLP, TP53, EGFR"),
                                                 hr(),
                                                 h5(bs_icon("diagram-3"), " View mode"),
                                                 radioButtons("gene_search_mode", NULL,
                                                              choices = c("Across diseases" = "disease",
                                                                          "Across datasets"  = "dataset"),
                                                              selected = "disease"),
                                                 hr(),
                                                 h5(bs_icon("sliders"), " Filters"),
                                                 checkboxInput("gene_sig_only", "Significant hits only", TRUE),
                                                 numericInput("global_padj_cutoff",  "Adjusted p-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
                                                 numericInput("global_logfc_cutoff", "|log2FC| cutoff", 0.5, min = 0, step = 0.1),
                                                 numericInput("gene_heatmap_top_n",  "Max rows to display (ranked)", 75, min = 10, step = 5),
                                                 selectInput("gene_rank_by", "Rank by",
                                                             choices = c("log2FC (most up-regulated)"   = "desc_fc",
                                                                         "log2FC (most down-regulated)" = "asc_fc",
                                                                         "Significance (min padj)"      = "min_padj")),
                                                 hr(),
                                                 actionButton("run_gene_search", "Generate heatmap", class = "btn-primary w-100"),
                                                 br(), br(),
                                                 downloadButton("download_global_gene_search", "Download results", class = "btn-sm btn-outline-primary w-100")
                                    ),
                                    mainPanel(width = 9,
                                              br(),
                                              uiOutput("gene_search_banner"),
                                              br(),
                                              h5("log2FC heatmap — genes (X) × ranked diseases/datasets (Y)"),
                                              helpText("Red = up-regulated, Blue = down-regulated. In 'Across diseases' mode, multiple datasets per disease are averaged."),
                                              plotlyOutput("global_gene_heatmap", height = "1100px", width = "100%")
                                    )
                                  )
                        ),
                        nav_panel("Forest plot (single gene)",
                                  sidebarLayout(
                                    sidebarPanel(width = 3,
                                                 h5(bs_icon("bullseye"), " Target Gene"),
                                                 textInput("target_gene", "Gene symbol", ""),
                                                 helpText("Enter any gene symbol to explore its expression across all diseases."),
                                                 hr(),
                                                 h5(bs_icon("sliders"), " Forest plot options"),
                                                 numericInput("target_padj_cutoff",  "Highlight padj ≤", 0.05, min = 0, max = 1, step = 0.01),
                                                 numericInput("target_logfc_cutoff", "Highlight |log2FC| ≥", 0.5, min = 0, step = 0.1),
                                                 checkboxInput("target_sig_only",   "Show only significant entries", FALSE),
                                                 selectInput("target_sort", "Sort diseases by",
                                                             choices = c("log2FC (descending)" = "desc_fc",
                                                                         "log2FC (ascending)"  = "asc_fc",
                                                                         "Alphabetical"        = "alpha")),
                                                 numericInput("target_top_n", "Max diseases to display", 60, min = 10, step = 10),
                                                 hr(),
                                                 actionButton("run_target", "Generate forest plot", class = "btn-primary w-100"),
                                                 br(), br(),
                                                 downloadButton("download_target", "Download data", class = "btn-sm btn-outline-primary w-100")
                                    ),
                                    mainPanel(width = 9,
                                              br(),
                                              uiOutput("target_banner"),
                                              br(),
                                              plotlyOutput("target_forest", height = "900px", width = "100%"),
                                              br(),
                                              h5("Data table"),
                                              DTOutput("target_table")
                                    )
                                  )
                        )
                      )
             ),
             
             # ── Disease Similarity ──────────────────────────────────────────────────────
             tabPanel("Disease Similarity",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     h5(bs_icon("diagram-2"), " Query Disease"),
                                     selectizeInput("corr_disease", "Disease",
                                                    choices = c("Choose a disease…" = "", disease_choices_named),
                                                    selected = "",
                                                    options = list(placeholder = "Search disease…")),
                                     uiOutput("corr_signature_selector"),
                                     hr(),
                                     h5(bs_icon("sliders"), " Correlation settings"),
                                     selectInput("corr_method", "Method", choices = c("Spearman" = "spearman", "Pearson" = "pearson")),
                                     numericInput("min_shared_genes", "Min shared genes", value = 500, min = 100, step = 500),
                                     numericInput("corr_top_n", "Top N to display (each direction)", value = 50, min = 5, step = 5),
                                     hr(),
                                     actionButton("run_corr", "Run correlation", class = "btn-primary w-100"),
                                     br(), br(),
                                     downloadButton("download_corr", "Download results", class = "btn-sm btn-outline-primary w-100")
                        ),
                        mainPanel(width = 9,
                                  br(),
                                  uiOutput("corr_query_banner"),
                                  br(),
                                  navset_tab(
                                    nav_panel("Bar chart",
                                              br(),
                                              plotlyOutput("corr_plot", height = "700px")
                                    ),
                                    nav_panel("Heatmap",
                                              br(),
                                              plotlyOutput("corr_heatmap", height = "1100px", width = "100%")
                                    ),
                                    nav_panel("Table",
                                              br(),
                                              DTOutput("corr_table")
                                    )
                                  )
                        )
                      )
             ),
             
             # ── Pathway Enrichment ──────────────────────────────────────────────────────
             tabPanel("Pathway Enrichment",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     h5(bs_icon("signpost-split"), " Disease"),
                                     selectizeInput("enrich_disease", "Disease",
                                                    choices = c("Choose a disease…" = "", disease_choices_named),
                                                    selected = "",
                                                    options = list(placeholder = "Search disease…")),
                                     radioButtons("enrich_mode", NULL,
                                                  choices = c("Across disease" = "disease",
                                                              "Select dataset"                                   = "dataset"),
                                                  selected = "disease"),
                                     conditionalPanel("input.enrich_mode == 'dataset'",
                                                      uiOutput("enrich_signature_selector")
                                     ),
                                     hr(),
                                     h5(bs_icon("sliders"), " Parameters"),
                                     numericInput("enrich_padj_cutoff",  "Adjusted p-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
                                     numericInput("enrich_logfc_cutoff", "|log2FC| cutoff", 1, min = 0, step = 0.1),
                                     selectInput("enrich_direction", "Gene set direction",
                                                 choices = c("Up-regulated" = "up", "Down-regulated" = "down")),
                                     selectInput("enrich_db", "Database",
                                                 choices = c(
                                                   "GO Biological Process 2023"  = "GO_Biological_Process_2023",
                                                   "KEGG 2021 Human"             = "KEGG_2021_Human",
                                                   "Reactome 2022"               = "Reactome_2022",
                                                   "WikiPathway 2023 Human"      = "WikiPathway_2023_Human",
                                                   "MSigDB Hallmark 2020"        = "MSigDB_Hallmark_2020"
                                                 )),
                                     numericInput("enrich_top_n", "Top pathways to plot", 20, min = 5, step = 5),
                                     hr(),
                                     actionButton("run_enrich", "Run enrichment", class = "btn-primary w-100"),
                                     br(), br(),
                                     downloadButton("download_enrich", "Download results", class = "btn-sm btn-outline-primary w-100")
                        ),
                        mainPanel(width = 9,
                                  br(),
                                  uiOutput("enrich_query_banner"),
                                  br(),
                                  navset_tab(
                                    nav_panel("Bar chart",    br(), plotlyOutput("enrich_plot",  height = "650px")),
                                    nav_panel("Result table", br(), DTOutput("enrich_table"))
                                  )
                        )
                      )
             )   # end Pathway Enrichment tabPanel
             
  ),  # end Disease Explorer navbarMenu
  
  # ── Drug Explorer menu ──────────────────────────────────────────────────────
  navbarMenu("Drug Explorer", icon = bs_icon("capsule"),
             
             # ── Drug Repurposing (RGES) ─────────────────────────────────────────────────
             tabPanel("Drug Repurposing",
                      sidebarLayout(
                        sidebarPanel(width = 3,
                                     h5(bs_icon("capsule"), " Signature source"),
                                     radioButtons("rges_source", NULL,
                                                  choices = c("Portal disease" = "portal", "Upload my own signature" = "upload"),
                                                  selected = "portal"),
                                     
                                     # Portal disease selectors
                                     conditionalPanel("input.rges_source == 'portal'",
                                                      selectizeInput("rges_disease", "Disease",
                                                                     choices = c("Choose a disease…" = "", disease_choices_named),
                                                                     selected = "",
                                                                     options = list(placeholder = "Search disease…")),
                                                      radioButtons("rges_disease_mode", NULL,
                                                                   choices = c("Meta-analysis signature (recommended)" = "meta",
                                                                               "Select one dataset"                    = "single"),
                                                                   selected = "meta"),
                                                      conditionalPanel("input.rges_disease_mode == 'single'",
                                                                       uiOutput("rges_signature_selector_single")
                                                      ),
                                                      uiOutput("rges_signature_selector")
                                     ),
                                     # Upload path (feature 7)
                                     conditionalPanel("input.rges_source == 'upload'",
                                                      fileInput("rges_upload", "DE table (CSV)", accept = ".csv"),
                                                      helpText("Must contain a gene-symbol column and a log2FoldChange column (padj optional)."),
                                                      uiOutput("rges_col_mapper")
                                     ),
                                     hr(),
                                     
                                     h5(bs_icon("sliders"), " Signature thresholds"),
                                     numericInput("rges_padj_cutoff",  "Adjusted p-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
                                     numericInput("rges_logfc_cutoff", "|log2FC| cutoff", 0.5, min = 0, step = 0.1),
                                     numericInput("rges_max_genes",    "Max genes per direction", 100, min = 10, step = 10),
                                     conditionalPanel("input.rges_disease_mode == 'meta' && input.rges_source == 'portal'",
                                                      hr(),
                                                      h5(bs_icon("bar-chart-steps"), " Meta-analysis filters"),
                                                      numericInput("rges_max_i2",     "Max heterogeneity (I², %)", 100, min = 0, max = 100, step = 5),
                                                      helpText("Lower I² excludes genes where datasets strongly disagree."),
                                                      numericInput("rges_min_k",      "Min datasets per gene", 1, min = 1, step = 1),
                                                      helpText("Require a gene to be measured in at least this many datasets.")
                                     ),
                                     hr(),
                                     
                                     # Cell-line context (feature 4)
                                     h5(bs_icon("diagram-3"), " Reference context"),
                                     selectizeInput("rges_cells", "LINCS cell lines (optional)",
                                                    choices = NULL, multiple = TRUE,
                                                    options = list(placeholder = "Default: all cell lines")),
                                     helpText("Restrict scoring to disease-relevant cell lines, or leave blank for all."),
                                     numericInput("rges_permutations", "Permutations", 10000, min = 1000, step = 1000),
                                     hr(),
                                     
                                     div(class = "text-muted small mb-2", bs_icon("clock"), " Typically takes 1–3 minutes, depending on signature size and gene count."),
                                     actionButton("run_rges", "Run RGES scoring", class = "btn-primary w-100"),
                                     br(), br(),
                                     downloadButton("download_rges", "Download drug candidates", class = "btn-sm btn-outline-primary w-100")
                        ),
                        mainPanel(width = 9,
                                  br(),
                                  uiOutput("rges_banner"),
                                  br(),
                                  navset_card_tab(
                                    # 1 + 2: ranked candidates + annotation
                                    nav_panel("Drug Candidates",
                                              br(),
                                              layout_columns(fill = FALSE,
                                                             value_box("Compounds scored", textOutput("rges_n_drugs"),   showcase = bs_icon("capsule"),        theme = "primary"),
                                                             value_box("Strong reversers", textOutput("rges_n_strong"),  showcase = bs_icon("arrow-down-circle"), theme = "success"),
                                                             value_box("FDA-approved hits", textOutput("rges_n_fda"),    showcase = bs_icon("patch-check"),     theme = "info")
                                              ),
                                              br(),
                                              checkboxInput("rges_fda_only", "Show FDA-approved / launched only", FALSE),
                                              helpText("sRGES ranges roughly −1 to 1. More negative = stronger reversal of the disease signature = better repurposing candidate."),
                                              DTOutput("rges_table")
                                    ),
                                    # 3: reversal plot
                                    nav_panel("Reversal Plot",
                                              br(),
                                              uiOutput("rges_drug_picker"),
                                              plotlyOutput("rges_reversal_plot", height = "600px")
                                    ),
                                    # 6: approved / clinical drugs among hits
                                    nav_panel("Approved & Clinical",
                                              br(),
                                              helpText("FDA-approved or clinically-advanced compounds (from octad.db's fda_drugs) that rank as reversers — a recovery / prioritization check."),
                                              DTOutput("rges_known_table")
                                    )
                                  ),
                                  br(),
                                  div(class = "alert alert-secondary d-flex align-items-start gap-2",
                                      bs_icon("rocket-takeoff", class = "mt-1"),
                                      div(
                                        strong("Looking for de novo compounds, not just approved drugs? "),
                                        "This tab screens the existing LINCS library of known compounds for repurposing. ",
                                        "For ", strong("virtual discovery of novel therapeutics"), " — predicting gene-expression effects directly ",
                                        "from a compound's chemical structure and designing new candidates — explore our lab's ",
                                        tags$a(href = "https://apps.octad.org/GPS/", target = "_blank", "GPS platform"),
                                        " (Gene expression profile Predictor on chemical Structures), published in ",
                                        tags$em("Cell"), " (2026). ",
                                        tags$a(href = "https://www.cell.com/cell/fulltext/S0092-8674(26)00223-0", target = "_blank", "Read the paper")
                                      )
                                  )
                        )
                      )
             )   # end Drug Repurposing tabPanel
             
  ), # end Drug Explorer navbarMenu
  
  # ── Disease Indication ──────────────────────────────────────────────────────
  tabPanel("Indication Expansion",
           sidebarLayout(
             sidebarPanel(width = 3,
                          h5(bs_icon("upload"), " Drug gene signature"),
                          selectInput("indic_id_type", "Gene ID type", choices = c("Gene symbol" = "symbol", "Ensembl ID" = "ensembl")),
                          textAreaInput("indic_up_genes",   "Up-regulated genes (one per line, or comma-separated)",   rows = 5, placeholder = "e.g.\nTP53\nEGFR\nMYC"),
                          textAreaInput("indic_down_genes", "Down-regulated genes (one per line, or comma-separated)", rows = 5, placeholder = "e.g.\nGLP1R\nINS"),
                          helpText("At least one of the two lists is required."),
                          div(class = "border rounded p-2 mb-2", style = "background: #f8f9fa;",
                              div(class = "small fw-bold mb-1", bs_icon("magic"), " Load an example (Semaglutide/GLP-1, Chen Lab)"),
                              selectInput("indic_example_tissue", "Tissue", choices = names(example_drug_files)),
                              selectInput("indic_example_topn", "Top N genes (each direction)",
                                          choices = c(10, 20, 50, 100), selected = 10),
                              actionButton("indic_load_example", "Load into gene lists", class = "btn-sm btn-outline-primary w-100")
                          ),
                          hr(),
                          
                          h5(bs_icon("geo-alt"), " Organ / Region"),
                          selectInput("indic_tissue", NULL, choices = c("All tissues" = "all", tissue_choices)),
                          helpText("Screening is per-dataset (per tissue), since meta-analysis pools tissues together."),
                          hr(),
                          
                          h5(bs_icon("sliders"), " Disease gene-set thresholds"),
                          numericInput("indic_padj_cutoff",  "Disease padj cutoff", 0.05, min = 0, max = 1, step = 0.01),
                          numericInput("indic_logfc_cutoff", "Disease |log2FC| cutoff", 0.5, min = 0, step = 0.1),
                          numericInput("indic_max_genes",    "Max disease genes per direction", 500, min = 20, step = 50),
                          hr(),
                          
                          h5(bs_icon("bar-chart-steps"), " Significance"),
                          numericInput("indic_p_cutoff", "Overlap p-value cutoff", 0.05, min = 0, max = 1, step = 0.01),
                          helpText("Fisher's-combined hypergeometric test on gene-list overlap. Below cutoff = significant reversal (Yes) or similarity (No)."),
                          hr(),
                          
                          div(class = "text-muted small mb-2", bs_icon("clock"), " Typically takes 1–2 minutes across all diseases (per-dataset, hypergeometric test — no permutation)."),
                          actionButton("run_indic", "Screen diseases", class = "btn-primary w-100"),
                          br(), br(),
                          downloadButton("download_indic", "Download results", class = "btn-sm btn-outline-primary w-100")
             ),
             mainPanel(width = 9,
                       br(),
                       uiOutput("indic_banner"),
                       br(),
                       navset_tab(
                         nav_panel("Ranked diseases",
                                   br(),
                                   helpText("Yes (1) = drug gene lists significantly reverse this disease-tissue signature. No (-1) = significantly same-direction (may worsen). Neutral (0) = no significant overlap."),
                                   DTOutput("indic_table")
                         ),
                         nav_panel("Bar chart",
                                   br(),
                                   numericInput("indic_top_n", "Diseases to plot", 30, min = 5, step = 5),
                                   plotlyOutput("indic_plot", height = "700px")
                         )
                       )
             )
           )
  ),   # end Disease Indication tabPanel
  
  # ── Documentation ────────────────────────────────────────────────────────────
  tabPanel("Documentation",
           div(class = "container mt-4",
               h3("Documentation"),
               p(class = "text-muted", "Brief description of each module and its key parameters."),
               br(),
               
               h4(bs_icon("bar-chart-line"), " Signature Explorer"),
               p("Visualizes one disease signature at a time: volcano plot, top differentially expressed genes, and summary statistics."),
               tags$ul(
                 tags$li(strong("Adjusted p-value cutoff"), " — padj threshold defining a gene as significant (default 0.05)."),
                 tags$li(strong("|log2FC| cutoff"), " — minimum absolute fold-change magnitude to call a gene up/down (default 1)."),
                 tags$li(strong("Label top N"), " — number of top up/down genes labeled on the volcano plot.")
               ),
               br(),
               
               h4(bs_icon("search"), " Gene Explorer"),
               p("Two sub-views: search a gene (or comma-separated list) across all disease signatures as a ranked heatmap, or trace a single gene's log2FC across all diseases as a forest plot."),
               tags$ul(
                 tags$li(strong("Across diseases vs across datasets"), " — averages multiple datasets per disease into one row, or shows each dataset separately."),
                 tags$li(strong("Rank by"), " — orders the heatmap/forest plot rows by strongest up, strongest down, or most significant.")
               ),
               br(),
               
               h4(bs_icon("diagram-2"), " Disease Similarity"),
               p("Correlates one query disease-dataset signature against every other dataset in the database (Spearman or Pearson, on shared genes)."),
               tags$ul(
                 tags$li(strong("Min shared genes"), " — minimum gene overlap required to report a correlation, avoiding spurious correlations from small overlaps."),
                 tags$li("Strongly ", strong("negative"), " correlations flag diseases with opposing transcriptional programs — useful as a proxy drug-context search.")
               ),
               br(),
               
               h4(bs_icon("signpost-split"), " Pathway Enrichment"),
               p("Runs enrichR against GO, KEGG, Reactome, WikiPathways, or MSigDB Hallmark for the up- or down-regulated gene set of a chosen disease/direction."),
               br(),
               
               h4(bs_icon("capsule"), " Drug Repurposing"),
               p("Scores compounds in the LINCS library by their capacity to reverse a disease's transcriptional signature, using the reverse gene expression score (sRGES) framework (", tags$a(href="https://doi.org/10.1038/s41467-018-06715-7", target="_blank", "Chen et al. 2018"), ")."),
               tags$ul(
                 tags$li(strong("Meta-analysis signature vs single dataset"), " — by default, diseases with multiple datasets use a random-effects meta-analyzed consensus signature; you can instead pick one raw dataset."),
                 tags$li(strong("Max heterogeneity (I²) / Min datasets per gene"), " — only apply in meta-analysis mode; filter out genes with high cross-dataset disagreement, or require a gene to be measured in enough datasets."),
                 tags$li(strong("|log2FC| / padj cutoffs"), " — define the disease's up/down gene sets used for scoring."),
                 tags$li(strong("Max genes per direction"), " — caps how many top up/down genes feed the scoring algorithm."),
                 tags$li(strong("Negative sRGES"), " = reversal = repurposing candidate. Results are annotated with clinical development phase where available.")
               ),
               br(),
               
               h4(bs_icon("bullseye"), " Indication Expansion"),
               p("The inverse query: given a drug's own gene signature (an arbitrary list of up- and/or down-regulated genes), screens every disease-tissue dataset in the database and classifies each as a significant reversal, significant concordance, or no significant association."),
               tags$ul(
                 tags$li(strong("Gene ID type"), " — whether your gene list uses gene symbols or Ensembl IDs."),
                 tags$li(strong("Organ/Region"), " — restrict the screen to a specific tissue, inferred from each dataset's source; datasets without a matched tissue label are always included under \"All tissues\"."),
                 tags$li(strong("Disease padj / |log2FC| cutoffs, max genes per direction"), " — define each disease's up/down gene sets, same logic as Drug Repurposing."),
                 tags$li(strong("Overlap p-value cutoff"), " — significance threshold for the Fisher's-combined hypergeometric enrichment test used to classify Yes/No/Neutral.")
               ),
               br(),
               
               h4(bs_icon("download"), " Download"),
               p("Download any individual signature (filtered by your chosen thresholds) or the full portal signature index as CSV."),
               br(),
               
               h4(bs_icon("bar-chart-steps"), " Meta-analysis"),
               p("Diseases with 2+ datasets are pooled with a DerSimonian-Laird random-effects model, run per gene across datasets. Reported metrics: pooled log2FC and p-value, ",
                 tags$strong("τ²"), " (between-dataset variance), ", tags$strong("I²"), " (percentage of variance due to heterogeneity rather than chance), and the number of contributing datasets per gene."),
               br()
           )
  ),
  
  # ── Download ────────────────────────────────────────────────────────────────
  tabPanel("Download",
           div(class = "container mt-4",
               h3("Download Signatures"),
               sidebarLayout(
                 sidebarPanel(width = 3,
                              selectizeInput("download_disease", "Disease",
                                             choices = c("Choose a disease…" = "", disease_choices_named),
                                             selected = "",
                                             options = list(placeholder = "Search disease…")),
                              uiOutput("download_signature_selector"),
                              br(),
                              downloadButton("download_selected_signature", "Download selected signature",
                                             class = "btn-primary w-100"),
                              br(), br(),
                              downloadButton("download_signature_index", "Download full index (CSV)",
                                             class = "btn-outline-secondary w-100"),
                              hr(),
                              h5(bs_icon("rocket-takeoff"), " GPS-ready format"),
                              checkboxInput("download_use_meta", "Use meta-analysis signature (recommended)", value = TRUE),
                              helpText("Two-column CSV (GeneSymbol, Value) matching the format expected by the ",
                                       tags$a(href = "https://apps.octad.org/GPS/", target = "_blank", "GPS platform"), "."),
                              downloadButton("download_gps_format", "Download for GPS", class = "btn-outline-success w-100")
                 ),
                 mainPanel(width = 9,
                           h5("Signature index"),
                           DTOutput("download_index")
                 )
               )
           )
  )
)  # end navbarPage

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # Enable URL bookmarking, but never restore disease/gene selections or file
  # uploads from a saved URL — those should always start blank for a fresh
  # user, even though the URL keeps updating with the current state (so
  # thresholds/parameters are still shareable/bookmarkable).
  session$setBookmarkExclude(c(
    "disease", "signature_file",
    "rges_disease", "rges_disease_mode", "rges_single_file", "rges_upload",
    "corr_disease", "corr_signature_file",
    "enrich_disease", "enrich_signature_file",
    "global_gene_search", "target_gene",
    "download_disease", "download_signature_file",
    "indic_up_genes", "indic_down_genes"
  ))
  
  observe({
    reactiveValuesToList(input)
    session$doBookmark()
  })
  onBookmarked(function(url) updateQueryString(url))
  
  # ── Disease Explorer helpers ─────────────────────────────────────────────────
  output$signature_selector <- renderUI({
    req(input$disease)
    files <- sig_index %>% filter(disease == input$disease) %>% pull(file)
    selectInput("signature_file", "Signature / dataset", choices = files)
  })
  
  sig_data <- reactive({
    req(input$disease, input$signature_file)
    get_sig(input$disease, input$signature_file)
  })
  
  
  
  filtered_data <- reactive({
    df <- sig_data() %>%
      filter(padj <= input$padj_cutoff, abs(log2FoldChange) >= input$logfc_cutoff)
    q <- trimws(input$gene_search)
    if (nchar(q) > 0) {
      df <- df %>% filter(
        str_detect(gene_symbol,      regex(q, ignore_case = TRUE)) |
          str_detect(identifier,       regex(q, ignore_case = TRUE)) |
          str_detect(gene_description, regex(q, ignore_case = TRUE))
      )
    }
    df
  })
  
  output$total_genes <- renderText({
    req(input$signature_file)
    nrow(sig_data())
  })
  output$sig_genes <- renderText({
    req(input$signature_file)
    sum(sig_data()$padj <= input$padj_cutoff, na.rm = TRUE)
  })
  output$up_genes <- renderText({
    req(input$signature_file)
    sum(sig_data()$padj <= input$padj_cutoff & sig_data()$log2FoldChange >= input$logfc_cutoff, na.rm = TRUE)
  })
  output$down_genes <- renderText({
    req(input$signature_file)
    sum(sig_data()$padj <= input$padj_cutoff & sig_data()$log2FoldChange <= -input$logfc_cutoff, na.rm = TRUE)
  })
  
  output$summary_table <- renderDT({
    df <- sig_data()
    datatable(data.frame(
      Disease          = pretty_disease(input$disease),
      Signature        = input$signature_file,
      Total_genes      = nrow(df),
      Significant      = sum(df$padj <= input$padj_cutoff, na.rm = TRUE),
      Up               = sum(df$padj <= input$padj_cutoff & df$log2FoldChange >= input$logfc_cutoff, na.rm = TRUE),
      Down             = sum(df$padj <= input$padj_cutoff & df$log2FoldChange <= -input$logfc_cutoff, na.rm = TRUE),
      Mean_abs_log2FC  = round(mean(abs(df$log2FoldChange), na.rm = TRUE), 4),
      Median_logCPM    = round(median(df$logCPM, na.rm = TRUE), 4)
    ), options = list(dom = "t", ordering = FALSE), rownames = FALSE)
  })
  
  output$direction_bar <- renderPlotly({
    df <- sig_data()
    summary_df <- data.frame(
      Direction = c("Up", "Down", "Not significant"),
      n = c(
        sum(df$padj <= input$padj_cutoff & df$log2FoldChange >= input$logfc_cutoff, na.rm = TRUE),
        sum(df$padj <= input$padj_cutoff & df$log2FoldChange <= -input$logfc_cutoff, na.rm = TRUE),
        sum(!(df$padj <= input$padj_cutoff & abs(df$log2FoldChange) >= input$logfc_cutoff), na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    )
    summary_df$Direction <- factor(summary_df$Direction, levels = c("Up","Down","Not significant"))
    ggplotly(
      ggplot(summary_df, aes(x = Direction, y = n, fill = Direction)) +
        geom_col(show.legend = FALSE) +
        scale_fill_manual(values = c("Up" = "#E74C3C", "Down" = "#2980B9", "Not significant" = "#BDC3C7")) +
        theme_bw(base_size = 13) +
        labs(x = NULL, y = "Number of genes", title = "Gene direction summary")
    )
  })
  
  output$volcano_plot <- renderPlotly({
    df <- sig_data() %>%
      mutate(
        neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
        status = case_when(
          padj <= input$padj_cutoff & log2FoldChange >=  input$logfc_cutoff ~ "Up",
          padj <= input$padj_cutoff & log2FoldChange <= -input$logfc_cutoff ~ "Down",
          TRUE ~ "NS"
        ),
        hover = paste0(
          "<b>", gene_symbol, "</b><br>",
          "Ensembl: ", identifier, "<br>",
          str_trunc(gene_description, 60), "<br>",
          "log2FC: ", round(log2FoldChange, 3), " | ",
          "-log10(padj): ", round(neg_log10_padj, 2), "<br>",
          "padj: ", signif(padj, 3)
        )
      )
    
    # Top-N labels
    top_label <- bind_rows(
      df %>% filter(status == "Up")   %>% arrange(padj, desc(log2FoldChange))  %>% head(input$label_n),
      df %>% filter(status == "Down") %>% arrange(padj,      log2FoldChange)   %>% head(input$label_n)
    )
    
    p <- ggplot(df, aes(log2FoldChange, neg_log10_padj, color = status, text = hover)) +
      geom_point(alpha = 0.6, size = 1.4) +
      geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(input$padj_cutoff), linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("Up" = "#E74C3C", "Down" = "#2980B9", "NS" = "#BDC3C7")) +
      theme_bw(base_size = 13) +
      labs(
        x = "log2 Fold Change", y = "-log10(adjusted p-value)",
        color = NULL,
        title = paste0(pretty_disease(input$disease), " — ", input$signature_file)
      )
    
    if (nrow(top_label) > 0) {
      p <- p + geom_text(
        data = top_label,
        aes(label = gene_symbol), color = "black", size = 3, hjust = -0.15, show.legend = FALSE
      )
    }
    
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", y = -0.1))
  })
  
  dt_opts <- list(pageLength = 10, scrollX = TRUE)
  
  output$top_up <- renderDT({
    sig_data() %>%
      filter(padj <= input$padj_cutoff) %>%
      arrange(desc(log2FoldChange)) %>%
      select(Gene = gene_symbol, Ensembl = identifier, Description = gene_description,
             log2FC = log2FoldChange, logCPM, LR, pvalue, padj) %>%
      head(input$top_n) %>%
      datatable(options = dt_opts, rownames = FALSE) %>%
      formatRound(c("log2FC", "logCPM", "LR"), 3) %>%
      formatSignif(c("pvalue", "padj"), 3) %>%
      formatStyle("log2FC", color = "#E74C3C")
  })
  
  output$top_down <- renderDT({
    sig_data() %>%
      filter(padj <= input$padj_cutoff) %>%
      arrange(log2FoldChange) %>%
      select(Gene = gene_symbol, Ensembl = identifier, Description = gene_description,
             log2FC = log2FoldChange, logCPM, LR, pvalue, padj) %>%
      head(input$top_n) %>%
      datatable(options = dt_opts, rownames = FALSE) %>%
      formatRound(c("log2FC", "logCPM", "LR"), 3) %>%
      formatSignif(c("pvalue", "padj"), 3) %>%
      formatStyle("log2FC", color = "#2980B9")
  })
  
  output$download_filtered <- downloadHandler(
    filename = function() paste0(tools::file_path_sans_ext(input$signature_file), "_filtered.csv"),
    content  = function(file) write.csv(filtered_data(), file, row.names = FALSE)
  )
  
  # ── Gene Search ─────────────────────────────────────────────────────────────
  global_gene_results <- eventReactive(input$run_gene_search, {
    req(input$global_gene_search)
    genes   <- str_trim(str_split(input$global_gene_search, ",")[[1]])
    genes   <- genes[nchar(genes) > 0]
    pattern <- paste(genes, collapse = "|")
    
    results <- list()
    withProgress(message = "Searching genes across signatures…", value = 0, {
      for (i in seq_len(nrow(sig_index))) {
        incProgress(1 / nrow(sig_index))
        df_i <- read_signature(sig_index$path[i]) %>%
          filter(
            str_detect(gene_symbol, regex(pattern, ignore_case = TRUE)) |
              str_detect(identifier,  regex(pattern, ignore_case = TRUE))
          )
        if (nrow(df_i) > 0) {
          df_i$disease        <- sig_index$disease[i]
          df_i$signature_file <- sig_index$file[i]
          df_i$signature_id   <- sig_index$signature_id[i]
          results[[length(results) + 1]] <- df_i
        }
      }
    })
    
    df <- bind_rows(results)
    if (input$gene_sig_only && nrow(df) > 0) {
      df <- df %>% filter(padj <= input$global_padj_cutoff, abs(log2FoldChange) >= input$global_logfc_cutoff)
    }
    
    # Across diseases mode: average log2FC per gene per disease
    if (!is.null(input$gene_search_mode) && input$gene_search_mode == "disease" && nrow(df) > 0) {
      df <- df %>%
        group_by(disease, gene_symbol) %>%
        summarise(
          log2FoldChange = mean(log2FoldChange, na.rm = TRUE),
          padj           = mean(padj, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(signature_id = disease)   # use disease as the row label
    }
    df
  })
  
  output$gene_search_banner <- renderUI({
    req(input$run_gene_search)
    df  <- global_gene_results()
    mode_label <- if (!is.null(input$gene_search_mode) && input$gene_search_mode == "disease")
      "averaged across datasets per disease" else "per dataset"
    div(class = "alert alert-info",
        bs_icon("info-circle"), " Searching for: ",
        strong(input$global_gene_search),
        " — found ", strong(nrow(df)), " entries (",
        mode_label, ")."
    )
  })
  
  output$global_gene_heatmap <- renderPlotly({
    df <- global_gene_results()
    validate(need(nrow(df) > 0, "No matching genes found with the current filters."))
    
    row_id <- if (!is.null(input$gene_search_mode) && input$gene_search_mode == "disease")
      "disease" else "signature_id"
    
    ranked_rows <- switch(input$gene_rank_by,
                          desc_fc  = df %>% group_by(.data[[row_id]]) %>%
                            summarise(rank_val = max(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
                            arrange(desc(rank_val)),
                          asc_fc   = df %>% group_by(.data[[row_id]]) %>%
                            summarise(rank_val = min(log2FoldChange, na.rm = TRUE), .groups = "drop") %>%
                            arrange(rank_val),
                          min_padj = df %>% group_by(.data[[row_id]]) %>%
                            summarise(rank_val = min(padj, na.rm = TRUE), .groups = "drop") %>%
                            arrange(rank_val)
    ) %>% head(input$gene_heatmap_top_n) %>% pull(.data[[row_id]])
    
    plot_df <- df %>%
      filter(.data[[row_id]] %in% ranked_rows) %>%
      select(row_label = all_of(row_id), gene_label = gene_symbol, log2FoldChange) %>%
      distinct() %>%
      pivot_wider(names_from = gene_label, values_from = log2FoldChange,
                  values_fn = mean)   # in case of duplicates, average
    
    mat <- plot_df %>% select(-row_label) %>% as.matrix()
    rownames(mat) <- plot_df$row_label
    mat <- mat[match(rev(ranked_rows[ranked_rows %in% plot_df$row_label]), rownames(mat)), , drop = FALSE]
    
    y_title <- if (!is.null(input$gene_search_mode) && input$gene_search_mode == "disease")
      "Disease (averaged across datasets, ranked)" else "Dataset signature (ranked)"
    
    plot_ly(
      z          = mat,
      x          = colnames(mat),
      y          = rownames(mat),
      type       = "heatmap",
      colorscale = list(c(0,"#2980B9"), c(0.5,"#FFFFFF"), c(1,"#E74C3C")),
      zmid       = 0,
      hovertemplate = "Gene: %{x}<br>Disease: %{y}<br>log2FC: %{z:.3f}<extra></extra>"
    ) %>%
      layout(
        xaxis  = list(title = "Gene", tickangle = -30, automargin = TRUE),
        yaxis  = list(title = y_title, automargin = TRUE),
        margin = list(l = 300, b = 80, r = 40, t = 40)
      )
  })
  
  output$download_global_gene_search <- downloadHandler(
    filename = function() "gene_search_results.csv",
    content  = function(file) write.csv(global_gene_results(), file, row.names = FALSE)
  )
  
  # ── Disease Similarity ──────────────────────────────────────────────────────
  output$corr_signature_selector <- renderUI({
    req(input$corr_disease)
    files <- sig_index %>% filter(disease == input$corr_disease) %>% pull(file)
    selectInput("corr_signature_file", "Signature / dataset", choices = files)
  })
  
  # Helper: get averaged signature for a disease across all its datasets
  get_averaged_sig <- function(disease_val) {
    files <- sig_index %>% filter(disease == disease_val) %>% pull(file)
    all_df <- lapply(files, function(f) {
      get_sig(disease_val, f) %>% select(gene_symbol, log2FoldChange)
    })
    bind_rows(all_df) %>%
      group_by(gene_symbol) %>%
      summarise(log2FoldChange = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")
  }
  
  output$corr_query_banner <- renderUI({
    req(input$corr_disease)
    div(class = "alert alert-secondary",
        bs_icon("diagram-2"), " Query: ",
        strong(pretty_disease(input$corr_disease)),
        if (!is.null(input$corr_signature_file)) paste0(" — ", input$corr_signature_file)
    )
  })
  
  corr_results <- eventReactive(input$run_corr, {
    req(input$corr_disease, input$corr_signature_file)
    
    query_df  <- get_sig(input$corr_disease, input$corr_signature_file) %>%
      select(identifier, query_log2FC = log2FoldChange)
    query_sid <- paste0(input$corr_disease, " | ", input$corr_signature_file)
    
    results <- list()
    withProgress(message = "Computing correlations across signatures…", value = 0, {
      for (i in seq_len(nrow(sig_index))) {
        incProgress(1 / nrow(sig_index))
        if (sig_index$signature_id[i] == query_sid) next
        df_i <- read_signature(sig_index$path[i]) %>%
          select(identifier, target_log2FC = log2FoldChange)
        merged   <- inner_join(query_df, df_i, by = "identifier")
        n_shared <- nrow(merged)
        if (n_shared < input$min_shared_genes) next
        corr_val <- suppressWarnings(
          cor(merged$query_log2FC, merged$target_log2FC,
              method = input$corr_method, use = "complete.obs")
        )
        if (is.na(corr_val)) next
        results[[length(results) + 1]] <- data.frame(
          target_disease = sig_index$disease[i],
          target_file    = sig_index$file[i],
          correlation    = corr_val, shared_genes = n_shared,
          stringsAsFactors = FALSE
        )
      }
    })
    
    res <- bind_rows(results)
    validate(need(nrow(res) > 0,
                  paste0("No matches found with 'Min shared genes' = ", input$min_shared_genes,
                         ". Try lowering this value in the sidebar.")))
    
    res %>%
      mutate(
        query_disease = input$corr_disease,
        query_file    = input$corr_signature_file,
        interpretation = case_when(
          correlation <= -0.3 ~ "Strong reversal (drug candidate context)",
          correlation <   0   ~ "Weak reversal",
          correlation >=  0.3 ~ "Strong similar biology",
          TRUE                ~ "Weak similar biology"
        )
      ) %>%
      arrange(desc(abs(correlation)))
  })
  
  output$corr_plot <- renderPlotly({
    req(corr_results())
    res <- corr_results()
    validate(need(nrow(res) > 0, "No diseases met the minimum shared genes threshold. Try lowering 'Min shared genes'."))
    plot_df <- bind_rows(
      res %>% arrange(desc(correlation)) %>% head(input$corr_top_n),
      res %>% arrange(correlation)       %>% head(input$corr_top_n)
    ) %>%
      distinct(target_disease, target_file, .keep_all = TRUE) %>%
      mutate(
        label     = paste0(pretty_disease(target_disease), "\n", target_file),
        direction = ifelse(correlation < 0, "Negative (reversal)", "Positive (similar)")
      )
    
    p <- ggplot(plot_df, aes(
      x    = reorder(str_trunc(paste(target_disease, target_file), 55), correlation),
      y    = correlation,
      fill = direction,
      text = paste0(
        "Disease: ", pretty_disease(target_disease), "<br>",
        "File: ", target_file, "<br>",
        "Correlation: ", round(correlation, 4), "<br>",
        "Shared genes: ", shared_genes, "<br>",
        interpretation
      )
    )) +
      geom_col() +
      coord_flip() +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
      scale_fill_manual(values = c("Negative (reversal)" = "#2980B9", "Positive (similar)" = "#E74C3C")) +
      theme_bw(base_size = 12) +
      labs(x = NULL, y = paste0(str_to_title(input$corr_method), " correlation"), fill = NULL,
           title = paste0("Similarity to: ", pretty_disease(input$corr_disease)))
    
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", y = -0.08))
  })
  
  output$corr_heatmap <- renderPlotly({
    req(corr_results())
    selected <- bind_rows(
      corr_results() %>% arrange(desc(correlation)) %>% head(input$corr_top_n),
      corr_results() %>% arrange(correlation)       %>% head(input$corr_top_n)
    ) %>% distinct(target_disease, target_file, .keep_all = TRUE)
    
    all_sids <- c(
      paste0(input$corr_disease, " | ", input$corr_signature_file),
      paste0(selected$target_disease, " | ", selected$target_file)
    )
    
    mat_list <- lapply(all_sids, function(sid) {
      parts <- strsplit(sid, " \\| ")[[1]]
      df <- get_sig(parts[1], parts[2]) %>%
        select(identifier, log2FoldChange) %>%
        as.data.frame()
      names(df)[names(df) == "log2FoldChange"] <- sid
      df
    })
    merged <- Reduce(function(x, y) full_join(x, y, by = "identifier"), mat_list)
    mat    <- as.matrix(merged[, -1])
    cm     <- cor(mat, use = "pairwise.complete.obs", method = input$corr_method)
    
    plot_ly(z = cm, x = colnames(cm), y = colnames(cm),
            type = "heatmap",
            colorscale = list(c(0,"#2980B9"), c(0.5,"#FFFFFF"), c(1,"#E74C3C")),
            zmin = -1, zmax = 1) %>%
      layout(
        height = max(800, 22 * ncol(cm)),
        margin = list(l = 320, b = 320, r = 40, t = 40),
        xaxis  = list(tickangle = 45, automargin = TRUE),
        yaxis  = list(automargin = TRUE)
      )
  })
  
  output$corr_table <- renderDT({
    corr_results() %>%
      select(Query_disease = query_disease, Query_file = query_file,
             Target_disease = target_disease, Target_file = target_file,
             Correlation = correlation, Shared_genes = shared_genes, Interpretation = interpretation) %>%
      datatable(options = list(pageLength = 25, scrollX = TRUE, filter = "top"), rownames = FALSE) %>%
      formatRound("Correlation", 4) %>%
      formatStyle("Correlation",
                  backgroundColor = styleInterval(c(-0.3, 0, 0.3),
                                                  c("#AED6F1", "#D6EAF8", "#FADBD8", "#F1948A")))
  })
  
  output$download_corr <- downloadHandler(
    filename = function() paste0(input$corr_disease, "_correlations.csv"),
    content  = function(file) write.csv(corr_results(), file, row.names = FALSE)
  )
  
  # ── Pathway Enrichment ──────────────────────────────────────────────────────
  output$enrich_signature_selector <- renderUI({
    req(input$enrich_disease)
    files <- sig_index %>% filter(disease == input$enrich_disease) %>% pull(file)
    selectInput("enrich_signature_file", "Signature / dataset", choices = files)
  })
  
  output$enrich_query_banner <- renderUI({
    req(input$enrich_disease)
    mode <- if (!is.null(input$enrich_mode)) input$enrich_mode else "disease"
    dir_label <- if (!is.null(input$enrich_direction))
      ifelse(input$enrich_direction == "up", "Up-regulated genes", "Down-regulated genes") else ""
    div(class = "alert alert-secondary",
        bs_icon("signpost-split"), " Enrichment for: ",
        strong(pretty_disease(input$enrich_disease)),
        if (mode == "disease") " — genes significant in ≥2 datasets"
        else if (!is.null(input$enrich_signature_file)) paste0(" — ", input$enrich_signature_file),
        " | ", strong(dir_label)
    )
  })
  
  enrich_results <- eventReactive(input$run_enrich, {
    req(input$enrich_disease)
    mode <- if (!is.null(input$enrich_mode)) input$enrich_mode else "disease"
    
    if (mode == "disease") {
      # Genes significant in at least 2 datasets for this disease
      files <- sig_index %>% filter(disease == input$enrich_disease) %>% pull(file)
      all_genes <- lapply(files, function(f) {
        df <- get_sig(input$enrich_disease, f)
        if (input$enrich_direction == "up") {
          df %>% filter(padj <= input$enrich_padj_cutoff, log2FoldChange >= input$enrich_logfc_cutoff) %>%
            pull(gene_symbol)
        } else {
          df %>% filter(padj <= input$enrich_padj_cutoff, log2FoldChange <= -input$enrich_logfc_cutoff) %>%
            pull(gene_symbol)
        }
      })
      gene_counts <- table(unlist(all_genes))
      min_datasets <- if (length(files) >= 2) 2 else 1
      genes <- names(gene_counts[gene_counts >= min_datasets])
    } else {
      req(input$enrich_signature_file)
      df <- get_sig(input$enrich_disease, input$enrich_signature_file)
      genes <- if (input$enrich_direction == "up") {
        df %>% filter(padj <= input$enrich_padj_cutoff, log2FoldChange >= input$enrich_logfc_cutoff) %>%
          pull(gene_symbol)
      } else {
        df %>% filter(padj <= input$enrich_padj_cutoff, log2FoldChange <= -input$enrich_logfc_cutoff) %>%
          pull(gene_symbol)
      }
    }
    
    genes <- unique(genes[!is.na(genes) & genes != "" & !str_detect(genes, "^ENSG")])
    validate(need(length(genes) >= 5, paste0(
      "Only ", length(genes), " genes found. Try relaxing thresholds or switching to single-dataset mode."
    )))
    
    tryCatch(
      enrichR::enrichr(genes, databases = input$enrich_db)[[input$enrich_db]],
      error = function(e) data.frame(Term = paste("Enrichment failed:", e$message),
                                     Overlap = NA, P.value = NA, Adjusted.P.value = NA,
                                     Odds.Ratio = NA, Combined.Score = NA, Genes = NA)
    )
  })
  
  output$enrich_plot <- renderPlotly({
    df <- enrich_results()
    validate(need(!all(is.na(df$Adjusted.P.value)), "No valid enrichment results."))
    plot_df <- df %>%
      arrange(Adjusted.P.value) %>% head(input$enrich_top_n) %>%
      mutate(neg_log10_fdr = -log10(pmax(Adjusted.P.value, .Machine$double.xmin)),
             Term_short = str_trunc(Term, 75))
    
    p <- ggplot(plot_df, aes(
      x    = reorder(Term_short, neg_log10_fdr),
      y    = neg_log10_fdr,
      fill = neg_log10_fdr,
      text = paste0("Pathway: ", Term, "<br>FDR: ", signif(Adjusted.P.value, 3),
                    "<br>Overlap: ", Overlap, "<br>Genes: ", str_trunc(Genes, 120))
    )) +
      geom_col(show.legend = FALSE) +
      scale_fill_gradient(low = "#AED6F1", high = "#1A5276") +
      coord_flip() +
      theme_bw(base_size = 12) +
      labs(x = NULL, y = "-log10(FDR)",
           title = paste0(pretty_disease(input$enrich_disease), " — ",
                          ifelse(input$enrich_direction == "up", "Up", "Down"), "-regulated | ",
                          input$enrich_db))
    
    ggplotly(p, tooltip = "text")
  })
  
  output$enrich_table <- renderDT({
    enrich_results() %>%
      select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes) %>%
      datatable(options = list(pageLength = 20, scrollX = TRUE, filter = "top"), rownames = FALSE) %>%
      formatSignif(c("P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score"), 3)
  })
  
  output$download_enrich <- downloadHandler(
    filename = function() paste0(input$enrich_disease, "_", input$enrich_direction, "_enrichment.csv"),
    content  = function(file) write.csv(enrich_results(), file, row.names = FALSE)
  )
  
  # ── GLP Explorer ────────────────────────────────────────────────────────────
  target_results <- eventReactive(input$run_target, {
    req(input$target_gene)
    gene    <- trimws(input$target_gene)
    pattern <- regex(paste0("^", gene, "$"), ignore_case = TRUE)
    
    results <- list()
    withProgress(message = "Scanning gene across signatures…", value = 0, {
      for (i in seq_len(nrow(sig_index))) {
        incProgress(1 / nrow(sig_index))
        df_i <- read_signature(sig_index$path[i]) %>%
          filter(str_detect(gene_symbol, pattern))
        if (nrow(df_i) > 0) {
          df_i$disease        <- sig_index$disease[i]
          df_i$signature_file <- sig_index$file[i]
          results[[length(results) + 1]] <- df_i
        }
      }
    })
    df <- bind_rows(results)
    
    validate(need(nrow(df) > 0, paste0("Gene '", gene, "' not found in any signature.")))
    
    if (input$target_sig_only) {
      df <- df %>% filter(padj <= input$target_padj_cutoff, abs(log2FoldChange) >= input$target_logfc_cutoff)
    }
    
    # Pick one row per disease (best padj)
    df <- df %>% group_by(disease, signature_file) %>%
      slice_min(padj, n = 1, with_ties = FALSE) %>% ungroup()
    
    df <- df %>% mutate(
      sig_hit  = padj <= input$target_padj_cutoff & abs(log2FoldChange) >= input$target_logfc_cutoff,
      label    = paste0(pretty_disease(disease), "\n(", signature_file, ")")
    )
    
    # Sort
    df <- switch(input$target_sort,
                 desc_fc = df %>% arrange(desc(log2FoldChange)),
                 asc_fc  = df %>% arrange(log2FoldChange),
                 alpha   = df %>% arrange(disease)
    )
    df %>% head(input$target_top_n)
  })
  
  output$target_banner <- renderUI({
    req(input$run_target)
    df <- target_results()
    n_sig <- sum(df$sig_hit, na.rm = TRUE)
    div(class = "alert alert-info",
        bs_icon("bullseye"), " Gene: ", strong(toupper(input$target_gene)),
        " found in ", strong(nrow(df)), " disease signatures. ",
        strong(n_sig), " pass significance thresholds."
    )
  })
  
  output$target_forest <- renderPlotly({
    df <- target_results()
    validate(need(nrow(df) > 0, "No data for this gene/filters."))
    
    # Compute approximate 95% CI from LR statistic (Wald-like: log2FC ± 1.96/sqrt(LR))
    # If LR column is missing or zero, CI width = 0
    df <- df %>% mutate(
      ci_half = ifelse(!is.na(LR) & LR > 0, 1.96 / sqrt(LR), 0),
      ci_lo   = log2FoldChange - ci_half,
      ci_hi   = log2FoldChange + ci_half,
      color   = case_when(
        sig_hit & log2FoldChange >  0 ~ "#E74C3C",
        sig_hit & log2FoldChange <  0 ~ "#2980B9",
        TRUE                          ~ "#95A5A6"
      ),
      shape = ifelse(sig_hit, "circle", "circle-open"),
      hover = paste0(
        "<b>", pretty_disease(disease), "</b><br>",
        "File: ", signature_file, "<br>",
        "log2FC: ", round(log2FoldChange, 3), "<br>",
        "95% CI: [", round(ci_lo, 3), ", ", round(ci_hi, 3), "]<br>",
        "padj: ", signif(padj, 3)
      )
    )
    
    fig <- plot_ly()
    
    # Error bars (CI)
    for (i in seq_len(nrow(df))) {
      fig <- fig %>% add_segments(
        x = df$ci_lo[i], xend = df$ci_hi[i],
        y = df$label[i], yend = df$label[i],
        line = list(color = df$color[i], width = 1.5),
        showlegend = FALSE, hoverinfo = "none"
      )
    }
    
    # Points
    fig <- fig %>%
      add_markers(
        data = df,
        x    = ~log2FoldChange,
        y    = ~label,
        color = ~sig_hit,
        colors = c("FALSE" = "#95A5A6", "TRUE" = "#E74C3C"),
        marker = list(size = 8),
        text   = ~hover,
        hovertemplate = "%{text}<extra></extra>"
      ) %>%
      add_segments(
        x = 0, xend = 0, y = 0.5, yend = nrow(df) + 0.5,
        line = list(color = "black", width = 1, dash = "dash"),
        showlegend = FALSE, hoverinfo = "none"
      ) %>%
      layout(
        title  = list(text = paste0(toupper(input$target_gene), " — log2FC forest plot across diseases")),
        xaxis  = list(title = "log2 Fold Change", zeroline = FALSE),
        yaxis  = list(title = "", automargin = TRUE, categoryorder = "array",
                      categoryarray = rev(df$label)),
        showlegend = FALSE,
        margin = list(l = 280, r = 40, t = 60, b = 60)
      )
    
    fig
  })
  
  output$target_table <- renderDT({
    target_results() %>%
      select(Disease = disease, Signature = signature_file,
             Gene = gene_symbol, log2FC = log2FoldChange,
             logCPM, LR, pvalue, padj, Significant = sig_hit) %>%
      datatable(options = list(pageLength = 15, scrollX = TRUE, filter = "top"), rownames = FALSE) %>%
      formatRound(c("log2FC", "logCPM", "LR"), 3) %>%
      formatSignif(c("pvalue", "padj"), 3) %>%
      formatStyle("Significant",
                  backgroundColor = styleEqual(c(TRUE, FALSE), c("#FADBD8", "#FDFEFE")))
  })
  
  output$download_target <- downloadHandler(
    filename = function() paste0(input$target_gene, "_across_diseases.csv"),
    content  = function(file) write.csv(target_results(), file, row.names = FALSE)
  )
  
  # ══ Drug Repurposing (RGES) ════════════════════════════════════════════════
  #
  # Powered by octad::runsRGES(). The disease signature (UPPERCASE Symbol +
  # log2FoldChange) is scored against the LINCS L1000 reference in octad.db;
  # sRGES < ~ -0.2 = meaningful reversal. Drug annotation (target, phase, MoA)
  # comes from octad.db's fda_drugs and cmpd_sets_mesh tables.
  #
  # NOTE: octad.db table/column names can vary by version. Helpers below use
  # flexible column lookup (pick_col) + tryCatch so a missing table degrades to
  # a friendly message instead of crashing.
  
  # Populate cell-line choices from octad.db's lincs_sig_info (EH7270)
  observe({
    cells <- tryCatch({
      tbl  <- load_octad_db("EH7270", "lincs_sig_info")
      if (is.null(tbl)) return(invisible(NULL))
      ccol <- pick_col(tbl, c("cell_id", "cell_line", "cell", "cell_iname"))
      if (is.na(ccol)) character(0) else sort(unique(toupper(as.character(tbl[[ccol]]))))
    }, error = function(e) character(0))
    if (!is.null(cells)) updateSelectizeInput(session, "rges_cells", choices = cells, server = TRUE)
  })
  
  # Flexible column picker: returns first matching column name present in df
  pick_col <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  # Load an octad.db object. Newer octad.db serves tables via ExperimentHub
  # (e.g. EH7269 = fda_drugs, EH7270 = lincs_sig_info), so try that first,
  # then fall back to data(); tolerate absence by returning NULL.
  load_octad_db <- function(eh_id = NULL, data_name = NULL) {
    obj <- NULL
    if (!is.null(eh_id)) {
      obj <- tryCatch(octad.db::get_ExperimentHub_data(eh_id), error = function(e) NULL)
    }
    if (is.null(obj) && !is.null(data_name)) {
      obj <- tryCatch({
        e <- new.env(); utils::data(list = data_name, package = "octad.db", envir = e)
        get(ls(e)[1], envir = e)
      }, error = function(e) NULL)
    }
    obj
  }
  
  output$rges_signature_selector_single <- renderUI({
    req(input$rges_disease)
    files <- sig_index %>% filter(disease == input$rges_disease) %>% pull(file)
    selectInput("rges_single_file", "Dataset", choices = files)
  })
  
  output$rges_signature_selector <- renderUI({
    req(input$rges_disease)
    files <- sig_index %>% filter(disease == input$rges_disease) %>% pull(file)
    n     <- length(files)
    mode  <- if (!is.null(input$rges_disease_mode)) input$rges_disease_mode else "meta"
    meta_exists <- !is.null(read_meta_signature(input$rges_disease))
    if (mode == "meta") {
      if (meta_exists) {
        div(class = "alert alert-info py-2 mb-1",
            bs_icon("info-circle"), " Meta-analysis signature (", strong(n), " dataset",
            if (n != 1) "s", " pooled via random-effects model)."
        )
      } else {
        div(class = "alert alert-warning py-2 mb-1",
            bs_icon("exclamation-triangle"), " No meta-analysis file found for this disease yet."
        )
      }
    } else if (n == 1) {
      helpText(paste("1 dataset:", files[1]))
    }
  })
  
  # ── Build RGES signature ────────────────────────────────────────────────────
  rges_raw_signature <- reactive({
    if (input$rges_source == "portal") {
      req(input$rges_disease)
      mode <- if (!is.null(input$rges_disease_mode)) input$rges_disease_mode else "meta"
      
      if (mode == "meta") {
        meta <- read_meta_signature(input$rges_disease)
        validate(need(!is.null(meta),
                      paste0("No meta-analysis signature found yet for ", pretty_disease(input$rges_disease),
                             ". Switch to 'Select one dataset', or generate/upload the meta signature.")))
        meta <- as.data.frame(meta)
        max_i2 <- if (!is.null(input$rges_max_i2)) input$rges_max_i2 else 100
        min_k  <- if (!is.null(input$rges_min_k))  input$rges_min_k  else 1
        out <- meta %>%
          filter(!is.na(gene_symbol), gene_symbol != "",
                 I2 <= max_i2, k_datasets >= min_k) %>%
          select(Symbol = gene_symbol, log2FoldChange = log2FoldChange_meta,
                 padj = padj_meta) %>%
          as.data.frame()
        out
      } else {
        req(input$rges_single_file)
        df <- get_sig(input$rges_disease, input$rges_single_file)
        data.frame(Symbol = df$gene_symbol, log2FoldChange = df$log2FoldChange,
                   padj = df$padj, stringsAsFactors = FALSE)
      }
    } else {
      req(input$rges_upload, input$rges_sym_col, input$rges_fc_col)
      up  <- fread(input$rges_upload$datapath) %>% as.data.frame()
      out <- data.frame(Symbol = up[[input$rges_sym_col]],
                        log2FoldChange = as.numeric(up[[input$rges_fc_col]]),
                        stringsAsFactors = FALSE)
      out$padj <- if (!is.null(input$rges_padj_col) && input$rges_padj_col %in% names(up))
        as.numeric(up[[input$rges_padj_col]]) else NA_real_
      out
    }
  })
  
  # Column mapper UI for uploaded files (feature 7)
  output$rges_col_mapper <- renderUI({
    req(input$rges_upload)
    cols <- names(fread(input$rges_upload$datapath, nrows = 0))
    tagList(
      selectInput("rges_sym_col",  "Gene symbol column",   choices = cols,
                  selected = cols[grep("symbol|gene", cols, ignore.case = TRUE)[1]]),
      selectInput("rges_fc_col",   "log2FC column",        choices = cols,
                  selected = cols[grep("log2|fc|fold", cols, ignore.case = TRUE)[1]]),
      selectInput("rges_padj_col", "padj column (optional)", choices = c("(none)", cols),
                  selected = cols[grep("padj|adj", cols, ignore.case = TRUE)[1]])
    )
  })
  
  # ── Construct up/down gene sets respecting thresholds + max genes ───────────
  # octad::runsRGES REQUIRES UPPERCASE HGNC symbols overlapping LINCS 978 genes.
  # When averaging across many datasets, padj is also averaged — use a relaxed
  # padj threshold or skip padj filtering if all values are NA.
  build_dz_sig <- function(sig_df) {
    sig <- sig_df %>%
      mutate(Symbol = toupper(trimws(as.character(Symbol)))) %>%
      filter(!is.na(Symbol), Symbol != "", !str_detect(Symbol, "^ENSG"),
             !is.na(log2FoldChange), is.finite(log2FoldChange))
    
    # Only apply padj filter if padj column has non-NA values
    has_padj <- "padj" %in% names(sig) && sum(!is.na(sig$padj)) > 0
    if (has_padj) {
      sig <- sig[sig$padj <= input$rges_padj_cutoff | is.na(sig$padj), ]
    }
    sig <- sig[abs(sig$log2FoldChange) >= input$rges_logfc_cutoff, ]
    
    up <- sig %>% filter(log2FoldChange > 0) %>% arrange(desc(log2FoldChange)) %>%
      head(input$rges_max_genes)
    dn <- sig %>% filter(log2FoldChange < 0) %>% arrange(log2FoldChange) %>%
      head(input$rges_max_genes)
    result <- as.data.frame(bind_rows(up, dn) %>% select(Symbol, log2FoldChange) %>% distinct())
    
    # Validate LINCS overlap — warn if too few genes will match
    if (nrow(result) < 10) {
      stop(paste0(
        "Only ", nrow(result), " genes passed filtering (need ≥10). ",
        "Try relaxing |log2FC| cutoff (currently ", input$rges_logfc_cutoff, ") ",
        "or padj cutoff (currently ", input$rges_padj_cutoff, ")."
      ))
    }
    result
  }
  
  # ── Core scoring wrapper ────────────────────────────────────────────────────
  # IMPORTANT: octad::runsRGES tests `missing(cells)` internally. Passing
  # cells = NULL is NOT the same as omitting it — runsRGES would then subset
  # lincs_sig_info to zero rows and crash ("replacement has 1 row, data has 0").
  # So we only include `cells` in the call when the user actually selected some,
  # and uppercase them (runsRGES uppercases its own cell_id before matching).
  run_srges <- function(dz_sig) {
    if (!.octad_available) {
      stop("The 'octad' package is not installed on this server. Install octad + octad.db (Bioconductor) to enable RGES scoring.")
    }
    if (nrow(dz_sig) < 4 || length(unique(sign(dz_sig$log2FoldChange))) < 2) {
      stop("Signature needs both up- and down-regulated genes after filtering. Relax the padj/log2FC thresholds.")
    }
    args <- list(
      dz_signature     = dz_sig,
      choose_fda_drugs = FALSE,
      max_gene_size    = input$rges_max_genes,
      permutations     = input$rges_permutations,
      output           = FALSE
    )
    if (length(input$rges_cells) > 0) {
      args$cells <- toupper(input$rges_cells)   # only when user selected cell lines
    }
    res <- do.call(octad::runsRGES, args)
    as.data.frame(res)
  }
  
  # ── Annotate drug hits with target / phase / MoA (features 2, 5, 6, 9) ──────
  # sRGES output keys drugs by `pert_iname`. octad.db ships `fda_drugs`
  # (name, target, clinical phase) and `cmpd_sets_mesh` (MeSH/MoA-like terms).
  annotate_drugs <- function(res) {
    if (is.null(res) || nrow(res) == 0) return(res)
    name_col <- pick_col(res, c("pert_iname", "Drug", "drug", "pert_name"))
    if (is.na(name_col)) { res$drug <- NA_character_; return(res) }
    res$drug <- as.character(res[[name_col]])
    
    # FDA drug table (EH7269): target + clinical phase
    fda <- load_octad_db("EH7269", "fda_drugs")
    
    if (!is.null(fda)) {
      fname  <- pick_col(fda, c("pert_iname", "name", "Drug", "drug"))
      ftgt   <- pick_col(fda, c("target", "Target", "gene_target"))
      fphase <- pick_col(fda, c("phase", "clinical_phase", "max_phase", "final.label", "status"))
      if (!is.na(fname)) {
        cols <- c(fname, ftgt, fphase); cols <- cols[!is.na(cols)]
        fda2 <- as.data.frame(fda)[, cols, drop = FALSE]
        names(fda2)[names(fda2) == fname]  <- "drug"
        if (!is.na(ftgt))   names(fda2)[names(fda2) == ftgt]   <- "target"
        if (!is.na(fphase)) names(fda2)[names(fda2) == fphase] <- "phase"
        fda2$drug <- toupper(as.character(fda2$drug))
        res$drug_key <- toupper(res$drug)
        res <- left_join(res, distinct(fda2, drug, .keep_all = TRUE),
                         by = c("drug_key" = "drug"))
        res$drug_key <- NULL
      }
    }
    
    # MeSH pharmacological terms as a MoA-like annotation (optional)
    moa <- load_octad_db(NULL, "cmpd_sets_mesh")
    if (!is.null(moa) && is.data.frame(moa)) {
      mname <- pick_col(moa, c("pert_iname", "name", "Drug", "drug"))
      mterm <- pick_col(moa, c("mesh", "MeSH", "term", "MOA", "moa"))
      if (!is.na(mname) && !is.na(mterm)) {
        m2 <- as.data.frame(moa)[, c(mname, mterm)]
        names(m2) <- c("drug", "MoA")
        m2$drug <- toupper(as.character(m2$drug))
        m2 <- m2 %>% group_by(drug) %>%
          summarise(MoA = paste(unique(MoA), collapse = "; "), .groups = "drop")
        res <- left_join(res, m2, by = c("drug" = "drug"))
      }
    }
    
    if (!"MoA"    %in% names(res)) res$MoA    <- NA_character_
    if (!"target" %in% names(res)) res$target <- NA_character_
    if (!"phase"  %in% names(res)) res$phase  <- NA_character_
    res
  }
  
  # ── Primary reactive: score disease A ───────────────────────────────────────
  rges_results <- eventReactive(input$run_rges, {
    dz <- build_dz_sig(rges_raw_signature())
    validate(need(nrow(dz) >= 10, "Need at least ~10 genes after filtering. Relax thresholds."))
    withProgress(message = "Scoring compounds against LINCS L1000 …", value = 0.5, {
      res <- tryCatch(run_srges(dz), error = function(e) {
        validate(need(FALSE, paste("RGES scoring failed:", conditionMessage(e))))
      })
      incProgress(0.4)
      annotate_drugs(res)
    })
  })
  
  srges_col <- function(res) pick_col(res, c("sRGES", "srges", "RGES", "score"))
  
  output$rges_banner <- renderUI({
    req(input$run_rges)
    if (input$rges_source == "portal") {
      req(input$rges_disease)
      n_files <- nrow(sig_index %>% filter(disease == input$rges_disease))
      div(class = "alert alert-info",
          bs_icon("capsule"), " Drug repurposing for: ", strong(pretty_disease(input$rges_disease)),
          if (n_files > 1)
            span(" — log2FC ", strong("averaged"), " across ", strong(n_files), " datasets before scoring.")
          else
            span(" — 1 dataset."),
          " Negative sRGES = signature reversal = candidate repurposing compound."
      )
    } else {
      div(class = "alert alert-info",
          bs_icon("capsule"), " Drug repurposing for: ", strong("uploaded signature"),
          " — negative sRGES = signature reversal = candidate repurposing compound."
      )
    }
  })
  
  output$rges_n_drugs  <- renderText({
    req(input$run_rges > 0)
    nrow(rges_results())
  })
  output$rges_n_strong <- renderText({
    req(input$run_rges > 0)
    res <- rges_results(); sc <- srges_col(res)
    if (is.na(sc)) "—" else sum(res[[sc]] <= -0.2, na.rm = TRUE)
  })
  output$rges_n_fda <- renderText({
    req(input$run_rges > 0)
    res <- rges_results()
    if (!"phase" %in% names(res)) return("—")
    sum(str_detect(tolower(as.character(res$phase)), "launch|approved|4"), na.rm = TRUE)
  })
  
  output$rges_table <- renderDT({
    res <- rges_results()
    sc  <- srges_col(res)
    validate(need(!is.na(sc), "No sRGES column found in scoring output."))
    df <- res %>% arrange(.data[[sc]])
    if (isTRUE(input$rges_fda_only) && "phase" %in% names(df)) {
      df <- df %>% filter(str_detect(tolower(as.character(phase)), "launch|approved|4"))
    }
    show_cols <- intersect(c("drug", sc, "MoA", "target", "phase"), names(df))
    datatable(df[, show_cols, drop = FALSE],
              options = list(pageLength = 25, scrollX = TRUE, filter = "top"),
              rownames = FALSE) %>%
      formatRound(sc, 4) %>%
      formatStyle(sc, backgroundColor = styleInterval(c(-0.3, -0.1, 0.1),
                                                      c("#196F3D", "#ABEBC6", "#FADBD8", "#F1948A")),
                  color = styleInterval(-0.3, c("white", "black")))
  })
  
  output$download_rges <- downloadHandler(
    filename = function() {
      base <- if (input$rges_source == "portal") input$rges_disease else "uploaded"
      paste0(base, "_drug_candidates_sRGES.csv")
    },
    content = function(file) write.csv(rges_results(), file, row.names = FALSE)
  )
  
  # ── Reversal plot (feature 3) ───────────────────────────────────────────────
  output$rges_drug_picker <- renderUI({
    res <- rges_results(); sc <- srges_col(res)
    top <- res %>% arrange(.data[[sc]]) %>% head(30) %>% pull(drug)
    selectInput("rges_pick_drug", "Compound", choices = top, width = "400px")
  })
  
  output$rges_reversal_plot <- renderPlotly({
    req(input$rges_pick_drug)
    dz <- build_dz_sig(rges_raw_signature()) %>%
      arrange(desc(log2FoldChange)) %>%
      mutate(rank = row_number(),
             direction = ifelse(log2FoldChange > 0, "Disease UP", "Disease DOWN"))
    # Conceptual reversal view: disease genes ranked by disease log2FC, colored by
    # direction. A reversing drug pushes UP genes down and DOWN genes up.
    p <- ggplot(dz, aes(rank, log2FoldChange, color = direction,
                        text = paste0(Symbol, "<br>disease log2FC: ", round(log2FoldChange, 2)))) +
      geom_segment(aes(xend = rank, yend = 0), alpha = 0.5) +
      geom_point(size = 1.6) +
      scale_color_manual(values = c("Disease UP" = "#E74C3C", "Disease DOWN" = "#2980B9")) +
      theme_bw(base_size = 12) +
      labs(x = "Disease signature gene rank", y = "Disease log2FC", color = NULL,
           title = paste0("Reversal target genes for ", input$rges_pick_drug))
    ggplotly(p, tooltip = "text")
  })
  
  # ── Known drugs cross-reference ─────────────────────────────────────────────
  output$rges_known_table <- renderDT({
    res <- rges_results(); sc <- srges_col(res)
    validate(need("phase" %in% names(res) && any(!is.na(res$phase)),
                  "No clinical-phase annotation available from octad.db's fda_drugs table."))
    out <- res %>%
      filter(!is.na(phase),
             str_detect(tolower(as.character(phase)), "launch|approved|phase|4|3")) %>%
      arrange(.data[[sc]]) %>%
      select(drug, all_of(sc), any_of(c("target", "phase", "MoA")))
    validate(need(nrow(out) > 0, "No clinically-advanced compounds among the scored hits."))
    datatable(out, options = list(pageLength = 15, scrollX = TRUE, filter = "top"),
              rownames = FALSE) %>%
      formatRound(sc, 4)
  })
  
  # ── Disease Indication (drug -> disease screening) ──────────────────────────
  # Ports the CMap KS-statistic connectivity score exactly as provided (RTF
  # script, Feb 2026 update). Faithful port: `rank` = -1 * rank(log2FC), i.e.
  # the drug's own gene ranking, tested for disease up/down gene enrichment.
  # ── Disease Indication: hypergeometric overlap test (binary: reversal/similar/neutral) ──
  # Fisher's-combined hypergeometric enrichment on gene-list membership only (no
  # fold-change magnitude needed from the drug side). Tissue-specific: operates
  # per-dataset (per tissue), since meta-analysis pools tissues together.
  hyper_p <- function(overlap, listA_size, listB_size, universe_size) {
    if (listA_size == 0 || listB_size == 0 || universe_size == 0) return(1)
    phyper(overlap - 1, listA_size, universe_size - listA_size, listB_size, lower.tail = FALSE)
  }
  combine_p <- function(p_vec) {
    p_vec <- p_vec[!is.na(p_vec)]
    if (length(p_vec) == 0) return(1)
    p_vec <- pmax(p_vec, 1e-300)
    if (length(p_vec) == 1) return(p_vec[1])
    chi <- -2 * sum(log(p_vec))
    pchisq(chi, df = 2 * length(p_vec), lower.tail = FALSE)
  }
  
  # Parse a manual gene-list textarea into a clean character vector
  parse_gene_list <- function(txt) {
    if (is.null(txt) || trimws(txt) == "") return(character(0))
    genes <- unlist(strsplit(txt, "[,\n\r\t ]+"))
    genes <- trimws(genes)
    unique(genes[genes != ""])
  }
  
  # Optional example: top N up / top N down genes from any Semaglutide
  # (GLP-1 agonist) tissue signature (Chen Lab), user-selected tissue and N.
  # Only populates the text areas when the user explicitly clicks the button
  # — never a default.
  observeEvent(input$indic_load_example, {
    req(input$indic_example_tissue, input$indic_example_topn)
    ex_path <- example_drug_files[[input$indic_example_tissue]]
    if (is.null(ex_path) || !file.exists(ex_path)) {
      showNotification("Example file not found on the server yet.", type = "warning")
      return()
    }
    top_n <- as.integer(input$indic_example_topn)
    df <- fread(ex_path) %>% as.data.frame()
    validate_cols <- all(c("Human_Symbol", "log2FoldChange", "padj") %in% names(df))
    if (!validate_cols) {
      showNotification("Example file is missing expected columns.", type = "error")
      return()
    }
    sig <- df %>% filter(!is.na(Human_Symbol), Human_Symbol != "", !is.na(padj))
    up   <- sig %>% filter(log2FoldChange > 0) %>% arrange(padj) %>% pull(Human_Symbol) %>% head(top_n)
    down <- sig %>% filter(log2FoldChange < 0) %>% arrange(padj) %>% pull(Human_Symbol) %>% head(top_n)
    updateTextAreaInput(session, "indic_up_genes",   value = paste(up, collapse = "\n"))
    updateTextAreaInput(session, "indic_down_genes", value = paste(down, collapse = "\n"))
    updateSelectInput(session, "indic_id_type", selected = "symbol")
    showNotification(paste0("Loaded top ", length(up), " up / ", length(down), " down genes from ",
                            input$indic_example_tissue, "."), type = "message")
  })
  
  # Build the drug's up/down gene lists (character vectors) + the ID type used
  indic_drug_lists <- reactive({
    up   <- parse_gene_list(input$indic_up_genes)
    down <- parse_gene_list(input$indic_down_genes)
    validate(need(length(up) > 0 || length(down) > 0,
                  "Enter at least one gene in the up-regulated or down-regulated list."))
    id_type <- if (!is.null(input$indic_id_type)) input$indic_id_type else "symbol"
    if (id_type == "symbol") { up <- toupper(up); down <- toupper(down) }
    list(up = up, down = down, id_type = id_type)
  })
  
  output$indic_banner <- renderUI({
    req(input$run_indic)
    dl <- indic_drug_lists()
    tiss <- if (!is.null(input$indic_tissue) && input$indic_tissue != "all")
      input$indic_tissue else "all tissues"
    div(class = "alert alert-info",
        bs_icon("upload"), " Screening ", strong(length(dl$up)), " up-regulated and ",
        strong(length(dl$down)), " down-regulated drug genes against disease signatures in ", strong(tiss), "."
    )
  })
  
  indic_results <- eventReactive(input$run_indic, {
    dl <- indic_drug_lists()
    id_col <- if (dl$id_type == "symbol") "gene_symbol" else "identifier"
    
    targets <- sig_index
    if (!is.null(input$indic_tissue) && input$indic_tissue != "all") {
      targets <- targets %>% filter(tissue == input$indic_tissue)
    }
    validate(need(nrow(targets) > 0, "No datasets match the selected tissue."))
    
    results <- list()
    withProgress(message = "Screening diseases…", value = 0, {
      n_total <- nrow(targets)
      for (i in seq_len(n_total)) {
        incProgress(1 / n_total)
        df <- get_sig(targets$disease[i], targets$file[i])
        id_vec <- df[[id_col]]
        if (dl$id_type == "symbol") id_vec <- toupper(id_vec)
        
        universe <- unique(id_vec[!is.na(id_vec) & id_vec != ""])
        universe_size <- length(universe)
        if (universe_size < 50) next
        
        sig_df <- df %>% mutate(.id = id_vec) %>%
          filter(!is.na(.id), .id != "", !is.na(padj), padj <= input$indic_padj_cutoff,
                 abs(log2FoldChange) >= input$indic_logfc_cutoff)
        dz_up   <- head(sig_df$.id[sig_df$log2FoldChange > 0], input$indic_max_genes)
        dz_down <- head(sig_df$.id[sig_df$log2FoldChange < 0], input$indic_max_genes)
        if (length(dz_up) < 5 && length(dz_down) < 5) next
        
        drug_up   <- intersect(dl$up,   universe)
        drug_down <- intersect(dl$down, universe)
        if (length(drug_up) == 0 && length(drug_down) == 0) next
        
        # Reversal: drug-up hits disease-down, drug-down hits disease-up
        p_rev <- combine_p(c(
          if (length(drug_up) > 0 && length(dz_down) > 0)
            hyper_p(length(intersect(drug_up, dz_down)), length(dz_down), length(drug_up), universe_size),
          if (length(drug_down) > 0 && length(dz_up) > 0)
            hyper_p(length(intersect(drug_down, dz_up)), length(dz_up), length(drug_down), universe_size)
        ))
        # Similarity: drug-up hits disease-up, drug-down hits disease-down
        p_sim <- combine_p(c(
          if (length(drug_up) > 0 && length(dz_up) > 0)
            hyper_p(length(intersect(drug_up, dz_up)), length(dz_up), length(drug_up), universe_size),
          if (length(drug_down) > 0 && length(dz_down) > 0)
            hyper_p(length(intersect(drug_down, dz_down)), length(dz_down), length(drug_down), universe_size)
        ))
        
        p_cut <- if (!is.null(input$indic_p_cutoff)) input$indic_p_cutoff else 0.05
        score <- if (p_rev <= p_cut && p_rev < p_sim) 1L
        else if (p_sim <= p_cut && p_sim < p_rev) -1L
        else 0L
        label <- c("-1" = "No", "0" = "Neutral", "1" = "Yes")[as.character(score)]
        
        results[[length(results) + 1]] <- data.frame(
          Disease = targets$disease[i], Tissue = targets$tissue[i], Dataset = targets$file[i],
          Indication = label, Score = score,
          P_reversal = p_rev, P_similarity = p_sim,
          Disease_up_genes = length(dz_up), Disease_down_genes = length(dz_down),
          Drug_up_matched = length(drug_up), Drug_down_matched = length(drug_down),
          stringsAsFactors = FALSE
        )
      }
    })
    
    res <- bind_rows(results)
    validate(need(nrow(res) > 0, "No datasets passed the gene-set thresholds. Try relaxing padj/log2FC cutoffs."))
    res %>%
      mutate(Disease_pretty = pretty_disease(Disease)) %>%
      arrange(pmin(P_reversal, P_similarity))
  })
  
  output$indic_table <- renderDT({
    df <- indic_results()
    datatable(df %>% select(Disease = Disease_pretty, Tissue, Dataset, Indication, Score,
                            P_reversal, P_similarity, Disease_up_genes, Disease_down_genes,
                            Drug_up_matched, Drug_down_matched),
              options = list(pageLength = 20, scrollX = TRUE, filter = "top"),
              rownames = FALSE) %>%
      formatRound(c("P_reversal", "P_similarity"), 4) %>%
      formatStyle("Indication", backgroundColor = styleEqual(
        c("Yes", "No", "Neutral"), c("#D5F5E3", "#FADBD8", "#F4F6F6")))
  })
  
  output$indic_plot <- renderPlotly({
    df <- indic_results()
    top_n <- if (!is.null(input$indic_top_n)) input$indic_top_n else 30
    plot_df <- df %>%
      mutate(neglog10p = -log10(pmax(pmin(P_reversal, P_similarity), 1e-10)),
             row_label = paste0(Disease_pretty, " (", Tissue, ")")) %>%
      arrange(pmin(P_reversal, P_similarity)) %>%
      head(top_n)
    validate(need(nrow(plot_df) > 0, "No results to plot."))
    
    p <- ggplot(plot_df, aes(x = reorder(row_label, neglog10p), y = neglog10p, fill = Indication,
                             text = paste0(row_label, "<br>", Indication,
                                           "<br>-log10(p): ", round(neglog10p, 2)))) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c("Yes" = "#2980B9", "No" = "#E74C3C", "Neutral" = "#BDC3C7")) +
      theme_bw(base_size = 11) +
      labs(x = NULL, y = "-log10(p)", fill = "Indication", title = "Strongest disease-tissue matches")
    ggplotly(p, tooltip = "text")
  })
  
  output$download_indic <- downloadHandler(
    filename = function() "disease_indication_results.csv",
    content  = function(file) write.csv(indic_results(), file, row.names = FALSE)
  )
  
  # ── Download tab ─────────────────────────────────────────────────────────────
  output$download_signature_selector <- renderUI({
    req(input$download_disease)
    files <- sig_index %>% filter(disease == input$download_disease) %>% pull(file)
    selectInput("download_signature_file", "Signature / dataset", choices = files)
  })
  
  output$download_index <- renderDT({
    datatable(
      sig_index %>%
        mutate(Disease = pretty_disease(disease)) %>%
        select(Disease, File = file, Signature_ID = signature_id),
      options = list(pageLength = 25, scrollX = TRUE, filter = "top"),
      rownames = FALSE
    )
  })
  
  output$download_selected_signature <- downloadHandler(
    filename = function() input$download_signature_file,
    content  = function(file) {
      req(input$download_disease, input$download_signature_file)
      path <- sig_index %>%
        filter(disease == input$download_disease, file == input$download_signature_file) %>%
        pull(path) %>% .[1]
      file.copy(path, file, overwrite = TRUE)
    }
  )
  
  output$download_signature_index <- downloadHandler(
    filename = function() "signature_index.csv",
    content  = function(file) {
      write.csv(sig_index %>% select(disease, file), file, row.names = FALSE)
    }
  )
  
  output$download_gps_format <- downloadHandler(
    filename = function() {
      req(input$download_disease)
      paste0("DZSIG__", input$download_disease, ".csv")
    },
    content = function(file) {
      req(input$download_disease)
      use_meta <- if (!is.null(input$download_use_meta)) input$download_use_meta else TRUE
      
      df <- if (use_meta) {
        meta <- read_meta_signature(input$download_disease)
        validate(need(!is.null(meta),
                      paste0("No meta-analysis signature available for ", pretty_disease(input$download_disease),
                             ". Uncheck 'Use meta-analysis signature' to export a single dataset instead.")))
        as.data.frame(meta)
      } else {
        req(input$download_signature_file)
        get_sig(input$download_disease, input$download_signature_file)
      }
      
      out <- df %>%
        transmute(GeneSymbol = gene_symbol, Value = log2FoldChange) %>%
        filter(!is.na(GeneSymbol), GeneSymbol != "", !is.na(Value)) %>%
        distinct(GeneSymbol, .keep_all = TRUE)
      write.csv(out, file, row.names = FALSE)
    }
  )
}

# ── Launch ─────────────────────────────────────────────────────────────────────
shinyApp(ui, server, enableBookmarking = "server")