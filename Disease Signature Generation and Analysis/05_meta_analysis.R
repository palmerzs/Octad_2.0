# =============================================================================
# REM Meta-Analysis: Precompute one meta-signature per disease
# =============================================================================
# For every disease with 2+ datasets, this runs a DerSimonian-Laird random-
# effects meta-analysis PER GENE across all its datasets, and writes one
# meta-signature CSV per disease. Diseases with only 1 dataset are copied
# through unchanged (with meta columns filled in as k=1, tau2=0, I2=0).
#
# Run this on the server (or your Mac) with: Rscript meta_analysis.R
# Requires only base R + data.table + dplyr (no octad/Bioconductor needed).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# ---- CONFIG -----------------------------------------------------------------
signature_dir <- "Desktop/Research/binchenlab/disease_signature/output_files/signature_for_drug_discovery/signature_may20/"
output_dir    <- "Desktop/Research/binchenlab/disease_signature/output_files/signature_for_drug_discovery/meta_signatures/"
# ------------------------------------------------------------------------------

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Same file-scanning logic as the Shiny app's build_signature_index()
files <- list.files(signature_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
sig_index <- data.frame(
  path    = files,
  disease = basename(dirname(files)),
  file    = basename(files),
  stringsAsFactors = FALSE
)
diseases <- sort(unique(sig_index$disease))
cat("Found", length(diseases), "diseases,", nrow(sig_index), "signature files.\n")

# ---- Per-gene SE approximation from log2FC + pvalue --------------------------
# Signature files have: identifier, gene_symbol, gene_description,
# log2FoldChange, logCPM, LR, pvalue, padj  (no lfcSE column available)
# SE is back-derived from the two-sided p-value: SE = |log2FC| / |qnorm(1 - p/2)|
add_se <- function(df) {
  p <- pmax(df$pvalue, 1e-300)                 # avoid p=0 -> Inf z
  p <- pmin(p, 1 - 1e-16)
  z <- qnorm(1 - p / 2)
  z[z <= 0 | !is.finite(z)] <- NA_real_
  df$se <- abs(df$log2FoldChange) / z
  # genes with log2FC == 0 or non-finite SE get a large (uninformative) SE
  # instead of NA, so they still contribute negligibly rather than being dropped
  bad <- !is.finite(df$se) | df$se <= 0
  if (any(bad)) df$se[bad] <- max(df$se[!bad], na.rm = TRUE) * 10
  df
}

# ---- DerSimonian-Laird random-effects meta-analysis, vectorized per gene ----
# y  = effect sizes (log2FC) for one gene across k datasets
# v  = variances (se^2) for one gene across k datasets
rem_dl <- function(y, v) {
  k <- length(y)
  if (k == 1) {
    return(list(est = y[1], se = sqrt(v[1]), z = NA_real_, p = NA_real_,
                tau2 = 0, I2 = 0, k = 1L))
  }
  w_fe   <- 1 / v
  y_fe   <- sum(w_fe * y) / sum(w_fe)
  Q      <- sum(w_fe * (y - y_fe)^2)
  df     <- k - 1
  C      <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
  tau2   <- max(0, (Q - df) / C)
  w_re   <- 1 / (v + tau2)
  y_re   <- sum(w_re * y) / sum(w_re)
  se_re  <- sqrt(1 / sum(w_re))
  z      <- y_re / se_re
  p      <- 2 * (1 - pnorm(abs(z)))
  I2     <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  list(est = y_re, se = se_re, z = z, p = p, tau2 = tau2, I2 = I2, k = k)
}

# ---- Main loop: one meta-signature CSV per disease ---------------------------
for (d in diseases) {
  d_files <- sig_index$path[sig_index$disease == d]
  out_path <- file.path(output_dir, paste0(d, "_meta.csv"))
  
  if (length(d_files) == 1) {
    # Single dataset: copy through, fill meta columns trivially
    df <- fread(d_files[1])
    df[, `:=`(
      log2FoldChange_meta = log2FoldChange,
      se_meta   = NA_real_,
      z_meta    = NA_real_,
      pval_meta = pvalue,
      padj_meta = padj,
      tau2      = 0,
      I2        = 0,
      k_datasets = 1L
    )]
    fwrite(df, out_path)
    cat(d, ": 1 dataset, copied through.\n")
    next
  }
  
  # Multi-dataset: REM per gene_symbol
  all_ds <- lapply(d_files, function(f) {
    df <- fread(f)
    df[, gene_symbol := fifelse(is.na(gene_symbol) | gene_symbol == "", identifier, gene_symbol)]
    df <- as.data.frame(df)
    df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "" & !is.na(df$log2FoldChange) & !is.na(df$pvalue), ]
    add_se(df[, c("identifier", "gene_symbol", "gene_description", "log2FoldChange", "pvalue")])
  })
  combined <- bind_rows(all_ds)
  
  meta_rows <- combined %>%
    group_by(gene_symbol) %>%
    filter(n() >= 1) %>%
    group_modify(~ {
      r <- rem_dl(.x$log2FoldChange, .x$se^2)
      data.frame(
        identifier          = .x$identifier[1],
        gene_description    = .x$gene_description[1],
        log2FoldChange_meta = r$est,
        se_meta             = r$se,
        z_meta              = r$z,
        pval_meta           = r$p,
        tau2                = r$tau2,
        I2                  = r$I2,
        k_datasets          = r$k
      )
    }) %>%
    ungroup()
  
  # BH-adjust the pooled p-values across genes for this disease
  meta_rows$padj_meta <- p.adjust(meta_rows$pval_meta, method = "BH")
  meta_rows$log2FoldChange <- meta_rows$log2FoldChange_meta  # convenience alias
  meta_rows$pvalue <- meta_rows$pval_meta
  meta_rows$padj   <- meta_rows$padj_meta
  
  fwrite(meta_rows, out_path)
  cat(d, ":", length(d_files), "datasets ->", nrow(meta_rows), "genes meta-analyzed.\n")
}

cat("\nDone. Meta-signatures written to:", output_dir, "\n")