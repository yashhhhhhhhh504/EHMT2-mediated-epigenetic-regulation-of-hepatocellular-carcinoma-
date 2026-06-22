# degs.R
#
# Preliminary Analysis of Gene Expression
#
# Input:  salmon.merged.gene_counts_scaled.tsv  (lengthScaledTPM method)
# Requires a sample metadata CSV with at minimum columns:
#   sample_id  – must match column names in the counts file
#   condition  – the primary contrast variable (e.g. "treated" / "control")
#   Any additional covariates you wish to include in the design
#
# Outputs (written to OUT_DIR):
#   library_sizes.png
#   sample_distance_heatmap.png
#   pca_plot.png
#   dispersion_estimates.png
#   ma_plot.png
#   volcano_plot.png
#   heatmap_top50_degs.png
#   deseq2_results_all.csv
#   deseq2_results_sig.csv


# ── Notes on input choice ─────────────────────────────────────────────────────
#
# | File                                  | Use for DESeq2? | Use for TPM viz? |
# |---------------------------------------|-----------------|------------------|
# | salmon.merged.gene_counts.tsv         | BEST            | No               |
# | salmon.merged.gene_counts_scaled.tsv  | Good *          | No               |
# | salmon.merged.gene_tpm.tsv            | No              | BEST             |
# | salmon.merged.gene_tpm_scaled.tsv     | No              | Acceptable       |
#
# * lengthScaledTPM corrects for transcript-length bias across samples —
#   recommended when isoform usage differs between conditions.
# WARNING: NEVER feed raw or scaled TPM to DESeq2.
# ──────────────────────────────────────────────────────────────────────────────


# 1  Setup -------------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)       # non-overlapping volcano labels
library(pheatmap)      # sample-distance & top-DEG heatmaps
library(RColorBrewer)
library(sva)

set.seed(42)

# ── User-defined parameters (edit before running) ───────────────────────────

METADATA_FILE  <- "sample_metadata.csv"   # see header comment for required cols
COUNTS_FILE    <- "salmon.merged.gene_counts.tsv"

# DESeq2 design — ComBat takes on batch effects
DESIGN_FORMULA <- ~ condition

# Which levels to contrast  (numerator vs denominator)
CONTRAST_FACTOR    <- "condition"
CONTRAST_NUMERATOR <- "treatment"     # the group of interest
CONTRAST_DENOM     <- "control"     # the reference group
BATCH_VARIABLE <- "batch"

# Thresholds
PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD  <- 1        # |log2FoldChange| >= 1  (i.e. >= 2-fold)
MIN_COUNT      <- 10       # pre-filter: keep genes with rowSum >= this

# Output directory
OUT_DIR <- "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# ── Helper: save a ggplot ────────────────────────────────────────────────────

save_plot <- function(plot, filename, width = 8, height = 6, ...) {
  ggsave(
    filename = file.path(OUT_DIR, filename),
    plot = plot, width = width, height = height, dpi = 300, ...
  )
  message("Saved: ", file.path(OUT_DIR, filename))
}


# 2  Load data ---------------------------------------------------------------

message("Loading counts ...")
counts_raw <- read.delim(COUNTS_FILE, sep = "\t", header = TRUE, check.names = FALSE)

# nf-core salmon output has gene_id + gene_name as the first two columns
gene_info <- counts_raw[, c("gene_id", "gene_name")]
count_mat  <- counts_raw[, !(colnames(counts_raw) %in% c("gene_id", "gene_name"))]
rownames(count_mat) <- counts_raw$gene_id

# Round to integers (lengthScaledTPM values are near-integers; DESeq2 requires int)
count_mat <- round(count_mat) |> as.matrix()
storage.mode(count_mat) <- "integer"

# Pre-filter low-count genes (speeds up analysis, reduces multiple-testing burden)
keep       <- rowSums(count_mat) >= MIN_COUNT
count_mat  <- count_mat[keep, ]
message(sprintf("Retained %d / %d genes after filtering (rowSum >= %d)",
                sum(keep), length(keep), MIN_COUNT))

message("Loading metadata ...")
meta <- read.csv(METADATA_FILE, row.names = 1, stringsAsFactors = TRUE)

# Ensure sample order matches count matrix columns
meta <- meta[colnames(count_mat), , drop = FALSE]
stopifnot(all(rownames(meta) == colnames(count_mat)))

# Make sure the reference level is set correctly
meta[[CONTRAST_FACTOR]] <- relevel(meta[[CONTRAST_FACTOR]], ref = CONTRAST_DENOM)

# Apply ComBat-seq to correct batch effects in the count matrix
message("Applying ComBat-seq batch correction ...")
batch        <- meta[[BATCH_VARIABLE]]
mod          <- model.matrix(~ condition, data = meta)
count_mat    <- ComBat_seq(count_mat, batch = batch, group = NULL, covar_mod = mod)
storage.mode(count_mat) <- "integer"
message("Batch correction complete.")


# 3  Build DESeqDataSet ------------------------------------------------------

message("Building DESeqDataSet ...")
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = meta,
  design    = DESIGN_FORMULA
)


# 4  Quality-control / Exploratory analysis ----------------------------------

# 4a  VST for visualisation (blind = TRUE -> unsupervised, for QC only)
vst <- vst(dds, blind = TRUE)

# ── Sample-to-sample distance heatmap ───────────────────────────────────────

sample_dists  <- dist(t(assay(vst)))
sample_mat    <- as.matrix(sample_dists)
colours       <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

ann_col <- meta[, CONTRAST_FACTOR, drop = FALSE]   # colour sidebar by condition

png(file.path(OUT_DIR, "sample_distance_heatmap.png"),
    width = 2400, height = 2000, res = 300)
pheatmap(
  sample_mat,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col            = colours,
  annotation_col = ann_col,
  main           = "Sample-to-sample distances (VST)"
)
dev.off()
message("Saved: ", file.path(OUT_DIR, "sample_distance_heatmap.png"))

# ── PCA ─────────────────────────────────────────────────────────────────────

pca_data <- plotPCA(vst, intgroup = CONTRAST_FACTOR, returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                 colour = .data[[CONTRAST_FACTOR]],
                                 label  = name)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    title  = "PCA — VST-normalised counts",
    x      = paste0("PC1 (", pct_var[1], "% variance)"),
    y      = paste0("PC2 (", pct_var[2], "% variance)"),
    colour = CONTRAST_FACTOR
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

save_plot(pca_plot, "pca_plot.png")

# ── Library size bar chart ───────────────────────────────────────────────────

lib_sizes <- colSums(count_mat) |>
  tibble::enframe(name = "sample", value = "total_counts") |>
  left_join(tibble::rownames_to_column(meta, "sample"), by = "sample")

libsize_plot <- ggplot(lib_sizes,
                       aes(x    = reorder(sample, total_counts),
                           y    = total_counts / 1e6,
                           fill = .data[[CONTRAST_FACTOR]])) +
  geom_col() +
  coord_flip() +
  labs(title = "Library sizes", x = NULL,
       y = "Total counts (millions)", fill = CONTRAST_FACTOR) +
  theme_bw(base_size = 12)

save_plot(libsize_plot, "library_sizes.png", width = 7, height = 5)


# 5  Differential expression -------------------------------------------------

message("Running DESeq2 ...")
dds <- DESeq(dds)

# Extract results for the contrast of interest
res <- results(
  dds,
  contrast = c(CONTRAST_FACTOR, CONTRAST_NUMERATOR, CONTRAST_DENOM),
  alpha    = PADJ_THRESHOLD
)

# Shrink LFC estimates (apeglm is recommended; falls back to ashr if unavailable)
res_shrunk <- tryCatch(
  lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm"),
  error = function(e) {
    message("apeglm not available, falling back to ashr: ", e$message)
    lfcShrink(dds,
              contrast = c(CONTRAST_FACTOR, CONTRAST_NUMERATOR, CONTRAST_DENOM),
              type     = "ashr")
  }
)

# Annotate with gene symbols
res_df <- as.data.frame(res_shrunk) |>
  tibble::rownames_to_column("gene_id") |>
  left_join(gene_info, by = "gene_id") |>
  arrange(padj)

# Classify significance
res_df <- res_df |>
  mutate(
    sig = case_when(
      padj < PADJ_THRESHOLD & log2FoldChange >=  LFC_THRESHOLD ~ "Up",
      padj < PADJ_THRESHOLD & log2FoldChange <= -LFC_THRESHOLD ~ "Down",
      TRUE                                                       ~ "NS"
    ),
    sig = factor(sig, levels = c("Up", "Down", "NS"))
  )

n_up   <- sum(res_df$sig == "Up",   na.rm = TRUE)
n_down <- sum(res_df$sig == "Down", na.rm = TRUE)
message(sprintf("DEGs: %d up, %d down (padj < %.2f, |LFC| >= %d)",
                n_up, n_down, PADJ_THRESHOLD, LFC_THRESHOLD))

# Write full results table
write.csv(res_df, file.path(OUT_DIR, "deseq2_results_all.csv"), row.names = FALSE)

# Write significant DEGs only
degs <- filter(res_df, sig != "NS")
write.csv(degs, file.path(OUT_DIR, "deseq2_results_sig.csv"), row.names = FALSE)


# 6  Visualisations ----------------------------------------------------------

# ── MA plot ──────────────────────────────────────────────────────────────────

ma_plot <- ggplot(res_df |> filter(!is.na(padj)),
                  aes(x      = log10(baseMean + 1),
                      y      = log2FoldChange,
                      colour = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70")) +
  labs(title  = "MA plot",
       x      = "log10(mean expression + 1)",
       y      = "log2 fold change",
       colour = NULL) +
  theme_bw(base_size = 13)

save_plot(ma_plot, "ma_plot.png")

# ── Volcano plot ─────────────────────────────────────────────────────────────

# Top genes to label on the volcano
top_labels <- res_df |>
  filter(sig != "NS") |>
  slice_min(padj, n = 20)

volcano_plot <- ggplot(res_df |> filter(!is.na(padj)),
                       aes(x      = log2FoldChange,
                           y      = -log10(padj),
                           colour = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_text_repel(
    data           = top_labels,
    aes(label      = gene_name),
    size           = 2.5,
    max.overlaps   = 25,
    segment.colour = "grey50"
  ) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(PADJ_THRESHOLD),
             linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(Up = "#d73027", Down = "#4575b4", NS = "grey70")) +
  labs(
    title   = sprintf("Volcano: %s vs %s", CONTRAST_NUMERATOR, CONTRAST_DENOM),
    x       = "log2 fold change",
    y       = expression(-log[10](p[adj])),
    colour  = NULL,
    caption = sprintf("%d up  |  %d down  (padj < %.2f, |LFC| >= %d)",
                      n_up, n_down, PADJ_THRESHOLD, LFC_THRESHOLD)
  ) +
  theme_bw(base_size = 13)

save_plot(volcano_plot, "volcano_plot.png")

# ── Heatmap of top 50 DEGs ───────────────────────────────────────────────────

top50_ids <- degs |>
  slice_min(padj, n = 50) |>
  pull(gene_id)

if (length(top50_ids) > 0) {
  # Use VST (blind = FALSE) for supervised visualisation
  vst_sup   <- vst(dds, blind = FALSE)
  mat_top50 <- assay(vst_sup)[top50_ids, , drop = FALSE]
  
  # Z-score across samples
  mat_scaled <- t(scale(t(mat_top50)))
  
  # Row labels: gene symbol where available, else gene_id
  row_labels <- gene_info$gene_name[match(rownames(mat_scaled), gene_info$gene_id)]
  row_labels[is.na(row_labels) | row_labels == ""] <-
    rownames(mat_scaled)[is.na(row_labels) | row_labels == ""]
  rownames(mat_scaled) <- row_labels
  
  ann_col_sup <- meta[, CONTRAST_FACTOR, drop = FALSE]
  
  png(file.path(OUT_DIR, "heatmap_top50_degs.png"),
      width = 2800, height = 3600, res = 300)
  pheatmap(
    mat_scaled,
    annotation_col = ann_col_sup,
    show_colnames  = TRUE,
    fontsize_row   = 8,
    color          = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
    main           = paste("Top 50 DEGs — Z-scored VST counts\n",
                           CONTRAST_NUMERATOR, "vs", CONTRAST_DENOM)
  )
  dev.off()
  message("Saved: ", file.path(OUT_DIR, "heatmap_top50_degs.png"))
} else {
  message("Fewer than 1 significant DEG found — heatmap skipped.")
}

# ── Dispersion & size-factor diagnostics ─────────────────────────────────────

png(file.path(OUT_DIR, "dispersion_estimates.png"),
    width = 2000, height = 1600, res = 300)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()
message("Saved: ", file.path(OUT_DIR, "dispersion_estimates.png"))


# 7  Summary -----------------------------------------------------------------

message("\n── Analysis complete ─────────────────────────────────────────────────────")
message(sprintf("  Contrast  : %s  (numerator: %s  |  denominator: %s)",
                CONTRAST_FACTOR, CONTRAST_NUMERATOR, CONTRAST_DENOM))
message(sprintf("  Thresholds: padj < %.2f  |  |LFC| >= %d", PADJ_THRESHOLD, LFC_THRESHOLD))
message(sprintf("  Up-regulated  : %d genes", n_up))
message(sprintf("  Down-regulated: %d genes", n_down))
message(sprintf("  Results saved to: %s/", OUT_DIR))
message("──────────────────────────────────────────────────────────────────────────")

writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt"))