## 1. filtering_16sample_2nd_run

# 2nd run 16 sample multiome (Ob-CARD)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(DoubletFinder)
library(Matrix)
library(tidyr)
library(ggplot2)
library(dplyr)
set.seed(123)
library(DoubletFinder) # Version 2.0.6
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'  
genome(annotations) <- "hg38"

samples_df <- read.csv("/home/n/nb443/Documents/Ob_CARD2_snRNA-seq_analysis/16_sample_2nd_run_locations", header = TRUE)
# Rename columns to standard names
colnames(samples_df) <- c("sample_id", "path")
# Remove extra whitespace in paths
samples_df$path <- gsub(" ", "", samples_df$path)
# Remove any empty rows (optional safety check)
samples_df <- na.omit(samples_df)
samples_df <- samples_df[samples_df$sample_id != "" & samples_df$path != "", ]

#  Prepare TSS positions
tss.positions <- promoters(genes(EnsDb.Hsapiens.v86), upstream = 0, downstream = 1)
tss.positions <- keepStandardChromosomes(tss.positions, pruning.mode = "coarse")
seqlevelsStyle(tss.positions) <- "UCSC"

# --- Clean paths from CSV and keep only existing folders ---
samples_df <- read.csv("/home/n/nb443/Documents/Ob_CARD2_snRNA-seq_analysis/16_sample_2nd_run_locations",
                       header = FALSE, stringsAsFactors = FALSE)
colnames(samples_df) <- c("sample_id", "path")

samples_df$path <- gsub("\r", "", samples_df$path)          # remove stray CRs
samples_df$path <- trimws(samples_df$path)                  # trim whitespace
samples_df$path <- sub("/+$", "", samples_df$path)          # drop trailing slashes
samples_df <- na.omit(samples_df)
samples_df <- samples_df[samples_df$sample_id != "" & samples_df$path != "", ]
samples_df$path <- normalizePath(samples_df$path, mustWork = FALSE)
samples_df <- samples_df[dir.exists(samples_df$path), , drop = FALSE]
raw_objects <- list()

# make a named list to hold each unfiltered Seurat object
raw_objects <- vector("list", nrow(samples_df))
names(raw_objects) <- samples_df$sample_id
for (i in seq_len(nrow(samples_df))) {
  sid   <- samples_df$sample_id[i]
  dpath <- samples_df$path[i] 
  x <- Read10X(data.dir = dpath)
  outs_dir <- dirname(dpath)
  frag_candidates <- c(
    file.path(outs_dir, "atac_fragments.tsv.gz"),
    file.path(outs_dir, "fragments.tsv.gz"))
  fragfile <- frag_candidates[file.exists(frag_candidates)][1]
  chrom_assay <- CreateChromatinAssay(
    counts     = x[["Peaks"]],
    sep        = c(":", "-"),
    genome     = "hg38",
    fragments  = fragfile,
    annotation = annotations)
  so <- CreateSeuratObject(
    counts       = x[["Gene Expression"]],
    assay        = "RNA",
    project      = sid,
    min.cells    = 1,
    min.features = 0)
  so[["ATAC"]] <- chrom_assay
  so$sample <- sid
  DefaultAssay(so) <- "RNA"
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  DefaultAssay(so) <- "ATAC"
  so <- NucleosomeSignal(so)
  so <- TSSEnrichment(so, tss.positions = tss.positions)
  DefaultAssay(so) <- "RNA"
  raw_objects[[sid]] <- so}
# optional: merge all unfiltered samples
unfiltered_merged <- Reduce(function(x, y) merge(x, y), raw_objects)

# Violin plots for QC
sample_levels <- unique(samples_df$sample_id)
ord <- order(as.integer(sub(".*_", "", sample_levels)))
unfiltered_merged$sample <- factor(unfiltered_merged$sample, levels = sample_levels[ord])
# nFeature_RNA
VlnPlot(unfiltered_merged, features = "nFeature_RNA", group.by = "sample", pt.size = 0.1, alpha = 0.3) +
  theme(legend.position = "none") +
  ggtitle("nFeature_RNA") +
  ylab("Detected genes per cell") + xlab(NULL)
# percent.mt
VlnPlot(unfiltered_merged, features = "percent.mt", group.by = "sample", pt.size = 0.1, alpha = 0.3) +
  theme(legend.position = "none") +
  ggtitle("percent.mt") +
  ylab("Percent mitochondrial reads") + xlab(NULL)
# nCount_ATAC
VlnPlot(unfiltered_merged, features = "nCount_ATAC", group.by = "sample", pt.size = 0.1, alpha = 0.3) +
  theme(legend.position = "none") +
  ggtitle("nCount_ATAC") +
  ylab("Number of ATAC fragments per cell")  + xlab(NULL)
# nucleosome_signal
VlnPlot(unfiltered_merged, features = "nucleosome_signal", group.by = "sample", pt.size = 0.1, alpha = 0.3) +
  theme(legend.position = "none") +
  ggtitle("nucleosome_signal") +
  ylab("Nucleosome signal")  + xlab(NULL)
# TSS.enrichment
VlnPlot(unfiltered_merged, features = "TSS.enrichment", group.by = "sample", pt.size = 0.1, alpha = 0.3) +
  theme(legend.position = "none") +
  ggtitle("TSS.enrichment") +
  ylab("TSS enrichment score")+ xlab(NULL)

# keep the same sample order as your plots
sample_order <- if (is.factor(unfiltered_merged$sample)) levels(unfiltered_merged$sample) else sort(unique(unfiltered_merged$sample))
qc_means <- unfiltered_merged@meta.data %>%
  dplyr::select(sample,
                nFeature_RNA,
                percent.mt,
                nCount_ATAC,
                nucleosome_signal,
                TSS.enrichment) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_nFeature_RNA   = mean(nFeature_RNA,   na.rm = TRUE),
    mean_percent_mt     = mean(percent.mt,     na.rm = TRUE),
    mean_nCount_ATAC    = mean(nCount_ATAC,    na.rm = TRUE),
    mean_nucleosome_sig = mean(nucleosome_signal, na.rm = TRUE),
    mean_TSS_enrich     = mean(TSS.enrichment, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(sample = factor(sample, levels = sample_order)) %>%
  arrange(sample) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))  # optional rounding
qc_means

# --- cumulative filtering + per-sample counts at each step ---
# convenience: all sample IDs (so we include zeros explicitly)
all_samples <- sort(unique(unfiltered_merged$sample))
count_by_sample <- function(obj) {
  data.frame(
    sample = all_samples,
    cells_remaining = as.integer(table(factor(obj$sample, levels = all_samples))),
    stringsAsFactors = FALSE)}
# stepwise cumulative subsets
so0 <- unfiltered_merged
so1 <- subset(so0, subset = nFeature_RNA > 200)
so2 <- subset(so1, subset = nFeature_RNA < 12500)
so3 <- subset(so2, subset = percent.mt < 5)
so4 <- subset(so3, subset = nCount_ATAC > 13)
so5 <- subset(so4, subset = nCount_ATAC < 100000)
so6 <- subset(so5, subset = nucleosome_signal < 1)
so7 <- subset(so6, subset = TSS.enrichment > 1)  # final filtered
# assemble per-step counts (long table)
qc_counts <- rbind(
  transform(count_by_sample(so0), step = "Start"),
  transform(count_by_sample(so1), step = "nFeature_RNA > 200"),
  transform(count_by_sample(so2), step = "nFeature_RNA < 12500"),
  transform(count_by_sample(so3), step = "percent.mt < 5"),
  transform(count_by_sample(so4), step = "nCount_ATAC > 13"),
  transform(count_by_sample(so5), step = "nCount_ATAC < 100000"),
  transform(count_by_sample(so6), step = "nucleosome_signal < 1"),
  transform(count_by_sample(so7), step = "TSS.enrichment > 1"))
# tidy ordering of steps
qc_counts$step <- factor(
  qc_counts$step,
  levels = c(
    "Start",
    "nFeature_RNA > 200",
    "nFeature_RNA < 12500",
    "percent.mt < 5",
    "nCount_ATAC > 13",
    "nCount_ATAC < 100000",
    "nucleosome_signal < 1",
    "TSS.enrichment > 1"))
# view in console
qc_counts <- qc_counts[order(qc_counts$sample, qc_counts$step), ]
print(qc_counts)
qc_counts_wide <- pivot_wider(qc_counts, id_cols = sample, names_from = step, values_from = cells_remaining)
print(qc_counts_wide)
saveRDS(so7, file = "/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_2_filtered.RDS")

# DoubletFinder doublet removal
# split the filtered merged object by sample
seurat_objects <- SplitObject(seurat_objects, split.by = "sample")
seurat_objects <- readRDS("/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_2_filtered.RDS")

# obcard_146
obcard_146 <- seurat_objects[["obcard_146"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_146, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_146@meta.data))) * (1 - (modelHomotypic(obcard_146@meta.data$seurat_clusters))))
obcard_146 <- doubletFinder(obcard_146, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_146@meta.data)

obcard_146 <- subset(obcard_146, DF.classifications_0.25_0.29_6 != 'Doublet')
seurat_objects[["obcard_146"]] <- obcard_146
doubletest1 <- round((round(0.01 * nrow(obcard_146@meta.data))) * (1 - (modelHomotypic(obcard_146@meta.data$seurat_clusters))))
print(doubletest1)
table(obcard_146@meta.data$DF.classifications_0.25_0.29_6)

# obcard_178
obcard_178 <- seurat_objects[["obcard_178"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_178, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_178@meta.data))) * (1 - (modelHomotypic(obcard_178@meta.data$seurat_clusters))))
obcard_178 <- doubletFinder(obcard_178, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_178@meta.data)

obcard_178 <- subset(obcard_178, DF.classifications_0.25_0.02_25 != 'Doublet')
seurat_objects[["obcard_178"]] <- obcard_178
doubletest2 <- round((round(0.01 * nrow(obcard_178@meta.data))) * (1 - (modelHomotypic(obcard_178@meta.data$seurat_clusters))))
print(doubletest2)
table(obcard_178@meta.data$DF.classifications_0.25_0.02_25)

# obcard_180
obcard_180 <- seurat_objects[["obcard_180"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_180, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_180@meta.data))) * (1 - (modelHomotypic(obcard_180@meta.data$seurat_clusters))))
obcard_180 <- doubletFinder(obcard_180, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_180@meta.data)

obcard_180 <- subset(obcard_180, DF.classifications_0.25_0.21_13 != 'Doublet')
seurat_objects[["obcard_180"]] <- obcard_180
doubletest3 <- round((round(0.01 * nrow(obcard_180@meta.data))) * (1 - (modelHomotypic(obcard_180@meta.data$seurat_clusters))))
print(doubletest3)
table(obcard_180@meta.data$DF.classifications_0.25_0.21_13)

# obcard_181
obcard_181 <- seurat_objects[["obcard_181"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_181, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_181@meta.data))) * (1 - (modelHomotypic(obcard_181@meta.data$seurat_clusters))))
obcard_181 <- doubletFinder(obcard_181, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_181@meta.data)

obcard_181 <- subset(obcard_181, DF.classifications_0.25_0.01_16 != 'Doublet')
seurat_objects[["obcard_181"]] <- obcard_181
doubletest4 <- round((round(0.01 * nrow(obcard_181@meta.data))) * (1 - (modelHomotypic(obcard_181@meta.data$seurat_clusters))))
print(doubletest4)
table(obcard_181@meta.data$DF.classifications_0.25_0.01_16)

# obcard_183
obcard_183 <- seurat_objects[["obcard_183"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_183, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_183@meta.data))) * (1 - (modelHomotypic(obcard_183@meta.data$seurat_clusters))))
obcard_183 <- doubletFinder(obcard_183, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_183@meta.data)

obcard_183 <- subset(obcard_183, DF.classifications_0.25_0.1_7 != 'Doublet')
seurat_objects[["obcard_183"]] <- obcard_183
doubletest5 <- round((round(0.01 * nrow(obcard_183@meta.data))) * (1 - (modelHomotypic(obcard_183@meta.data$seurat_clusters))))
print(doubletest5)
table(obcard_183@meta.data$DF.classifications_0.25_0.1_7)

# obcard_186 - not enough cells

# obcard_187
obcard_187 <- seurat_objects[["obcard_187"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_187, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_187@meta.data))) * (1 - (modelHomotypic(obcard_187@meta.data$seurat_clusters))))
obcard_187 <- doubletFinder(obcard_187, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_187@meta.data)

obcard_187 <- subset(obcard_187, DF.classifications_0.25_0.12_4 != 'Doublet')
seurat_objects[["obcard_187"]] <- obcard_187
doubletest7 <- round((round(0.01 * nrow(obcard_187@meta.data))) * (1 - (modelHomotypic(obcard_187@meta.data$seurat_clusters))))
print(doubletest7)
table(obcard_187@meta.data$DF.classifications_0.25_0.12_4)

# obcard_188
obcard_188 <- seurat_objects[["obcard_188"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_188, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_188@meta.data))) * (1 - (modelHomotypic(obcard_188@meta.data$seurat_clusters))))
obcard_188 <- doubletFinder(obcard_188, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_188@meta.data)

obcard_188 <- subset(obcard_188, DF.classifications_0.25_0.16_3 != 'Doublet')
seurat_objects[["obcard_188"]] <- obcard_188
doubletest8 <- round((round(0.01 * nrow(obcard_188@meta.data))) * (1 - (modelHomotypic(obcard_188@meta.data$seurat_clusters))))
print(doubletest8)
table(obcard_188@meta.data$DF.classifications_0.25_0.16_3)

# obcard_190 
obcard_190 <- seurat_objects[["obcard_190"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_190, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_190@meta.data))) * (1 - (modelHomotypic(obcard_190@meta.data$seurat_clusters))))
obcard_190 <- doubletFinder(obcard_190, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_190@meta.data)

obcard_190 <- subset(obcard_190, DF.classifications_0.25_0.13_4 != 'Doublet')
seurat_objects[["obcard_190"]] <- obcard_190
doubletest9 <- round((round(0.01 * nrow(obcard_190@meta.data))) * (1 - (modelHomotypic(obcard_190@meta.data$seurat_clusters))))
print(doubletest9)
table(obcard_190@meta.data$DF.classifications_0.25_0.13_4)

# obcard_192
obcard_192 <- seurat_objects[["obcard_192"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_192, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_192@meta.data))) * (1 - (modelHomotypic(obcard_192@meta.data$seurat_clusters))))
obcard_192 <- doubletFinder(obcard_192, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_192@meta.data)

obcard_192 <- subset(obcard_192, DF.classifications_0.25_0.03_8 != 'Doublet')
seurat_objects[["obcard_192"]] <- obcard_192
doubletest10 <- round((round(0.01 * nrow(obcard_192@meta.data))) * (1 - (modelHomotypic(obcard_192@meta.data$seurat_clusters))))
print(doubletest10)
table(obcard_192@meta.data$DF.classifications_0.25_0.03_8)

# obcard_193
obcard_193 <- seurat_objects[["obcard_193"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_193, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_193@meta.data))) * (1 - (modelHomotypic(obcard_193@meta.data$seurat_clusters))))
obcard_193 <- doubletFinder(obcard_193, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_193@meta.data)

obcard_193 <- subset(obcard_193, DF.classifications_0.25_0.25_6 != 'Doublet')
seurat_objects[["obcard_193"]] <- obcard_193
doubletest11 <- round((round(0.01 * nrow(obcard_193@meta.data))) * (1 - (modelHomotypic(obcard_193@meta.data$seurat_clusters))))
print(doubletest11)
table(obcard_193@meta.data$DF.classifications_0.25_0.25_6)

# obcard_194
obcard_194 <- seurat_objects[["obcard_194"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_194, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_194@meta.data))) * (1 - (modelHomotypic(obcard_194@meta.data$seurat_clusters))))
obcard_194 <- doubletFinder(obcard_194, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_194@meta.data)

obcard_194 <- subset(obcard_194, DF.classifications_0.25_0.02_12 != 'Doublet')
seurat_objects[["obcard_194"]] <- obcard_194
doubletest12 <- round((round(0.01 * nrow(obcard_194@meta.data))) * (1 - (modelHomotypic(obcard_194@meta.data$seurat_clusters))))
print(doubletest12)
table(obcard_194@meta.data$DF.classifications_0.25_0.02_12)

# obcard_195
obcard_195 <- seurat_objects[["obcard_195"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_195, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_195@meta.data))) * (1 - (modelHomotypic(obcard_195@meta.data$seurat_clusters))))
obcard_195 <- doubletFinder(obcard_195, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_195@meta.data)

obcard_195 <- subset(obcard_195, DF.classifications_0.25_0.03_4 != 'Doublet')
seurat_objects[["obcard_195"]] <- obcard_195
doubletest13 <- round((round(0.01 * nrow(obcard_195@meta.data))) * (1 - (modelHomotypic(obcard_195@meta.data$seurat_clusters))))
print(doubletest13)
table(obcard_195@meta.data$DF.classifications_0.25_0.03_4)

# obcard_196
obcard_196 <- seurat_objects[["obcard_196"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_196, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_196@meta.data))) * (1 - (modelHomotypic(obcard_196@meta.data$seurat_clusters))))
obcard_196 <- doubletFinder(obcard_196, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_196@meta.data)

obcard_196 <- subset(obcard_196, DF.classifications_0.25_0.15_2 != 'Doublet')
seurat_objects[["obcard_196"]] <- obcard_196
doubletest14 <- round((round(0.01 * nrow(obcard_196@meta.data))) * (1 - (modelHomotypic(obcard_196@meta.data$seurat_clusters))))
print(doubletest14)
table(obcard_196@meta.data$DF.classifications_0.25_0.15_2)

# obcard_197 - not enough cells

# obcard_198 - not enough cells

# Save the list of all samples (after doublet removal)
saveRDS(seurat_objects, file = "/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_2_filtered_noDoublets.RDS")
