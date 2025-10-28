## 1. filtering_16sample_1st_run

# 1st run 16 sample multiome (Ob-CARD)
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
packageVersion("DoubletFinder") 

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'  # make sure chromosome names match your data
genome(annotations) <- "hg38"

samples_df <- read.csv("/home/n/nb443/Documents/Ob_CARD2_snRNA-seq_analysis/16_sample_1st_run_locations", header = TRUE)
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
samples_df <- read.csv("/home/n/nb443/Documents/Ob_CARD2_snRNA-seq_analysis/16_sample_1st_run_locations",
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
    .groups = "drop"
  ) %>%
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
saveRDS(so7, file = "/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_1_filtered.RDS")

# DoubletFinder doublet removal
# split the filtered merged object by sample
seurat_objects <- SplitObject(so7, split.by = "sample")
seurat_objects <- readRDS("/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_1_filtered.RDS")

# obcard_92
obcard_92 <- seurat_objects[["obcard_92"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_92, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_92@meta.data))) * (1 - (modelHomotypic(obcard_92@meta.data$seurat_clusters))))
obcard_92 <- doubletFinder(obcard_92, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_92@meta.data)

obcard_92 <- subset(obcard_92, DF.classifications_0.25_0.1_4 != 'Doublet')
seurat_objects[["obcard_92"]] <- obcard_92
doubletest1 <- round((round(0.01 * nrow(obcard_92@meta.data))) * (1 - (modelHomotypic(obcard_92@meta.data$seurat_clusters))))
print(doubletest1)
table(obcard_92@meta.data$DF.classifications_0.25_0.1_4)

# obcard_113
obcard_113 <- seurat_objects[["obcard_113"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_113, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_113@meta.data))) * (1 - (modelHomotypic(obcard_113@meta.data$seurat_clusters))))
obcard_113 <- doubletFinder(obcard_113, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_113@meta.data)

obcard_113 <- subset(obcard_113, DF.classifications_0.25_0.25_1 != 'Doublet')
seurat_objects[["obcard_113"]] <- obcard_113
doubletest2 <- round((round(0.01 * nrow(obcard_113@meta.data))) * (1 - (modelHomotypic(obcard_113@meta.data$seurat_clusters))))
print(doubletest2)
table(obcard_113@meta.data$DF.classifications_0.25_0.25_1)

# obcard_140
obcard_140 <- seurat_objects[["obcard_140"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_140, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_140@meta.data))) * (1 - (modelHomotypic(obcard_140@meta.data$seurat_clusters))))
obcard_140 <- doubletFinder(obcard_140, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_140@meta.data)

obcard_140 <- subset(obcard_140, DF.classifications_0.25_0.3_3 != 'Doublet')
seurat_objects[["obcard_140"]] <- obcard_140
doubletest3 <- round((round(0.01 * nrow(obcard_140@meta.data))) * (1 - (modelHomotypic(obcard_140@meta.data$seurat_clusters))))
print(doubletest3)
table(obcard_140@meta.data$DF.classifications_0.25_0.3_3)

# obcard_147
obcard_147 <- seurat_objects[["obcard_147"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_147, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_147@meta.data))) * (1 - (modelHomotypic(obcard_147@meta.data$seurat_clusters))))
obcard_147 <- doubletFinder(obcard_147, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_147@meta.data)

obcard_147 <- subset(obcard_147, DF.classifications_0.25_0.02_10 != 'Doublet')
seurat_objects[["obcard_147"]] <- obcard_147
doubletest4 <- round((round(0.01 * nrow(obcard_147@meta.data))) * (1 - (modelHomotypic(obcard_147@meta.data$seurat_clusters))))
print(doubletest4)
table(obcard_147@meta.data$DF.classifications_0.25_0.02_10)

# obcard_148
obcard_148 <- seurat_objects[["obcard_148"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_148, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_148@meta.data))) * (1 - (modelHomotypic(obcard_148@meta.data$seurat_clusters))))
obcard_148 <- doubletFinder(obcard_148, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_148@meta.data)

obcard_148 <- subset(obcard_148, DF.classifications_0.25_0.03_6 != 'Doublet')
seurat_objects[["obcard_148"]] <- obcard_148
doubletest5 <- round((round(0.01 * nrow(obcard_148@meta.data))) * (1 - (modelHomotypic(obcard_148@meta.data$seurat_clusters))))
print(doubletest5)
table(obcard_148@meta.data$DF.classifications_0.25_0.03_6)

# obcard_151
obcard_151 <- seurat_objects[["obcard_151"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_151, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_151@meta.data))) * (1 - (modelHomotypic(obcard_151@meta.data$seurat_clusters))))
obcard_151 <- doubletFinder(obcard_151, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_151@meta.data)

obcard_151 <- subset(obcard_151, DF.classifications_0.25_0.12_13 != 'Doublet')
seurat_objects[["obcard_151"]] <- obcard_151
doubletest6 <- round((round(0.01 * nrow(obcard_151@meta.data))) * (1 - (modelHomotypic(obcard_151@meta.data$seurat_clusters))))
print(doubletest6)
table(obcard_151@meta.data$DF.classifications_0.25_0.12_13)

# obcard_154
obcard_154 <- seurat_objects[["obcard_154"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_154, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_154@meta.data))) * (1 - (modelHomotypic(obcard_154@meta.data$seurat_clusters))))
obcard_154 <- doubletFinder(obcard_154, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_154@meta.data)

obcard_154 <- subset(obcard_154, DF.classifications_0.25_0.13_1 != 'Doublet')
seurat_objects[["obcard_154"]] <- obcard_154
doubletest7 <- round((round(0.01 * nrow(obcard_154@meta.data))) * (1 - (modelHomotypic(obcard_154@meta.data$seurat_clusters))))
print(doubletest7)
table(obcard_154@meta.data$DF.classifications_0.25_0.13_1)

# obcard_156
obcard_156 <- seurat_objects[["obcard_156"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_156, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_156@meta.data))) * (1 - (modelHomotypic(obcard_156@meta.data$seurat_clusters))))
obcard_156 <- doubletFinder(obcard_156, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_156@meta.data)

obcard_156 <- subset(obcard_156, DF.classifications_0.25_0.14_13 != 'Doublet')
seurat_objects[["obcard_156"]] <- obcard_156
doubletest8 <- round((round(0.01 * nrow(obcard_156@meta.data))) * (1 - (modelHomotypic(obcard_156@meta.data$seurat_clusters))))
print(doubletest8)
table(obcard_156@meta.data$DF.classifications_0.25_0.14_13)

# obcard_157 - unable to perform 

# obcard_159
obcard_159 <- seurat_objects[["obcard_159"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_159, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_159@meta.data))) * (1 - (modelHomotypic(obcard_159@meta.data$seurat_clusters))))
obcard_159 <- doubletFinder(obcard_159, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_159@meta.data)

obcard_159 <- subset(obcard_159, DF.classifications_0.25_0.01_14 != 'Doublet')
seurat_objects[["obcard_159"]] <- obcard_159
doubletest10 <- round((round(0.01 * nrow(obcard_159@meta.data))) * (1 - (modelHomotypic(obcard_159@meta.data$seurat_clusters))))
print(doubletest10)
table(obcard_159@meta.data$DF.classifications_0.25_0.01_14)

# obcard_166
obcard_166 <- seurat_objects[["obcard_166"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_166, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_166@meta.data))) * (1 - (modelHomotypic(obcard_166@meta.data$seurat_clusters))))
obcard_166 <- doubletFinder(obcard_166, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_166@meta.data)

obcard_166 <- subset(obcard_166, DF.classifications_0.25_0.27_0 != 'Doublet')
seurat_objects[["obcard_166"]] <- obcard_166
doubletest11 <- round((round(0.01 * nrow(obcard_166@meta.data))) * (1 - (modelHomotypic(obcard_166@meta.data$seurat_clusters))))
print(doubletest11)
table(obcard_166@meta.data$DF.classifications_0.25_0.27_0)

# obcard_167 - not enough cells

# obcard_168
obcard_168 <- seurat_objects[["obcard_168"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_168, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_168@meta.data))) * (1 - (modelHomotypic(obcard_168@meta.data$seurat_clusters))))
obcard_168 <- doubletFinder(obcard_168, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_168@meta.data)

obcard_168 <- subset(obcard_168, DF.classifications_0.25_0.27_0 != 'Doublet')
seurat_objects[["obcard_168"]] <- obcard_168
doubletest12 <- round((round(0.01 * nrow(obcard_168@meta.data))) * (1 - (modelHomotypic(obcard_168@meta.data$seurat_clusters))))
print(doubletest12)
table(obcard_168@meta.data$DF.classifications_0.25_0.27_0)

# obcard_171
obcard_171 <- seurat_objects[["obcard_171"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_171, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_171@meta.data))) * (1 - (modelHomotypic(obcard_171@meta.data$seurat_clusters))))
obcard_171 <- doubletFinder(obcard_171, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_171@meta.data)

obcard_171 <- subset(obcard_171, DF.classifications_0.25_0.1_5 != 'Doublet')
seurat_objects[["obcard_171"]] <- obcard_171
doubletest13 <- round((round(0.01 * nrow(obcard_171@meta.data))) * (1 - (modelHomotypic(obcard_171@meta.data$seurat_clusters))))
print(doubletest13)
table(obcard_171@meta.data$DF.classifications_0.25_0.1_5)

# obcard_175
obcard_175 <- seurat_objects[["obcard_175"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_175, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_175@meta.data))) * (1 - (modelHomotypic(obcard_175@meta.data$seurat_clusters))))
obcard_175 <- doubletFinder(obcard_175, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_175@meta.data)

obcard_175 <- subset(obcard_175, DF.classifications_0.25_0.24_36 != 'Doublet')
seurat_objects[["obcard_175"]] <- obcard_175
doubletest14 <- round((round(0.01 * nrow(obcard_175@meta.data))) * (1 - (modelHomotypic(obcard_175@meta.data$seurat_clusters))))
print(doubletest14)
table(obcard_175@meta.data$DF.classifications_0.25_0.24_36)

# obcard_176
obcard_176 <- seurat_objects[["obcard_176"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_176, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_176@meta.data))) * (1 - (modelHomotypic(obcard_176@meta.data$seurat_clusters))))
obcard_176 <- doubletFinder(obcard_176, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_176@meta.data)

obcard_176 <- subset(obcard_176, DF.classifications_0.25_0.09_5 != 'Doublet')
seurat_objects[["obcard_176"]] <- obcard_176
doubletest15 <- round((round(0.01 * nrow(obcard_176@meta.data))) * (1 - (modelHomotypic(obcard_176@meta.data$seurat_clusters))))
print(doubletest15)
table(obcard_176@meta.data$DF.classifications_0.25_0.09_5)

# Save the list of all samples (after doublet removal)
saveRDS(seurat_objects, file = "/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_1_filtered_noDoublets.RDS")

