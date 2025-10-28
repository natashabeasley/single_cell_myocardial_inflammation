## 1. filtering_18sample_valcard

# valcard multiome
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

samples_df <- read.csv("/home/n/nb443/Documents/Val-CARD/file_locations.csv", header = TRUE)
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
samples_df <- read.csv("/home/n/nb443/Documents/Val-CARD/file_locations.csv",
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
saveRDS(so7, file = "/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_valcard_filtered.RDS")

# DoubletFinder doublet removal
# split the filtered merged object by sample
seurat_objects <- SplitObject(so7, split.by = "sample")
seurat_objects <- readRDS("/scratch/pigblast/nb443/GITHUB_FINALS/16_sample_valcard_filtered.RDS")

# V003
V003 <- seurat_objects[["V003"]]
params <- find.pK(summarizeSweep(paramSweep(V003, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V003@meta.data))) * (1 - (modelHomotypic(V003@meta.data$seurat_clusters))))
V003 <- doubletFinder(V003, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V003@meta.data)

V003 <- subset(V003, DF.classifications_0.25_0.07_5 != 'Doublet')
seurat_objects[["V003"]] <- V003
doubletest1 <- round((round(0.01 * nrow(V003@meta.data))) * (1 - (modelHomotypic(V003@meta.data$seurat_clusters))))
print(doubletest1)
table(V003@meta.data$DF.classifications_0.25_0.07_5)

# V004
V004 <- seurat_objects[["V004"]]
params <- find.pK(summarizeSweep(paramSweep(V004, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V004@meta.data))) * (1 - (modelHomotypic(V004@meta.data$seurat_clusters))))
V004 <- doubletFinder(V004, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V004@meta.data)

V004 <- subset(V004, DF.classifications_0.25_0.13_8 != 'Doublet')
seurat_objects[["V004"]] <- V004
doubletest2 <- round((round(0.01 * nrow(V004@meta.data))) * (1 - (modelHomotypic(V004@meta.data$seurat_clusters))))
print(doubletest2)
table(V004@meta.data$DF.classifications_0.25_0.13_8)

# V007
V007 <- seurat_objects[["V007"]]
params <- find.pK(summarizeSweep(paramSweep(V007, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V007@meta.data))) * (1 - (modelHomotypic(V007@meta.data$seurat_clusters))))
V007 <- doubletFinder(V007, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V007@meta.data)

V007 <- subset(V007, DF.classifications_0.25_0.2_6 != 'Doublet')
seurat_objects[["V007"]] <- V007
doubletest3 <- round((round(0.01 * nrow(V007@meta.data))) * (1 - (modelHomotypic(V007@meta.data$seurat_clusters))))
print(doubletest3)
table(V007@meta.data$DF.classifications_0.25_0.2_6)

# V010
V010 <- seurat_objects[["V010"]]
params <- find.pK(summarizeSweep(paramSweep(V010, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V010@meta.data))) * (1 - (modelHomotypic(V010@meta.data$seurat_clusters))))
V010 <- doubletFinder(V010, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V010@meta.data)

V010 <- subset(V010, DF.classifications_0.25_0.28_1 != 'Doublet')
seurat_objects[["V010"]] <- V010
doubletest4 <- round((round(0.01 * nrow(V010@meta.data))) * (1 - (modelHomotypic(V010@meta.data$seurat_clusters))))
print(doubletest4)
table(V010@meta.data$DF.classifications_0.25_0.28_1)

# V011
V011 <- seurat_objects[["V011"]]
params <- find.pK(summarizeSweep(paramSweep(V011, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V011@meta.data))) * (1 - (modelHomotypic(V011@meta.data$seurat_clusters))))
V011 <- doubletFinder(V011, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V011@meta.data)

V011 <- subset(V011, DF.classifications_0.25_0.03_6 != 'Doublet')
seurat_objects[["V011"]] <- V011
doubletest5 <- round((round(0.01 * nrow(V011@meta.data))) * (1 - (modelHomotypic(V011@meta.data$seurat_clusters))))
print(doubletest5)
table(V011@meta.data$DF.classifications_0.25_0.03_6)

# V014
V014 <- seurat_objects[["V014"]]
params <- find.pK(summarizeSweep(paramSweep(V014, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V014@meta.data))) * (1 - (modelHomotypic(V014@meta.data$seurat_clusters))))
V014 <- doubletFinder(V014, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V014@meta.data)

V014 <- subset(V014, DF.classifications_0.25_0.16_7 != 'Doublet')
seurat_objects[["V014"]] <- V014
doubletest6 <- round((round(0.01 * nrow(V014@meta.data))) * (1 - (modelHomotypic(V014@meta.data$seurat_clusters))))
print(doubletest6)
table(V014@meta.data$DF.classifications_0.25_0.16_7)

# V017
V017 <- seurat_objects[["V017"]]
params <- find.pK(summarizeSweep(paramSweep(V017, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V017@meta.data))) * (1 - (modelHomotypic(V017@meta.data$seurat_clusters))))
V017 <- doubletFinder(V017, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V017@meta.data)

V017 <- subset(V017, DF.classifications_0.25_0.07_4 != 'Doublet')
seurat_objects[["V017"]] <- V017
doubletest7 <- round((round(0.01 * nrow(V017@meta.data))) * (1 - (modelHomotypic(V017@meta.data$seurat_clusters))))
print(doubletest7)
table(V017@meta.data$DF.classifications_0.25_0.07_4)

# V026
V026 <- seurat_objects[["V026"]]
params <- find.pK(summarizeSweep(paramSweep(V026, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V026@meta.data))) * (1 - (modelHomotypic(V026@meta.data$seurat_clusters))))
V026 <- doubletFinder(V026, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V026@meta.data)

V026 <- subset(V026, DF.classifications_0.25_0.26_2 != 'Doublet')
seurat_objects[["V026"]] <- V026
doubletest8 <- round((round(0.01 * nrow(V026@meta.data))) * (1 - (modelHomotypic(V026@meta.data$seurat_clusters))))
print(doubletest8)
table(V026@meta.data$DF.classifications_0.25_0.26_2)

# V027 
V027 <- seurat_objects[["V027"]]
params <- find.pK(summarizeSweep(paramSweep(V027, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V027@meta.data))) * (1 - (modelHomotypic(V027@meta.data$seurat_clusters))))
V027 <- doubletFinder(V027, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V027@meta.data)

V027 <- subset(V027, DF.classifications_0.25_0.03_7 != 'Doublet')
seurat_objects[["V027"]] <- V027
doubletest9 <- round((round(0.01 * nrow(V027@meta.data))) * (1 - (modelHomotypic(V027@meta.data$seurat_clusters))))
print(doubletest9)
table(V027@meta.data$DF.classifications_0.25_0.03_7)

# V028
V028 <- seurat_objects[["V028"]]
params <- find.pK(summarizeSweep(paramSweep(V028, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V028@meta.data))) * (1 - (modelHomotypic(V028@meta.data$seurat_clusters))))
V028 <- doubletFinder(V028, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V028@meta.data)

V028 <- subset(V028, DF.classifications_0.25_0.03_5 != 'Doublet')
seurat_objects[["V028"]] <- V028
doubletest10 <- round((round(0.01 * nrow(V028@meta.data))) * (1 - (modelHomotypic(V028@meta.data$seurat_clusters))))
print(doubletest10)
table(V028@meta.data$DF.classifications_0.25_0.03_5)

# V030
V030 <- seurat_objects[["V030"]]
params <- find.pK(summarizeSweep(paramSweep(V030, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V030@meta.data))) * (1 - (modelHomotypic(V030@meta.data$seurat_clusters))))
V030 <- doubletFinder(V030, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V030@meta.data)

V030 <- subset(V030, DF.classifications_0.25_0.11_5 != 'Doublet')
seurat_objects[["V030"]] <- V030
doubletest11 <- round((round(0.01 * nrow(V030@meta.data))) * (1 - (modelHomotypic(V030@meta.data$seurat_clusters))))
print(doubletest11)
table(V030@meta.data$DF.classifications_0.25_0.11_5)

# V031
V031 <- seurat_objects[["V031"]]
params <- find.pK(summarizeSweep(paramSweep(V031, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V031@meta.data))) * (1 - (modelHomotypic(V031@meta.data$seurat_clusters))))
V031 <- doubletFinder(V031, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V031@meta.data)

V031 <- subset(V031, DF.classifications_0.25_0.06_2 != 'Doublet')
seurat_objects[["V031"]] <- V031
doubletest12 <- round((round(0.01 * nrow(V031@meta.data))) * (1 - (modelHomotypic(V031@meta.data$seurat_clusters))))
print(doubletest12)
table(V031@meta.data$DF.classifications_0.25_0.06_2)

# V032
V032 <- seurat_objects[["V032"]]
params <- find.pK(summarizeSweep(paramSweep(V032, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V032@meta.data))) * (1 - (modelHomotypic(V032@meta.data$seurat_clusters))))
V032 <- doubletFinder(V032, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V032@meta.data)

V032 <- subset(V032, DF.classifications_0.25_0.27_2 != 'Doublet')
seurat_objects[["V032"]] <- V032
doubletest13 <- round((round(0.01 * nrow(V032@meta.data))) * (1 - (modelHomotypic(V032@meta.data$seurat_clusters))))
print(doubletest13)
table(V032@meta.data$DF.classifications_0.25_0.27_2)

# V034
V034 <- seurat_objects[["V034"]]
params <- find.pK(summarizeSweep(paramSweep(V034, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V034@meta.data))) * (1 - (modelHomotypic(V034@meta.data$seurat_clusters))))
V034 <- doubletFinder(V034, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V034@meta.data)

V034 <- subset(V034, DF.classifications_0.25_0.26_3 != 'Doublet')
seurat_objects[["V034"]] <- V034
doubletest14 <- round((round(0.01 * nrow(V034@meta.data))) * (1 - (modelHomotypic(V034@meta.data$seurat_clusters))))
print(doubletest14)
table(V034@meta.data$DF.classifications_0.25_0.26_3)

# V035
V035 <- seurat_objects[["V035"]]
params <- find.pK(summarizeSweep(paramSweep(V035, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V035@meta.data))) * (1 - (modelHomotypic(V035@meta.data$seurat_clusters))))
V035 <- doubletFinder(V035, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V035@meta.data)

V035 <- subset(V035, DF.classifications_0.25_0.25_1 != 'Doublet')
seurat_objects[["V035"]] <- V035
doubletest15 <- round((round(0.01 * nrow(V035@meta.data))) * (1 - (modelHomotypic(V035@meta.data$seurat_clusters))))
print(doubletest15)
table(V035@meta.data$DF.classifications_0.25_0.25_1)

# V041
V041 <- seurat_objects[["V041"]]
params <- find.pK(summarizeSweep(paramSweep(V041, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V041@meta.data))) * (1 - (modelHomotypic(V041@meta.data$seurat_clusters))))
V041 <- doubletFinder(V041, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V041@meta.data)

V041 <- subset(V041, DF.classifications_0.25_0.11_1 != 'Doublet')
seurat_objects[["V041"]] <- V041
doubletest16 <- round((round(0.01 * nrow(V041@meta.data))) * (1 - (modelHomotypic(V041@meta.data$seurat_clusters))))
print(doubletest16)
table(V041@meta.data$DF.classifications_0.25_0.11_1)

# V043
V043 <- seurat_objects[["V043"]]
params <- find.pK(summarizeSweep(paramSweep(V043, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V043@meta.data))) * (1 - (modelHomotypic(V043@meta.data$seurat_clusters))))
V043 <- doubletFinder(V043, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V043@meta.data)

V043 <- subset(V043, DF.classifications_0.25_0.2_3 != 'Doublet')
seurat_objects[["V043"]] <- V043
doubletest17 <- round((round(0.01 * nrow(V043@meta.data))) * (1 - (modelHomotypic(V043@meta.data$seurat_clusters))))
print(doubletest17)
table(V043@meta.data$DF.classifications_0.25_0.2_3)

# V045
V045 <- seurat_objects[["V045"]]
params <- find.pK(summarizeSweep(paramSweep(V045, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(V045@meta.data))) * (1 - (modelHomotypic(V045@meta.data$seurat_clusters))))
V045 <- doubletFinder(V045, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(V045@meta.data)

V045 <- subset(V045, DF.classifications_0.25_0.11_1 != 'Doublet')
seurat_objects[["V045"]] <- V045
doubletest18 <- round((round(0.01 * nrow(V043@meta.data))) * (1 - (modelHomotypic(V045@meta.data$seurat_clusters))))
print(doubletest18)
table(V045@meta.data$DF.classifications_0.25_0.11_1)

# Save the list of all samples (after doublet removal)
saveRDS(seurat_objects, file = "/scratch/pigblast/nb443/GITHUB_FINALS/18_sample_valcard_filtered_noDoublets.RDS")
