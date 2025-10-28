# 1. filtering_mito+feat_21_samples
## 2. doublet_removal_21_samples

set.seed(123)
library(DoubletFinder) # Version 2.0.6
packageVersion("DoubletFinder") 

seurat_objects <- readRDS("/scratch/pigblast/nb443/GITHUB_FINALS/21_mito+feat_filtered.rds")
seurat_objects <- SplitObject(seurat_objects, split.by = "sample")

# obcard_69
obcard_69 <- seurat_objects[["obcard_69"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_69, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_69@meta.data))) * (1 - (modelHomotypic(obcard_69@meta.data$seurat_clusters))))
obcard_69 <- doubletFinder(obcard_69, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_69@meta.data)

obcard_69 <- subset(obcard_69, DF.classifications_0.25_0.05_15 != 'Doublet')
seurat_objects[["obcard_69"]] <- obcard_69
doubletest1 <- round((round(0.01 * nrow(obcard_69@meta.data))) * (1 - (modelHomotypic(obcard_69@meta.data$seurat_clusters))))
print(doubletest1)
table(obcard_69@meta.data$DF.classifications_0.25_0.05_15)

# obcard_74
obcard_74 <- seurat_objects[["obcard_74"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_74, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_74@meta.data))) * (1 - (modelHomotypic(obcard_74@meta.data$seurat_clusters))))
obcard_74 <- doubletFinder(obcard_74, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_74@meta.data)

obcard_74 <- subset(obcard_74, DF.classifications_0.25_0.06_21 != 'Doublet')
seurat_objects[["obcard_74"]] <- obcard_74
doubletest2 <- round((round(0.01 * nrow(obcard_74@meta.data))) * (1 - (modelHomotypic(obcard_74@meta.data$seurat_clusters))))
print(doubletest2)
table(obcard_74@meta.data$DF.classifications_0.25_0.06_21)

# obcard_77
obcard_77 <- seurat_objects[["obcard_77"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_77, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_77@meta.data))) * (1 - (modelHomotypic(obcard_77@meta.data$seurat_clusters))))
obcard_77 <- doubletFinder(obcard_77, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_77@meta.data)

obcard_77 <- subset(obcard_77, DF.classifications_0.25_0.02_8 != 'Doublet')
seurat_objects[["obcard_77"]] <- obcard_77
doubletest3 <- round((round(0.01 * nrow(obcard_77@meta.data))) * (1 - (modelHomotypic(obcard_77@meta.data$seurat_clusters))))
print(doubletest3)
table(obcard_77@meta.data$DF.classifications_0.25_0.02_8)

# obcard_85
obcard_85 <- seurat_objects[["obcard_85"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_85, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_85@meta.data))) * (1 - (modelHomotypic(obcard_85@meta.data$seurat_clusters))))
obcard_85 <- doubletFinder(obcard_85, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_85@meta.data)

obcard_85 <- subset(obcard_85, DF.classifications_0.25_0.28_2 != 'Doublet')
seurat_objects[["obcard_85"]] <- obcard_85
doubletest4 <- round((round(0.01 * nrow(obcard_85@meta.data))) * (1 - (modelHomotypic(obcard_85@meta.data$seurat_clusters))))
print(doubletest4)
table(obcard_85@meta.data$DF.classifications_0.25_0.28_2)

# obcard_88
obcard_88 <- seurat_objects[["obcard_88"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_88, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_88@meta.data))) * (1 - (modelHomotypic(obcard_88@meta.data$seurat_clusters))))
obcard_88 <- doubletFinder(obcard_88, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_88@meta.data)

obcard_88 <- subset(obcard_88, DF.classifications_0.25_0.14_4 != 'Doublet')
seurat_objects[["obcard_88"]] <- obcard_88
doubletest5 <- round((round(0.01 * nrow(obcard_88@meta.data))) * (1 - (modelHomotypic(obcard_88@meta.data$seurat_clusters))))
print(doubletest5)
table(obcard_88@meta.data$DF.classifications_0.25_0.14_4)

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

obcard_92 <- subset(obcard_92, DF.classifications_0.25_0.25_3 != 'Doublet')
seurat_objects[["obcard_92"]] <- obcard_92
doubletest6 <- round((round(0.01 * nrow(obcard_92@meta.data))) * (1 - (modelHomotypic(obcard_92@meta.data$seurat_clusters))))
print(doubletest6)
table(obcard_92@meta.data$DF.classifications_0.25_0.25_3)

# obcard_96
obcard_96 <- seurat_objects[["obcard_96"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_96, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_96@meta.data))) * (1 - (modelHomotypic(obcard_96@meta.data$seurat_clusters))))
obcard_96 <- doubletFinder(obcard_96, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_96@meta.data)

obcard_96 <- subset(obcard_96, DF.classifications_0.25_0.04_4 != 'Doublet')
seurat_objects[["obcard_96"]] <- obcard_96
doubletest7 <- round((round(0.01 * nrow(obcard_96@meta.data))) * (1 - (modelHomotypic(obcard_96@meta.data$seurat_clusters))))
print(doubletest7)
table(obcard_96@meta.data$DF.classifications_0.25_0.04_4)

# obcard_97
obcard_97 <- seurat_objects[["obcard_97"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_97, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_97@meta.data))) * (1 - (modelHomotypic(obcard_97@meta.data$seurat_clusters))))
obcard_97 <- doubletFinder(obcard_97, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_97@meta.data)

obcard_97 <- subset(obcard_97, DF.classifications_0.25_0.23_3 != 'Doublet')
seurat_objects[["obcard_97"]] <- obcard_97
doubletest8 <- round((round(0.01 * nrow(obcard_97@meta.data))) * (1 - (modelHomotypic(obcard_97@meta.data$seurat_clusters))))
print(doubletest8)
table(obcard_97@meta.data$DF.classifications_0.25_0.23_3)

# obcard_126
obcard_126 <- seurat_objects[["obcard_126"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_126, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_126@meta.data))) * (1 - (modelHomotypic(obcard_126@meta.data$seurat_clusters))))
obcard_126 <- doubletFinder(obcard_126, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_126@meta.data)

obcard_126 <- subset(obcard_126, DF.classifications_0.25_0.28_12 != 'Doublet')
seurat_objects[["obcard_126"]] <- obcard_126
doubletest9 <- round((round(0.01 * nrow(obcard_126@meta.data))) * (1 - (modelHomotypic(obcard_126@meta.data$seurat_clusters))))
print(doubletest9)
table(obcard_126@meta.data$DF.classifications_0.25_0.28_12)

# obcard_129
obcard_129 <- seurat_objects[["obcard_129"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_129, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_129@meta.data))) * (1 - (modelHomotypic(obcard_129@meta.data$seurat_clusters))))
obcard_129 <- doubletFinder(obcard_129, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_129@meta.data)

obcard_129 <- subset(obcard_129, DF.classifications_0.25_0.03_18 != 'Doublet')
seurat_objects[["obcard_129"]] <- obcard_129
doubletest10 <- round((round(0.01 * nrow(obcard_129@meta.data))) * (1 - (modelHomotypic(obcard_129@meta.data$seurat_clusters))))
print(doubletest10)
table(obcard_129@meta.data$DF.classifications_0.25_0.03_18)

# obcard_132
obcard_132 <- seurat_objects[["obcard_132"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_132, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_132@meta.data))) * (1 - (modelHomotypic(obcard_132@meta.data$seurat_clusters))))
obcard_132 <- doubletFinder(obcard_132, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_132@meta.data)

obcard_132 <- subset(obcard_132, DF.classifications_0.25_0.12_9 != 'Doublet')
seurat_objects[["obcard_132"]] <- obcard_132
doubletest11 <- round((round(0.01 * nrow(obcard_132@meta.data))) * (1 - (modelHomotypic(obcard_132@meta.data$seurat_clusters))))
print(doubletest11)
table(obcard_132@meta.data$DF.classifications_0.25_0.12_9)

# obcard_134
obcard_134 <- seurat_objects[["obcard_134"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_134, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_134@meta.data))) * (1 - (modelHomotypic(obcard_134@meta.data$seurat_clusters))))
obcard_134 <- doubletFinder(obcard_134, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_134@meta.data)

obcard_134 <- subset(obcard_134, DF.classifications_0.25_0.16_6 != 'Doublet')
seurat_objects[["obcard_134"]] <- obcard_134
doubletest12 <- round((round(0.01 * nrow(obcard_134@meta.data))) * (1 - (modelHomotypic(obcard_134@meta.data$seurat_clusters))))
print(doubletest12)
table(obcard_134@meta.data$DF.classifications_0.25_0.16_6)

# obcard_135
obcard_135 <- seurat_objects[["obcard_135"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_135, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_135@meta.data))) * (1 - (modelHomotypic(obcard_135@meta.data$seurat_clusters))))
obcard_135 <- doubletFinder(obcard_135, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_135@meta.data)

obcard_135 <- subset(obcard_135, DF.classifications_0.25_0.04_3 != 'Doublet')
seurat_objects[["obcard_135"]] <- obcard_135
doubletest13 <- round((round(0.01 * nrow(obcard_135@meta.data))) * (1 - (modelHomotypic(obcard_135@meta.data$seurat_clusters))))
print(doubletest13)
table(obcard_135@meta.data$DF.classifications_0.25_0.04_3)

# obcard_136
obcard_136 <- seurat_objects[["obcard_136"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_136, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_136@meta.data))) * (1 - (modelHomotypic(obcard_136@meta.data$seurat_clusters))))
obcard_136 <- doubletFinder(obcard_136, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_136@meta.data)

obcard_136 <- subset(obcard_136, DF.classifications_0.25_0.06_6 != 'Doublet')
seurat_objects[["obcard_136"]] <- obcard_136
doubletest14 <- round((round(0.01 * nrow(obcard_136@meta.data))) * (1 - (modelHomotypic(obcard_136@meta.data$seurat_clusters))))
print(doubletest14)
table(obcard_136@meta.data$DF.classifications_0.25_0.06_6)

# obcard_137
obcard_137 <- seurat_objects[["obcard_137"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_137, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_137@meta.data))) * (1 - (modelHomotypic(obcard_137@meta.data$seurat_clusters))))
obcard_137 <- doubletFinder(obcard_137, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_137@meta.data)

obcard_137 <- subset(obcard_137, DF.classifications_0.25_0.16_3 != 'Doublet')
seurat_objects[["obcard_137"]] <- obcard_137
doubletest15 <- round((round(0.01 * nrow(obcard_137@meta.data))) * (1 - (modelHomotypic(obcard_137@meta.data$seurat_clusters))))
print(doubletest15)
table(obcard_137@meta.data$DF.classifications_0.25_0.16_3)

# obcard_138
obcard_138 <- seurat_objects[["obcard_138"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_138, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_138@meta.data))) * (1 - (modelHomotypic(obcard_138@meta.data$seurat_clusters))))
obcard_138 <- doubletFinder(obcard_138, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_138@meta.data)

obcard_138 <- subset(obcard_138, DF.classifications_0.25_0.19_2 != 'Doublet')
seurat_objects[["obcard_138"]] <- obcard_138
doubletest16 <- round((round(0.01 * nrow(obcard_138@meta.data))) * (1 - (modelHomotypic(obcard_138@meta.data$seurat_clusters))))
print(doubletest16)
table(obcard_138@meta.data$DF.classifications_0.25_0.19_2)

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

obcard_140 <- subset(obcard_140, DF.classifications_0.25_0.3_5 != 'Doublet')
seurat_objects[["obcard_140"]] <- obcard_140
doubletest17 <- round((round(0.01 * nrow(obcard_140@meta.data))) * (1 - (modelHomotypic(obcard_140@meta.data$seurat_clusters))))
print(doubletest17)
table(obcard_140@meta.data$DF.classifications_0.25_0.3_5)

# obcard_141
obcard_141 <- seurat_objects[["obcard_141"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_141, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_141@meta.data))) * (1 - (modelHomotypic(obcard_141@meta.data$seurat_clusters))))
obcard_141 <- doubletFinder(obcard_141, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_141@meta.data)

obcard_141 <- subset(obcard_141, DF.classifications_0.25_0.23_3 != 'Doublet')
seurat_objects[["obcard_141"]] <- obcard_141
doubletest18 <- round((round(0.01 * nrow(obcard_141@meta.data))) * (1 - (modelHomotypic(obcard_141@meta.data$seurat_clusters))))
print(doubletest18)
table(obcard_141@meta.data$DF.classifications_0.25_0.23_3)

# obcard_142
obcard_142 <- seurat_objects[["obcard_142"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_142, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_142@meta.data))) * (1 - (modelHomotypic(obcard_142@meta.data$seurat_clusters))))
obcard_142 <- doubletFinder(obcard_142, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_142@meta.data)

obcard_142 <- subset(obcard_142, DF.classifications_0.25_0.05_4 != 'Doublet')
seurat_objects[["obcard_142"]] <- obcard_142
doubletest19 <- round((round(0.01 * nrow(obcard_142@meta.data))) * (1 - (modelHomotypic(obcard_142@meta.data$seurat_clusters))))
print(doubletest19)
table(obcard_142@meta.data$DF.classifications_0.25_0.05_4)

# obcard_143
obcard_143 <- seurat_objects[["obcard_143"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_143, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_143@meta.data))) * (1 - (modelHomotypic(obcard_143@meta.data$seurat_clusters))))
obcard_143 <- doubletFinder(obcard_143, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_143@meta.data)

obcard_143 <- subset(obcard_143, DF.classifications_0.25_0.25_2 != 'Doublet')
seurat_objects[["obcard_143"]] <- obcard_143
doubletest20 <- round((round(0.01 * nrow(obcard_143@meta.data))) * (1 - (modelHomotypic(obcard_143@meta.data$seurat_clusters))))
print(doubletest20)
table(obcard_143@meta.data$DF.classifications_0.25_0.25_2)

# obcard_144
obcard_144 <- seurat_objects[["obcard_144"]]
params <- find.pK(summarizeSweep(paramSweep(obcard_144, sct = TRUE), GT = FALSE))
pK <- as.numeric(as.character(params$pK))  # Optimal pK value
BCmetric <- params$BCmetric
peak <- pK[which(BCmetric %in% max(BCmetric))]  # Find the optimal pK based on the BCmetric
plot(x = pK, y = BCmetric, pch = 16, type = "b", col = "#41b6c4", lty = 1, cex = 0.8)
abline(v = peak, lwd = 2, col = 'brown', lty = 2)
text(peak, max(BCmetric), as.character(peak), pos = 4, col = "red")
doubletest <- round((round(0.01 * nrow(obcard_144@meta.data))) * (1 - (modelHomotypic(obcard_144@meta.data$seurat_clusters))))
obcard_144 <- doubletFinder(obcard_144, pK = peak, nExp = doubletest, reuse.pANN = NULL, sct = TRUE, PCs = 1:10)
colnames(obcard_144@meta.data)

obcard_144 <- subset(obcard_144, DF.classifications_0.25_0.06_2 != 'Doublet')
seurat_objects[["obcard_144"]] <- obcard_144
doubletest21 <- round((round(0.01 * nrow(obcard_144@meta.data))) * (1 - (modelHomotypic(obcard_144@meta.data$seurat_clusters))))
print(doubletest21)
table(obcard_144@meta.data$DF.classifications_0.25_0.06_2)


#  Save object with doublets removed
saveRDS(seurat_objects, file = "/scratch/pigblast/nb443/GITHUB_FINALS/21_mito+feat_filtered_doublets_removed.rds")
