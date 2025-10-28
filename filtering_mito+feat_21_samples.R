## 1. filtering_mito+feat_21_samples
# 2. doublet_removal_21_samples

### 21 sample filtering, with figure + details ###
library(Seurat) # Version 5.2.1
library(remotes) 
install_version("ggplot2", version = "3.5.2", repos = "https://cran.r-project.org") # for violin plots (BUG)
install_version("shadowtext", version = "0.1.4", repos = "http://cran.r-project.org") # for violin plots (BUG)
library(ggplot2) 
library(shadowtext)
library(dplyr)
library(DoubletFinder)
set.seed(123)
samples <- read.csv("/home/n/nb443/Documents/Ob_CARD1_snRNA-seq_analysis/Scripts/sophia_samples.csv",
                    stringsAsFactors = FALSE)
# Build named list of count matrices with checks
samples_data <- setNames(vector("list", nrow(samples)), samples$name)
for (i in seq_len(nrow(samples))) {
  sample_name <- samples$name[i]
  data_dir    <- samples$location[i]
  if (!dir.exists(data_dir)) {
    warning(sprintf("Directory not found: %s (sample: %s)", data_dir, sample_name))
    next}
  samples_data[[sample_name]] <- Read10X(data.dir = data_dir)}
# Create Seurat objects and add metadata/QC
seurat_objects <- lapply(names(samples_data), function(sample_name) {
  mat <- samples_data[[sample_name]]
  stopifnot(!is.null(mat))
  so <- CreateSeuratObject(counts = mat, project = sample_name)
  so$sample <- sample_name
  # adjust "^MT-" vs "^mt-" depending on organism
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so})
names(seurat_objects) <- names(samples_data)
# Look at one sample by name
print(seurat_objects[[1]])
all_merged <- merge(seurat_objects[[1]], seurat_objects[-1])
table(all_merged$sample)  # counts per sample from the merged object

### Mitochondrial content ###
# Violin plots of QC metrics per sample
# Extract numeric part of sample names
sample_order <- sort(as.numeric(gsub(".*_(\\d+)$", "\\1", names(table(all_merged$sample)))))
# Rebuild full sample names in order
ordered_samples <- paste0("obcard_", sample_order)
# Set the factor order in the Seurat object
all_merged$sample <- factor(all_merged$sample, levels = ordered_samples)
# Plot ordered by sample name, no legend
# Percent mitochondrial content
VlnPlot(all_merged, features = "percent.mt",
  group.by = "sample", pt.size = 0.01,
  alpha = 0.1) + theme(legend.position = "none") +
  xlab(NULL) + ylab("Percent mitochondrial genes")
# Total RNA counts
VlnPlot(all_merged, features = "nCount_RNA",
  group.by = "sample", pt.size = 0.01,
  alpha = 0.1) + theme(legend.position = "none") +
  xlab(NULL) + ylab("Total RNA counts per cell")
# Number of detected genes
VlnPlot(all_merged, features = "nFeature_RNA",
  group.by = "sample", pt.size = 0.01,
  alpha = 0.1) + theme(legend.position = "none") +
  xlab(NULL) + ylab("Detected genes per cell")
qc <- FetchData(
  all_merged,
  vars = c("sample", "percent.mt", "nCount_RNA", "nFeature_RNA"))
summary_stats <- qc %>%
  group_by(sample) %>%
  summarise(
    cells            = n(),
    mean_percent_mt  = mean(percent.mt,  na.rm = TRUE),
    mean_nCount_RNA  = mean(nCount_RNA,  na.rm = TRUE),
    mean_nFeature_RNA= mean(nFeature_RNA,  na.rm = TRUE)) %>%
  arrange(sample)  
print(summary_stats, n = Inf)

# Filtering mitochondrial (remove >5% mitochondrial reads)
mito_filtered <- subset(all_merged, subset = percent.mt <= 5)
# 2) Count cells per sample (before vs after)
before <- table(all_merged$sample)
after  <- table(mito_filtered$sample)
# Neat summary table
cell_summary <- data.frame(
  sample       = union(names(before), names(after)),
  before       = as.integer(before[union(names(before), names(after))]),
  after        = as.integer(after[union(names(before), names(after))]))
cell_summary$removed     <- cell_summary$before - cell_summary$after
print(cell_summary[order(-cell_summary$after), ], row.names = FALSE)

## Filtering nFeature_RNA
qc_filtered <- subset(all_merged, subset = percent.mt <= 5 & 
                        nFeature_RNA > 200 & nFeature_RNA < 12500)
table(qc_filtered$sample)

# Save object + move onto doublet removal
saveRDS(qc_filtered, file = "/scratch/pigblast/nb443/GITHUB_FINALS/21_mito+feat_filtered.rds")

