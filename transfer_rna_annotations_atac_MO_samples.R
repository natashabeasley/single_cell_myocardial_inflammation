# Run parts A, B and C on new R session on node (prevents crashing + reaching the 2 hour time limit)
## ===================== PART A =====================
## RNA MULTIOME-ONLY RPCA INTEGRATION + ANNOTATION
## =================================================
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(future)
plan(sequential)
options(future.globals.maxSize = 256 * 1024^2)
Sys.setenv("OMP_NUM_THREADS" = "1", "MKL_NUM_THREADS" = "1")
set.seed(123)
## ---------- Paths ----------
outdir <- "/scratch/pigblast/nb443/GITHUB_FINALS"
f_run1 <- file.path(outdir, "run1_multiome_merged_processed.rds")
f_run2 <- file.path(outdir, "run2_multiome_merged_processed.rds")
f_valc <- file.path(outdir, "4_control_valcards_filtered_merged_processed.RDS")
f_rna_integrated <- file.path(outdir, "MO_RNA_RPCA_integrated_keepALL_lowRAM.rds")
f_rna_annotated  <- file.path(outdir, "MO_RNA_RPCA_integrated_annotated.rds")
ref_path <- "/scratch/pigblast/nb443/heart_reference_with_sct.rds"
## ---------- Load multiome objects (RNA+ATAC) ----------
run1_obcard <- readRDS(f_run1)
run2_obcard <- readRDS(f_run2)
valcard_sub <- readRDS(f_valc)
# Keep per-sample IDs (important!)
run1_obcard$sample <- run1_obcard$orig.ident
run2_obcard$sample <- run2_obcard$orig.ident
valcard_sub$sample <- valcard_sub$orig.ident
# Dataset labels (batch-level)
run1_obcard$dataset <- "ObCARD_MO_run1"
run2_obcard$dataset <- "ObCARD_MO_run2"
valcard_sub$dataset <- "ValCARD_ctrl"
## ---------- RNA assay ----------
to_std_rna <- function(x) {
  get_rna <- function(obj) {
    if ("RNA" %in% Assays(obj)) {
      lyr <- try(Layers(obj[["RNA"]]), silent = TRUE)
      if (!inherits(lyr, "try-error") && length(lyr) > 0) {
        if ("counts" %in% lyr) return(GetAssayData(obj, assay = "RNA", layer = "counts"))
        if ("data"   %in% lyr) return(GetAssayData(obj, assay = "RNA", layer = "data"))} else {
        return(GetAssayData(obj, assay = "RNA"))}}
    if ("GeneExpression" %in% Assays(obj)) {
      lyr <- try(Layers(obj[["GeneExpression"]]), silent = TRUE)
      if (!inherits(lyr, "try-error") && length(lyr) > 0) {
        if ("counts" %in% lyr) return(GetAssayData(obj, assay = "GeneExpression", layer = "counts"))
        if ("data"   %in% lyr) return(GetAssayData(obj, assay = "GeneExpression", layer = "data"))} else {
        return(GetAssayData(obj, assay = "GeneExpression"))}}
    if ("SCT" %in% Assays(obj)) {
      return(GetAssayData(obj, assay = "SCT", layer = "data"))}
    stop("No usable RNA/GeneExpression/SCT assay found.")}
  mat <- get_rna(x)
  stopifnot(!is.null(rownames(mat)), !is.null(colnames(mat)))
  if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
  x[["RNA"]] <- CreateAssayObject(counts = mat)
  DefaultAssay(x) <- "RNA"
  # Drop all *other* assays/dimreds/graphs to save RAM
  for (a in setdiff(Assays(x), "RNA")) {
    x[[a]] <- NULL}
  x <- DietSeurat(x, assays = "RNA", dimreducs = NULL, graphs = NULL)
  x}
run1_rna_rpca <- to_std_rna(run1_obcard)
run2_rna_rpca <- to_std_rna(run2_obcard)
valc_rna_rpca <- to_std_rna(valcard_sub)
rm(run1_obcard, run2_obcard, valcard_sub); gc()
## ---------- Preprocess each (no per-sample split) ----------
prep <- function(x, nfeat = 2000, npcs = 30) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = nfeat, verbose = FALSE)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = npcs, verbose = FALSE)
  x}
run1_rna_rpca <- prep(run1_rna_rpca, nfeat = 2000, npcs = 30)
run2_rna_rpca <- prep(run2_rna_rpca, nfeat = 2000, npcs = 30)
valc_rna_rpca <- prep(valc_rna_rpca, nfeat = 2000, npcs = 30)
## ---------- Select features and re-PCA ----------
objs     <- list(run1_rna_rpca, run2_rna_rpca, valc_rna_rpca)
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 2000)
repc <- function(x, feats = features, npcs = 30) {
  DefaultAssay(x) <- "RNA"
  x <- ScaleData(x, features = feats, verbose = FALSE)
  x <- RunPCA(x, features = feats, npcs = npcs, verbose = FALSE)
  x}
run1_rna_rpca <- repc(run1_rna_rpca)
run2_rna_rpca <- repc(run2_rna_rpca)
valc_rna_rpca <- repc(valc_rna_rpca)
objs <- list(run1_rna_rpca, run2_rna_rpca, valc_rna_rpca)
min_cells <- min(vapply(objs, ncol, integer(1)))
dims_use  <- 1:max(2, min(30, min_cells - 1))
message("Using dims 1:", max(dims_use), " (min cells = ", min_cells, ")")
## ---------- Anchors (RPCA) ----------
anchors <- FindIntegrationAnchors(
  object.list     = objs,
  reduction       = "rpca",
  anchor.features = features,
  dims            = dims_use)
## ---------- Integrate ----------
integrated_rna <- tryCatch(
  IntegrateData(anchorset = anchors, dims = dims_use, k.weight = 5),
  error = function(e) {
    message("Direct integration failed: ", conditionMessage(e), " — retrying with reference.")
    sizes   <- vapply(objs, ncol, integer(1))
    ref_idx <- which.max(sizes)
    anchors2 <- FindIntegrationAnchors(
      object.list     = objs,
      reference       = ref_idx,
      reduction       = "rpca",
      anchor.features = features,
      dims            = dims_use)
    IntegrateData(anchorset = anchors2, dims = dims_use, k.weight = 5)})
DefaultAssay(integrated_rna) <- "integrated"
## ---------- Ensure 'integrated' has a proper 'data' layer ----------
get_mat <- function(obj) {
  m <- tryCatch(GetAssayData(obj, assay = "integrated", layer = "data"),   error = function(e) NULL)
  if (is.null(m)) m <- tryCatch(GetAssayData(obj, assay = "integrated", layer = "counts"), error = function(e) NULL)
  if (is.null(m)) m <- tryCatch(GetAssayData(obj, slot = "data"),          error = function(e) NULL)
  m}
mat <- get_mat(integrated_rna)
mat <- as.matrix(mat)
rownames(mat) <- rownames(integrated_rna)
colnames(mat) <- colnames(integrated_rna)
integrated_rna[["integrated"]] <- SetAssayData(integrated_rna[["integrated"]], layer = "data", new.data = mat)
## ---------- Downstream: UMAP / clusters ----------
nz       <- Matrix::rowSums(mat != 0)
features <- head(names(sort(nz, decreasing = TRUE)), 2000)
features <- intersect(features, rownames(integrated_rna))
VariableFeatures(integrated_rna) <- features
integrated_rna <- ScaleData(integrated_rna, features = features, verbose = FALSE)
integrated_rna <- RunPCA(integrated_rna, features = features, npcs = 30, verbose = FALSE)
integrated_rna <- RunUMAP(integrated_rna, dims = dims_use, verbose = FALSE)
integrated_rna <- FindNeighbors(integrated_rna, dims = dims_use, verbose = FALSE)
integrated_rna <- FindClusters(integrated_rna, resolution = 0.5, verbose = FALSE)
# Keep sample + dataset
if (!"sample" %in% colnames(integrated_rna@meta.data)) {
  integrated_rna$sample <- integrated_rna$orig.ident}
if (!"dataset" %in% colnames(integrated_rna@meta.data)) {
  integrated_rna$dataset <- "multiome"}
cat("Cells per dataset:\n")
print(table(integrated_rna$dataset))
cat("Cells per sample:\n")
print(table(integrated_rna$sample))
saveRDS(integrated_rna, f_rna_integrated)
cat("\nSaved RNA multiome-only integrated object to:\n  ", f_rna_integrated, "\n")
## ---------- Annotation vs heart reference ----------
obj       <- integrated_rna
reference <- readRDS(ref_path)
reference <- UpdateSeuratObject(reference)
DefaultAssay(reference) <- "SCT"
reference <- RunUMAP(
  reference,
  reduction    = "pca",
  dims         = 1:50,
  return.model = TRUE,
  verbose      = FALSE)
# Light SCT of the query (adds SCT assay)
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE)
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
anchors_ref <- FindTransferAnchors(
  reference            = reference,
  query                = obj,
  normalization.method = "SCT",
  reference.reduction  = "pca",
  dims                 = 1:50,
  k.filter             = NA,
  verbose              = FALSE)
obj <- MapQuery(
  anchorset           = anchors_ref,
  query               = obj,
  reference           = reference,
  refdata             = list(seurat_label = "cell_type_leiden0.5"),
  reference.reduction = "pca",
  reduction.model     = "umap")
# Use seurat_label as RNA cell-type field
DefaultAssay(obj) <- "integrated"
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)
if (!"dataset" %in% colnames(obj@meta.data)) {
  obj$dataset <- if ("orig.ident" %in% colnames(obj@meta.data)) obj$orig.ident else "batch"}
# Just to be explicit:
obj$seurat_label <- obj$seurat_label
saveRDS(obj, f_rna_annotated)
cat("\nSaved RNA multiome-only integrated + annotated object to:\n  ", f_rna_annotated, "\n")

## ===================== PART B1 =====================
## Build ATAC GA + RNA multiome-only & preprocess for CCA
## =====================================================
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(ggplot2)
set.seed(123)
outdir <- "/scratch/pigblast/nb443/GITHUB_FINALS"
# Inputs
f_run1      <- file.path(outdir, "run1_multiome_merged_processed.rds")
f_run2      <- file.path(outdir, "run2_multiome_merged_processed.rds")
f_valc      <- file.path(outdir, "4_control_valcards_filtered_merged_processed.RDS")
f_rna_annot <- file.path(outdir, "MO_RNA_RPCA_integrated_annotated.rds")
# Outputs from this session (for Session B2)
f_rna_mo_prepped  <- file.path(outdir, "MO_RNA_multiome_only_PREPPED_for_CCA.rds")
f_atac_co_prepped <- file.path(outdir, "MO_ATAC_GA_PREPPED_for_CCA.rds")
## ---------- Load annotated RNA and raw multiome objects ----------
rna   <- readRDS(f_rna_annot)
run1  <- readRDS(f_run1)
run2  <- readRDS(f_run2)
valc  <- readRDS(f_valc)
# Make sure RNA has per-sample IDs and labels
if (!"sample" %in% colnames(rna@meta.data)) {
  rna$sample <- rna$orig.ident}
stopifnot("seurat_label" %in% colnames(rna@meta.data))
## ---------- Helper functions for ATAC/GA ----------
get_atac_assay <- function(obj) {
  poss <- intersect(c("ATAC", "peaks", "ChromatinAccessibility"), Assays(obj))
  if (length(poss) == 0) stop("No ATAC/peaks/ChromatinAccessibility assay found.")
  poss[1]}
to_std_atac <- function(x, atac_assay) {
  if (!atac_assay %in% Assays(x)) stop("Assay ", atac_assay, " not found.")
  mat <- GetAssayData(x, assay = atac_assay, layer = "counts")
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
    x[[atac_assay]] <- SetAssayData(x[[atac_assay]], layer = "counts", new.data = mat)}
  DefaultAssay(x) <- atac_assay
  keep_assays <- c(atac_assay, intersect(c("GA", "GeneActivity"), Assays(x)))
  x <- DietSeurat(x, assays = keep_assays, dimreducs = NULL, graphs = NULL)
  x}
add_ga_assay <- function(obj, atac_assay) {
  if ("GA" %in% Assays(obj)) {
    DefaultAssay(obj) <- "GA"
  } else if ("GeneActivity" %in% Assays(obj)) {
    obj[["GA"]] <- obj[["GeneActivity"]]
    obj[["GeneActivity"]] <- NULL
    DefaultAssay(obj) <- "GA"} else {
    DefaultAssay(obj) <- atac_assay
    ga <- GeneActivity(obj, assay = atac_assay)
    ga <- as(ga, "dgCMatrix")
    obj[["GA"]] <- CreateAssayObject(counts = ga)
    DefaultAssay(obj) <- "GA"}
  obj}
## ---------- Slim to ATAC + GA and merge into one ATAC object ----------
r1.atac <- get_atac_assay(run1)
r2.atac <- get_atac_assay(run2)
vc.atac <- get_atac_assay(valc)
run1_atac <- to_std_atac(run1, r1.atac)
run2_atac <- to_std_atac(run2, r2.atac)
valc_atac <- to_std_atac(valc, vc.atac)
rm(run1, run2, valc); gc()
# Add GA
run1_atac <- add_ga_assay(run1_atac, r1.atac)
run2_atac <- add_ga_assay(run2_atac, r2.atac)
valc_atac <- add_ga_assay(valc_atac, vc.atac)
# sample + dataset meta
run1_atac$sample  <- run1_atac$orig.ident
run2_atac$sample  <- run2_atac$orig.ident
valc_atac$sample  <- valc_atac$orig.ident
run1_atac$dataset <- "ObCARD_MO_run1"
run2_atac$dataset <- "ObCARD_MO_run2"
valc_atac$dataset <- "ValCARD_ctrl"
# Merge GA objects 
objs_ga <- list(run1_atac, run2_atac, valc_atac)
atac <- merge(
  x = objs_ga[[1]],
  y = objs_ga[2:3],
  add.cell.ids = c("run1", "run2", "val"))
DefaultAssay(atac) <- "GA"
if (!"sample" %in% colnames(atac@meta.data))  atac$sample  <- atac$orig.ident
if (!"dataset" %in% colnames(atac@meta.data)) atac$dataset <- "multiome_ATAC"
cat("Merged ATAC GA cells: ", ncol(atac), "\n")
## ---------- Restrict RNA to multiome samples ----------
multiome_samples <- sort(unique(atac$orig.ident))
cat("Multiome samples in ATAC:\n")
print(multiome_samples)
rna_mo <- subset(rna, subset = sample %in% multiome_samples)
rna_mo$sample <- rna_mo$sample
cat("RNA multiome-only cells: ", ncol(rna_mo), "\n")
## ---------- Make GA look like 'RNA' and drop extra assays ----------
atac_co <- atac
ga_assay <- atac_co[["GA"]]
atac_co[["RNA"]] <- ga_assay
DefaultAssay(atac_co) <- "RNA"
for (a in setdiff(Assays(atac_co), "RNA")) {
  atac_co[[a]] <- NULL}
## ---------- Ensure common genes ----------
common_genes <- intersect(rownames(rna_mo), rownames(atac_co))
cat("Common genes between RNA and ATAC GA: ", length(common_genes), "\n")
rna_mo  <- subset(rna_mo,  features = common_genes)
atac_co <- subset(atac_co, features = common_genes)
## ---------- Rename cells to track modality ----------
rna_mo  <- RenameCells(rna_mo,  add.cell.id = "RNA")
atac_co <- RenameCells(atac_co, add.cell.id = "ATAC")
rna_mo$modality  <- "RNA"
atac_co$modality <- "ATAC"
## ---------- Preprocess (Normalize, HVGs, Scale, PCA) ----------
prep_fun <- function(x, nfeat = 2000, npcs = 30) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, nfeatures = nfeat, verbose = FALSE)
  x <- ScaleData(x, verbose = FALSE)
  x <- RunPCA(x, npcs = npcs, verbose = FALSE)
  x}
rna_mo  <- prep_fun(rna_mo)
atac_co <- prep_fun(atac_co)
## ---------- Save prepped objects for Session B2 ----------
saveRDS(rna_mo,  f_rna_mo_prepped)
saveRDS(atac_co, f_atac_co_prepped)
cat("\nSaved prepped objects for CCA coembedding to:\n")
cat("  RNA : ", f_rna_mo_prepped,  "\n")
cat("  ATAC: ", f_atac_co_prepped, "\n")

## ===================== PART C =====================
## CCA coembedding + FORCED RNA→ATAC labels
## ================================================
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(ggplot2)
library(tidyr)
set.seed(123)
outdir <- "/scratch/pigblast/nb443/GITHUB_FINALS"
# Inputs from PART B
f_rna_mo_prepped  <- file.path(outdir, "MO_RNA_multiome_only_PREPPED_for_CCA.rds")
f_atac_co_prepped <- file.path(outdir, "MO_ATAC_GA_PREPPED_for_CCA.rds")
# Final output
f_coembed <- file.path(outdir, "MO_RNA_ATAC_GA_coembed_CCA_forcedLabels.rds")
## ---------- Load prepped RNA & ATAC ----------
rna_mo  <- readRDS(f_rna_mo_prepped)
atac_co <- readRDS(f_atac_co_prepped)
cat("Prepped RNA cells:  ", ncol(rna_mo),  "\n")
cat("Prepped ATAC cells: ", ncol(atac_co), "\n")
# Modality, sample should exist
if (!"modality" %in% colnames(rna_mo@meta.data))  rna_mo$modality  <- "RNA"
if (!"modality" %in% colnames(atac_co@meta.data)) atac_co$modality <- "ATAC"
if (!"sample" %in% colnames(rna_mo@meta.data))  rna_mo$sample  <- rna_mo$orig.ident
if (!"sample" %in% colnames(atac_co@meta.data)) atac_co$sample <- atac_co$orig.ident
## ---------- CCA-based coembedding (as in v3.1 ATAC vignette) ----------
objs <- list(rna_mo, atac_co)
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 2000)
repc_fun <- function(x, feats, npcs = 30) {
  DefaultAssay(x) <- "RNA"
  x <- ScaleData(x, features = feats, verbose = FALSE)
  x <- RunPCA(x, features = feats, npcs = npcs, verbose = FALSE)
  x}
rna_mo  <- repc_fun(rna_mo,  features)
atac_co <- repc_fun(atac_co, features)
objs     <- list(rna_mo, atac_co)
dims_use <- 1:30
anchors_cca <- FindIntegrationAnchors(
  object.list     = objs,
  anchor.features = features,
  reduction       = "cca",
  dims            = dims_use)
coembed <- IntegrateData(
  anchorset = anchors_cca,
  dims      = dims_use)
DefaultAssay(coembed) <- "integrated"
coembed <- ScaleData(coembed, verbose = FALSE)
coembed <- RunPCA(coembed, npcs = 30, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = dims_use, verbose = FALSE)
coembed <- FindNeighbors(coembed, dims = dims_use, verbose = FALSE)
coembed <- FindClusters(coembed, resolution = 0.5, verbose = FALSE)
## ---------- Recover modality + sample into coembed ----------
coembed$modality <- ifelse(grepl("^RNA_",  colnames(coembed)), "RNA", "ATAC")
# Build sample vector using lookups from rna_mo / atac_co
rna_sample_lookup  <- rna_mo$sample
names(rna_sample_lookup) <- colnames(rna_mo)
atac_sample_lookup <- atac_co$sample
names(atac_sample_lookup) <- colnames(atac_co)
sample_vec <- rep(NA_character_, ncol(coembed))
names(sample_vec) <- colnames(coembed)
rna_idx  <- grepl("^RNA_",  colnames(coembed))
atac_idx <- grepl("^ATAC_", colnames(coembed))
sample_vec[rna_idx]  <- rna_sample_lookup[colnames(coembed)[rna_idx]]
sample_vec[atac_idx] <- atac_sample_lookup[colnames(coembed)[atac_idx]]
coembed$sample <- sample_vec
## ---------- Seed RNA labels into coembed ----------
# Ensure coembed has a seurat_label column
if (!"seurat_label" %in% colnames(coembed@meta.data)) {
  coembed$seurat_label <- NA_character_}
# Fill RNA cell labels from rna_mo
rna_label_lookup <- rna_mo$seurat_label
names(rna_label_lookup) <- colnames(rna_mo)
coembed$seurat_label[rna_idx] <- rna_label_lookup[colnames(coembed)[rna_idx]]
## ---------- Build RNA–ATAC 1:1 pairing from PREPPED objects ----------
rna_cells_pre  <- colnames(rna_mo)
atac_cells_pre <- colnames(atac_co)
# Extract the core 16bp barcode from RNA and ATAC cell names
rna_base <- sub("^RNA_([^ -]+)-.*", "\\1", rna_cells_pre)
atac_base <- sub("^ATAC_[^_]+_([^ -]+)-.*", "\\1", atac_cells_pre)
cat("Example RNA barcodes:\n")
print(head(rna_base))
cat("Example ATAC barcodes:\n")
print(head(atac_base))
cat("Unique RNA barcodes:  ", length(unique(rna_base)),  "\n")
cat("Unique ATAC barcodes: ", length(unique(atac_base)), "\n")
# Build keys using sample + base barcode
rna_key  <- paste(rna_mo$sample,  rna_base,  sep = "__")
atac_key <- paste(atac_co$sample, atac_base, sep = "__")
cat("Unique RNA keys (prepped):  ", length(unique(rna_key)),  "\n")
cat("Unique ATAC keys (prepped): ", length(unique(atac_key)), "\n")
common_keys <- intersect(rna_key, atac_key)
cat("Matched RNA–ATAC keys (prepped): ", length(common_keys), "\n")
# Map keys back to cell names
rna_key_to_cell_pre  <- setNames(rna_cells_pre,  rna_key)
atac_key_to_cell_pre <- setNames(atac_cells_pre, atac_key)
rna_cells_common_pre  <- rna_key_to_cell_pre[common_keys]
atac_cells_common_pre <- atac_key_to_cell_pre[common_keys]
stopifnot(!any(is.na(rna_cells_common_pre)))
stopifnot(!any(is.na(atac_cells_common_pre)))
stopifnot(length(rna_cells_common_pre) == length(atac_cells_common_pre))
cat("Final matched cell pairs (prepped): ", length(rna_cells_common_pre), "\n")
# Check samples match
rna_idx  <- match(rna_cells_common_pre,  rna_cells_pre)
atac_idx <- match(atac_cells_common_pre, atac_cells_pre)
same_sample <- all(rna_mo$sample[rna_idx] == atac_co$sample[atac_idx])
cat("Do RNA & ATAC samples match for all pairs? ", same_sample, "\n")
## ---------- Translate those pairs into coembed cell names ----------
# colnames in coembed are the same as in rna_mo / atac_co
rna_cells_common_coembed  <- rna_cells_common_pre
atac_cells_common_coembed <- atac_cells_common_pre
# Ensure they all exist in coembed
keep <- rna_cells_common_coembed  %in% colnames(coembed) &
  atac_cells_common_coembed %in% colnames(coembed)
rna_cells_common_coembed  <- rna_cells_common_coembed[keep]
atac_cells_common_coembed <- atac_cells_common_coembed[keep]
cat("Matched pairs present in coembed: ", length(rna_cells_common_coembed), "\n")
## ---------- FORCE ATAC labels = RNA labels ----------
rna_labels_for_pairs <- coembed$seurat_label[rna_cells_common_coembed]
# Overwrite ATAC labels with the RNA labels
coembed$seurat_label[atac_cells_common_coembed] <- rna_labels_for_pairs
cat("ATAC seurat_label has been overwritten to match RNA for all matched pairs.\n")
## ---------- Optional: per-sample RNA vs ATAC proportions (matched cells only) ----------
meta <- data.frame(
  cell     = colnames(coembed),
  sample   = coembed$sample,
  modality = coembed$modality,
  label    = coembed$seurat_label,
  stringsAsFactors = FALSE)
meta_matched <- meta[meta$cell %in% c(rna_cells_common_coembed, atac_cells_common_coembed), ]
counts <- meta_matched %>%
  group_by(sample, modality, label) %>%
  summarise(n = n(), .groups = "drop")
props <- counts %>%
  group_by(sample, modality) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
compare_props <- props %>%
  select(sample, modality, label, prop) %>%
  pivot_wider(
    names_from  = modality,
    values_from = prop,
    values_fill = 0)
cat("\nFirst few rows of per-sample label proportions (RNA vs ATAC, matched cells only):\n")
print(head(compare_props, 20))
## ---------- Plots ----------
DimPlot(coembed,
  reduction = "umap",
  group.by  = "modality") + ggtitle("RNA + ATAC GA coembedded UMAP (modality, CCA)")
DimPlot(coembed,
  reduction = "umap",
  group.by  = "seurat_label",
  label     = TRUE,
  repel     = TRUE) + ggtitle("RNA + ATAC GA coembedded UMAP (forced RNA labels)")
## ---------- Save final coembedded object ----------
saveRDS(coembed, f_coembed)
cat("\nSaved RNA+ATAC GA coembedded object with forced labels to:\n  ", f_coembed, "\n")
