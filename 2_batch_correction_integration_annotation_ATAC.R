# R 4.4.0 - run on node - 3/4 different 2 hour sessions 
# on bash - srun --account=preopenergy --partition=devel --x11 --cpus-per-task=4 --time=2:0:0 --mem=80g --pty /bin/bash
# ==== SESSION 1: GA prep + RPCA + anchors ====
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(future)
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)
Sys.setenv("OMP_NUM_THREADS"="1", "MKL_NUM_THREADS"="1")
set.seed(123)
# Paths
outdir <- "/scratch/pigblast/nb443/GITHUB_FINALS"
f_run1 <- file.path(outdir, "run1_multiome_merged_processed.rds")
f_run2 <- file.path(outdir, "run2_multiome_merged_processed.rds")
f_valc <- file.path(outdir, "4_control_valcards_filtered_merged_processed.RDS")
# Checkpoints produced here
f_ga_feats    <- file.path(outdir, "ATAC_GA_features_1500.rds")
f_ga_dims     <- file.path(outdir, "ATAC_GA_dims_1_20.rds")
f_ga_anchors  <- file.path(outdir, "ATAC_GA_RPCA_anchors.rds")
f_ga_objs_rpca<- file.path(outdir, "ATAC_GA_objs_afterRPCA.rds")  # optional, helps reuse PCs
# Load
run1_obcard <- readRDS(f_run1)
run2_obcard <- readRDS(f_run2)
valcard_sub <- readRDS(f_valc)
# Keep sample label
run1_obcard$sample <- run1_obcard$orig.ident
run2_obcard$sample <- run2_obcard$orig.ident
valcard_sub$sample <- valcard_sub$orig.ident
# Helpers
get_atac_assay <- function(obj) {
  poss <- intersect(c("ATAC","peaks","ChromatinAccessibility"), Assays(obj))
  if (length(poss) == 0) stop("No ATAC/peaks/ChromatinAccessibility assay found.")
  poss[1]}
to_std_atac <- function(x, atac_assay) {
  if (!atac_assay %in% Assays(x)) stop("Assay ", atac_assay, " not found.")
  mat <- GetAssayData(x, assay = atac_assay, layer = "counts")
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
    x[[atac_assay]] <- SetAssayData(x[[atac_assay]], layer = "counts", new.data = mat)}
  keep_assays <- c(atac_assay, intersect(c("GA","GeneActivity"), Assays(x)))
  DefaultAssay(x) <- atac_assay
  DietSeurat(x, assays = keep_assays, dimreducs = NULL, graphs = NULL)}
add_ga_assay <- function(obj, atac_assay) {
  if ("GA" %in% Assays(obj)) {
    DefaultAssay(obj) <- "GA"
  } else if ("GeneActivity" %in% Assays(obj)) {
    obj[["GA"]] <- obj[["GeneActivity"]]; obj[["GeneActivity"]] <- NULL
    DefaultAssay(obj) <- "GA"
  } else {
    DefaultAssay(obj) <- atac_assay
    ga <- GeneActivity(obj, assay = atac_assay)
    ga <- as(ga, "dgCMatrix")
    obj[["GA"]] <- CreateAssayObject(counts = ga)
    DefaultAssay(obj) <- "GA"}
  DefaultAssay(obj) <- "GA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj}
re_pca_ga <- function(x, feats, npcs = 20) {
  DefaultAssay(x) <- "GA"
  x <- ScaleData(x, features = feats, verbose = FALSE)
  x <- RunPCA(x, features = feats, npcs = npcs, verbose = FALSE)
  x}
# Slim to ATAC & GA
r1.atac <- get_atac_assay(run1_obcard)
r2.atac <- get_atac_assay(run2_obcard)
vc.atac <- get_atac_assay(valcard_sub)
run1_atac <- to_std_atac(run1_obcard, r1.atac)
run2_atac <- to_std_atac(run2_obcard, r2.atac)
valc_atac <- to_std_atac(valcard_sub,  vc.atac)
rm(run1_obcard, run2_obcard, valcard_sub); invisible(gc())
# Labels
run1_atac$dataset <- "ObCARD_MO_run1"
run2_atac$dataset <- "ObCARD_MO_run2"
valc_atac$dataset <- "ValCARD_ctrl"
# GA + initial PCA
run1_ga <- add_ga_assay(run1_atac, r1.atac)
run2_ga <- add_ga_assay(run2_atac, r2.atac)
valc_ga <- add_ga_assay(valc_atac,  vc.atac)
# Shared features + RPCA on shared features
objs_ga <- list(run1_ga, run2_ga, valc_ga)
for (i in seq_along(objs_ga)) DefaultAssay(objs_ga[[i]]) <- "GA"
ga_features <- SelectIntegrationFeatures(object.list = objs_ga, nfeatures = 1500)
run1_ga <- re_pca_ga(run1_ga, feats = ga_features, npcs = 20)
run2_ga <- re_pca_ga(run2_ga, feats = ga_features, npcs = 20)
valc_ga <- re_pca_ga(valc_ga, feats = ga_features, npcs = 20)
objs_ga <- list(run1_ga, run2_ga, valc_ga)
ga_dims  <- 1:20
anchors_ga <- FindIntegrationAnchors(
  object.list     = objs_ga,
  anchor.features = ga_features,
  reduction       = "rpca",
  dims            = ga_dims)
# SAVE & STOP
saveRDS(ga_features, f_ga_feats)
saveRDS(ga_dims,     f_ga_dims)
saveRDS(anchors_ga,  f_ga_anchors)
saveRDS(objs_ga,     f_ga_objs_rpca)   # optional but useful
cat("Saved:\n- ", f_ga_feats, "\n- ", f_ga_dims, "\n- ", f_ga_anchors, "\n- ", f_ga_objs_rpca, "\n")

# ==== SESSION 2: IntegrateData only (takes approx 2 hours) ====
library(Seurat)
library(Signac)
library(Matrix)
library(future)
plan(sequential)
options(future.globals.maxSize = 8 * 1024^3)
Sys.setenv("OMP_NUM_THREADS"="1", "MKL_NUM_THREADS"="1")
set.seed(123)
outdir <- "/scratch/pigblast/nb443/GITHUB_FINALS"
f_ga_dims    <- file.path(outdir, "ATAC_GA_dims_1_20.rds")
f_ga_anchors <- file.path(outdir, "ATAC_GA_RPCA_anchors.rds")
f_integr_ga  <- file.path(outdir, "ATAC_GA_RPCA_integrated.rds")
ga_dims    <- readRDS(f_ga_dims)
anchors_ga <- readRDS(f_ga_anchors)
# Lower k.weight to keep this step quick
integrated_ga <- IntegrateData(
  anchorset = anchors_ga,
  dims      = ga_dims,
  k.weight  = 25)
saveRDS(integrated_ga, f_integr_ga)
cat("Saved integrated object:\n- ", f_integr_ga, "\n")

# ==== SESSION 3: UMAP / neighbors / clusters ====
library(Seurat)
library(Matrix)
library(ggplot2)
set.seed(123)
outdir      <- "/scratch/pigblast/nb443/GITHUB_FINALS"
f_integr_ga <- file.path(outdir, "ATAC_GA_RPCA_integrated.rds")
f_umap_ga   <- file.path(outdir, "ATAC_GA_RPCA_integrated_postUMAP.rds")
integrated_ga <- readRDS(f_integr_ga)
DefaultAssay(integrated_ga) <- "integrated"
# Ensure features available
mat <- tryCatch(GetAssayData(integrated_ga, assay = "integrated", layer = "data"),
                error = function(e) NULL)
if (is.null(mat) || nrow(mat) == 0) {
  mat <- GetAssayData(integrated_ga, assay = "integrated", layer = "counts")}
nz <- Matrix::rowSums(mat != 0)
vf <- head(names(sort(nz, decreasing = TRUE)), 2000)
VariableFeatures(integrated_ga) <- vf
# Use the same dims as Session 1
ga_dims <- readRDS(file.path(outdir, "ATAC_GA_dims_1_20.rds"))
integrated_ga <- ScaleData(integrated_ga, features = vf, verbose = FALSE)
integrated_ga <- RunPCA(integrated_ga, features = vf, npcs = 30, verbose = FALSE)
integrated_ga <- RunUMAP(integrated_ga, dims = ga_dims, verbose = FALSE)
integrated_ga <- FindNeighbors(integrated_ga, dims = ga_dims, verbose = FALSE)
integrated_ga <- FindClusters(integrated_ga, resolution = 0.5, verbose = FALSE)
# Keep dataset/sample if present in meta
if (!"dataset" %in% colnames(integrated_ga@meta.data)) {
  integrated_ga$dataset <- if ("orig.ident" %in% colnames(integrated_ga@meta.data)) integrated_ga$orig.ident else "batch"}
if (!"sample" %in% colnames(integrated_ga@meta.data)) {
  integrated_ga$sample <- integrated_ga$orig.ident}
saveRDS(integrated_ga, f_umap_ga)
cat("Saved object with UMAP/clusters:\n- ", f_umap_ga, "\n")
DimPlot(integrated_ga, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + NoLegend()
DimPlot(integrated_ga, reduction = "umap", group.by = "dataset")

# ==== SESSION 4: Cross-modal label transfer (RNA ref -> GA) ====
library(Seurat)
library(Signac)
library(Matrix)
library(ggplot2)
set.seed(123)
outdir     <- "/scratch/pigblast/nb443/GITHUB_FINALS"
f_umap_ga  <- file.path(outdir, "ATAC_GA_RPCA_integrated_postUMAP.rds")
f_final_ga <- file.path(outdir, "ATAC_GA_RPCA_integrated_annotated.rds")
ref_path   <- "/scratch/pigblast/nb443/heart_reference_with_sct.rds"
# Load GA integrated object (with UMAP)
integrated_ga <- readRDS(f_umap_ga)
# Load RNA reference (SCT)
reference <- readRDS(ref_path)
reference <- UpdateSeuratObject(reference)
DefaultAssay(reference) <- "SCT"
reference <- RunPCA(reference, npcs = 50, verbose = FALSE)
reference <- RunUMAP(reference, reduction = "pca", dims = 1:50, return.model = TRUE, verbose = FALSE)
# Prepare GA query (you already did this)
ga_query <- integrated_ga
DefaultAssay(ga_query) <- "GA"
ga_query <- NormalizeData(ga_query, verbose = FALSE)
# Label column in reference
label_col_in_ref <- "cell_type_leiden0.5"  # change if needed
# Corrected: use SCT normalization, reference.assay = "SCT", query.assay = "GA"
anchors_xm <- FindTransferAnchors(
  reference           = reference,
  query               = ga_query,
  normalization.method= "SCT",
  reference.assay     = "SCT",
  query.assay         = "GA",
  reference.reduction = "pca",   
  dims                = 1:50,
  k.filter            = NA,
  verbose             = FALSE)
ga_query <- MapQuery(
  anchorset           = anchors_xm,
  query               = ga_query,
  reference           = reference,
  refdata             = list(seurat_label = label_col_in_ref),
  reference.reduction = "pca",
  reduction.model     = "umap")
# Copy labels back
integrated_ga$seurat_label <- ga_query$predicted.seurat_label
# Plots 
p1 <- DimPlot(integrated_ga, reduction = "umap", group.by = "seurat_label", label = TRUE, repel = TRUE) +
   ggtitle("ATAC (GA) — Transferred cell types")
p2 <- DimPlot(integrated_ga, reduction = "umap", group.by = "dataset") +
   ggtitle("ATAC (GA) — Dataset")
print(p1); print(p2)
saveRDS(integrated_ga, f_final_ga)
cat("Saved final annotated object:\n- ", f_final_ga, "\n")

