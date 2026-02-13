#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required to generate simulations.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required to generate simulations.")
  }
})

.usage <- function() {
  cat(
"simulate_inputs.R

Usage:
  Rscript simulate_inputs.R [options]

Options:
  --outdir <dir>    Output directory (default: ./simulated_inputs)
  --seed <int>      RNG seed (default: 123)
  --help            Print help and exit

Outputs:
  - s5f_gio.tab.uap.integrated.wt.lm_drop.rds
  - scMuscle_mm10_slim_v1-1.RData   (contains object named: seu)
  - SIMULATION_NOTES.md
"
  )
}

.parse_args <- function(args) {
  opts <- list(
    outdir = "simulated_inputs",
    seed = 123L
  )

  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("-h", "--help")) {
      opts$help <- TRUE
      i <- i + 1L
      next
    }
    if (grepl("^--", a)) {
      key <- sub("^--", "", a)
      if (key %in% c("outdir", "seed")) {
        if (i == length(args) || grepl("^-", args[[i + 1L]])) {
          stop("Missing value for --", key)
        }
        val <- args[[i + 1L]]
        i <- i + 1L
        if (key == "seed") val <- as.integer(val)
        opts[[key]] <- val
      } else {
        stop("Unknown option: ", a)
      }
      i <- i + 1L
      next
    }
    stop("Unrecognized argument: ", a)
  }

  opts
}

.make_counts <- function(n_genes, n_cells, lambda = 3, dropout = 0.85, gene_prefix = "Gene", cell_prefix = "Cell") {
  x <- matrix(stats::rpois(n_genes * n_cells, lambda = lambda), nrow = n_genes, ncol = n_cells)
  drop_mask <- matrix(stats::runif(n_genes * n_cells) < dropout, nrow = n_genes, ncol = n_cells)
  x[drop_mask] <- 0L

  rownames(x) <- sprintf("%s%04d", gene_prefix, seq_len(n_genes))
  colnames(x) <- sprintf("%s%04d", cell_prefix, seq_len(n_cells))
  x
}

.safe_scale <- function(mat) {
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

.add_feature_metadata <- function(seu, assay, variable_features, style = c("vst", "mvp")) {
  style <- match.arg(style)
  feats <- rownames(SeuratObject::GetAssayData(seu, assay = assay, layer = "counts"))

  if (style == "vst") {
    fmd <- data.frame(
      vst.mean = stats::runif(length(feats), min = 0, max = 8),
      vst.variance = stats::rgamma(length(feats), shape = 2, rate = 1),
      vst.variance.standardized = stats::rnorm(length(feats), mean = 1, sd = 0.5),
      highly_variable = feats %in% variable_features,
      gene_biotype = sample(c("protein_coding", "lncRNA", "pseudogene"), length(feats), replace = TRUE),
      row.names = feats,
      stringsAsFactors = FALSE
    )
  } else {
    fmd <- data.frame(
      mean = stats::runif(length(feats), min = 0, max = 6),
      variance = stats::rgamma(length(feats), shape = 1.8, rate = 1),
      dispersion = stats::runif(length(feats), min = 0, max = 5),
      dispersion.scaled = stats::rnorm(length(feats), mean = 0, sd = 1),
      highly_variable_rank = ifelse(feats %in% variable_features, match(feats, variable_features), NA_integer_),
      row.names = feats,
      stringsAsFactors = FALSE
    )
  }

  a <- seu[[assay]]
  a@meta.data <- fmd
  seu[[assay]] <- a
  seu
}

.add_reductions <- function(seu, assay, n_pcs = 30L) {
  ncells <- ncol(seu)

  pca <- matrix(stats::rnorm(ncells * n_pcs), nrow = ncells, ncol = n_pcs)
  rownames(pca) <- colnames(seu)
  colnames(pca) <- paste0("PC_", seq_len(n_pcs))
  seu[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = pca, key = "PC_", assay = assay)

  umap <- matrix(stats::rnorm(ncells * 2L), nrow = ncells, ncol = 2L)
  rownames(umap) <- colnames(seu)
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  seu[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = assay)

  seu
}

.build_uapinyoying_like <- function() {
  counts <- .make_counts(
    n_genes = 2200L,
    n_cells = 650L,
    lambda = 3,
    dropout = 0.82,
    gene_prefix = "UAPG",
    cell_prefix = "UAPC"
  )

  seu <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    assay = "RNA",
    project = "Uapinyoying2023_sim"
  )

  md <- seu@meta.data
  ncells <- ncol(seu)
  md$sample_id <- sample(sprintf("gio_%02d", 1:8), ncells, replace = TRUE)
  md$condition <- sample(c("WT", "LM_drop"), ncells, replace = TRUE, prob = c(0.55, 0.45))
  md$batch <- sample(c("batchA", "batchB", "batchC"), ncells, replace = TRUE)
  md$celltype <- factor(sample(c("FAP", "Mesenchymal", "Pericyte", "Immune"), ncells, replace = TRUE, prob = c(0.45, 0.3, 0.15, 0.1)))
  md$percent.mt <- round(stats::runif(ncells, min = 0, max = 18), 2)
  md$pct_counts_mt <- md$percent.mt
  md$qc_notes <- I(replicate(ncells, sample(c("ok", "high_mt", "low_genes"), 2, replace = TRUE), simplify = FALSE))
  seu@meta.data <- md
  SeuratObject::Idents(seu) <- md$celltype

  counts_mat <- as.matrix(SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts"))
  data_mat <- log1p(counts_mat)
  scale_mat <- .safe_scale(data_mat)

  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data", new.data = Matrix::Matrix(data_mat, sparse = TRUE))
  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "scale.data", new.data = scale_mat)

  # Add split-like layers to exercise join-layer code paths.
  counts_sparse <- SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts")
  data_sparse <- SeuratObject::GetAssayData(seu, assay = "RNA", layer = "data")
  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "counts.batchA", new.data = counts_sparse)
  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "counts.batchB", new.data = counts_sparse)
  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data.batchA", new.data = data_sparse)
  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data.batchB", new.data = data_sparse)

  vfs <- sample(rownames(seu), size = 1200L)
  SeuratObject::VariableFeatures(seu, assay = "RNA") <- vfs
  seu <- .add_feature_metadata(seu, assay = "RNA", variable_features = vfs, style = "vst")
  SeuratObject::VariableFeatures(seu, assay = "RNA") <- vfs

  # Add SCT assay to emulate SCTransform-present signaling.
  sct_counts <- counts_mat + matrix(stats::rpois(length(counts_mat), lambda = 1), nrow = nrow(counts_mat), ncol = ncol(counts_mat))
  rownames(sct_counts) <- rownames(counts_mat)
  colnames(sct_counts) <- colnames(counts_mat)
  seu[["SCT"]] <- SeuratObject::CreateAssay5Object(counts = Matrix::Matrix(sct_counts, sparse = TRUE))
  seu <- SeuratObject::SetAssayData(seu, assay = "SCT", layer = "data", new.data = Matrix::Matrix(log1p(sct_counts), sparse = TRUE))

  seu <- .add_reductions(seu, assay = "RNA", n_pcs = 30L)

  seu@commands <- list(
    NormalizeData.RNA = list(params = list(normalization.method = "LogNormalize")),
    FindVariableFeatures.RNA = list(params = list(selection.method = "vst")),
    ScaleData.RNA = list(params = list(features = "variable")),
    RunPCA.RNA = list(params = list(npcs = 30)),
    RunUMAP.RNA = list(params = list(dims = 1:30)),
    SCTransform.RNA = list(params = list(variable.features.n = 3000))
  )

  seu
}

.build_mckellard_like <- function() {
  counts <- .make_counts(
    n_genes = 1800L,
    n_cells = 520L,
    lambda = 2,
    dropout = 0.86,
    gene_prefix = "MMG",
    cell_prefix = "MMC"
  )

  seu <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    assay = "RNA",
    project = "mckellardw_sim"
  )

  md <- seu@meta.data
  ncells <- ncol(seu)
  md$sample_id <- sample(sprintf("mouse_%02d", 1:6), ncells, replace = TRUE)
  md$muscle_region <- factor(sample(c("proximal", "distal", "middle"), ncells, replace = TRUE))
  md$age_weeks <- sample(c(8, 12, 16, 20), ncells, replace = TRUE)
  md$condition <- sample(c("control", "injury"), ncells, replace = TRUE, prob = c(0.65, 0.35))
  md$cell_subclass <- factor(sample(c("satellite", "myoblast", "fibroblast", "endothelial"), ncells, replace = TRUE))
  md$percent.mt <- round(stats::runif(ncells, min = 0, max = 15), 2)
  md$pct_counts_mt <- md$percent.mt
  seu@meta.data <- md
  SeuratObject::Idents(seu) <- md$cell_subclass

  counts_mat <- as.matrix(SeuratObject::GetAssayData(seu, assay = "RNA", layer = "counts"))
  data_mat <- log1p(counts_mat)
  seu <- SeuratObject::SetAssayData(seu, assay = "RNA", layer = "data", new.data = Matrix::Matrix(data_mat, sparse = TRUE))

  vfs <- sample(rownames(seu), size = 900L)
  SeuratObject::VariableFeatures(seu, assay = "RNA") <- vfs
  seu <- .add_feature_metadata(seu, assay = "RNA", variable_features = vfs, style = "mvp")
  SeuratObject::VariableFeatures(seu, assay = "RNA") <- vfs

  seu <- .add_reductions(seu, assay = "RNA", n_pcs = 20L)
  seu@commands <- list(
    NormalizeData.RNA = list(params = list(normalization.method = "LogNormalize")),
    FindVariableFeatures.RNA = list(params = list(selection.method = "mean.var.plot")),
    RunPCA.RNA = list(params = list(npcs = 20))
  )

  seu
}

.write_notes <- function(path, rds_path, rdata_path) {
  notes <- c(
    "# Simulated Inputs",
    "",
    "These mock files are intended to approximate likely Seurat structures for conversion testing.",
    "",
    "## 1) Uapinyoying-like integrated object (RDS)",
    paste0("- File: `", basename(rds_path), "`"),
    "- Assays: RNA + SCT",
    "- RNA layers: counts, data, scale.data, plus split-like layers (counts.batchA/B, data.batchA/B)",
    "- Meta.data includes: nCount_RNA, nFeature_RNA, percent.mt, pct_counts_mt, condition, batch, sample_id, celltype",
    "- Includes a list-column (`qc_notes`) to test metadata sanitization/drop behavior",
    "- Feature metadata includes VST-like columns: vst.mean, vst.variance, vst.variance.standardized, highly_variable",
    "- Reductions: pca, umap",
    "- Commands slot includes common pipeline command names and SCTransform signal",
    "",
    "## 2) Mckellard-like reference object (RData)",
    paste0("- File: `", basename(rdata_path), "`"),
    "- Contains Seurat object named `seu` (plus a couple of non-Seurat helper objects)",
    "- Assays: RNA only",
    "- RNA layers: counts + data (no scale.data by default)",
    "- Meta.data includes: nCount_RNA, nFeature_RNA, percent.mt, pct_counts_mt, muscle_region, age_weeks, condition, cell_subclass",
    "- Feature metadata includes MVP-like columns: mean, variance, dispersion, dispersion.scaled, highly_variable_rank",
    "- Reductions: pca, umap",
    "- Commands slot includes NormalizeData / FindVariableFeatures / RunPCA"
  )
  writeLines(notes, con = path, useBytes = TRUE)
}

main <- function() {
  opts <- .parse_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(opts$help)) {
    .usage()
    quit(status = 0L)
  }

  set.seed(as.integer(opts$seed))
  outdir <- opts$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  uap <- .build_uapinyoying_like()
  rds_path <- file.path(outdir, "s5f_gio.tab.uap.integrated.wt.lm_drop.rds")
  saveRDS(uap, file = rds_path)

  seu <- .build_mckellard_like()
  gene_lookup <- data.frame(
    feature_id = rownames(seu),
    gene_symbol = rownames(seu),
    stringsAsFactors = FALSE
  )
  version_tag <- "simulated_scMuscle_mm10_slim_v1-1"
  rdata_path <- file.path(outdir, "scMuscle_mm10_slim_v1-1.RData")
  save(seu, gene_lookup, version_tag, file = rdata_path)

  notes_path <- file.path(outdir, "SIMULATION_NOTES.md")
  .write_notes(notes_path, rds_path, rdata_path)

  cat("Wrote simulated files:\n")
  cat(" - ", normalizePath(rds_path), "\n", sep = "")
  cat(" - ", normalizePath(rdata_path), "\n", sep = "")
  cat(" - ", normalizePath(notes_path), "\n", sep = "")
}

tryCatch(
  main(),
  error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(status = 1L)
  }
)
