#!/usr/bin/env Rscript

# Seurat inspection + conversion CLI
# - Reads: .rds OR .RData/.rda
# - Prints:
#     (1) cells + genes dimensions
#     (2) normalization / stats presence (counts/data/scale, HVGs, QC cols, commands, reductions)
# - Converts to: python-friendly .h5ad via SeuratDisk
#
# Requirements:
#   - SeuratObject
#   - SeuratDisk (for conversion)
# Optional:
#   - Seurat (for some command logs / convenience, but not strictly required here)

# ----------------------------- CLI plumbing -----------------------------------

.usage <- function() {
  cat(
"inspectconvert.R

Usage:
  Rscript inspectconvert.R --input <file.rds|file.RData|file.rda> [options]

Required:
  --input            Path to .rds OR .RData/.rda that contains a Seurat object.

Common options:
  --out              Output .h5ad path (default: same basename as input, .h5ad)
  --assay            Assay to export/inspect (default: RNA if present else DefaultAssay)
  --object           For .RData/.rda only: object name to load if multiple Seurat objects exist
  --overwrite        TRUE/FALSE (default: FALSE)
  --inspect-only     TRUE/FALSE (default: FALSE)  # if TRUE, only prints inspection report
  --quiet            TRUE/FALSE (default: FALSE)

Python-friendly export toggles:
  --join-layers          TRUE/FALSE (default: TRUE)   # Seurat v5: attempt JoinLayers() for split layers
  --prefer-data-as-X     TRUE/FALSE (default: TRUE)   # drop scale.data so AnnData.X uses data (not scale.data)
  --add-ident-to-obs     TRUE/FALSE (default: TRUE)   # adds seurat_ident column to obs (meta.data)
  --coerce-factors       TRUE/FALSE (default: TRUE)   # factor -> character in obs/var
  --drop-bad-metadata    TRUE/FALSE (default: TRUE)   # drop unsupported list/matrix columns in obs/var
  --keep-h5seurat        TRUE/FALSE (default: FALSE)  # keep intermediate .h5Seurat

Examples:
  Rscript inspectconvert.R --input obj.rds --out obj.h5ad --assay RNA --overwrite TRUE
  Rscript inspectconvert.R --input ws.RData --object seu --out seu.h5ad
  Rscript inspectconvert.R --input obj.rds --inspect-only TRUE
\n"
  )
}

.as_bool <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  if (x %in% c("true","t","1","yes","y","on")) return(TRUE)
  if (x %in% c("false","f","0","no","n","off")) return(FALSE)
  stop("Cannot parse boolean value from: ", x)
}

.parse_args <- function(args) {
  # supports:
  #   --key=value
  #   --key value
  #   --flag           (treated as TRUE)
  # shorthands:
  #   -i -o -a -n -h
  map_short <- list(i="input", o="out", a="assay", n="object", h="help")

  opts <- list(
    input = NULL,
    out = NULL,
    assay = NULL,
    object = NULL,
    overwrite = FALSE,
    inspect_only = FALSE,
    quiet = FALSE,
    join_layers = TRUE,
    prefer_data_as_X = TRUE,
    add_ident_to_obs = TRUE,
    coerce_factors = TRUE,
    drop_bad_metadata = TRUE,
    keep_h5seurat = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]

    if (a %in% c("-h","--help")) {
      opts$help <- TRUE
      i <- i + 1
      next
    }

    # short flags like -i value
    if (grepl("^-", a) && !grepl("^--", a)) {
      key <- sub("^-", "", a)
      if (!key %in% names(map_short)) stop("Unknown short flag: ", a)
      key <- map_short[[key]]
      if (key == "help") { opts$help <- TRUE; i <- i + 1; next }
      if (i == length(args) || grepl("^-", args[[i+1]])) {
        val <- "TRUE"
      } else {
        val <- args[[i+1]]
        i <- i + 1
      }
      opts[[key]] <- val
      i <- i + 1
      next
    }

    # long flags
    if (grepl("^--", a)) {
      if (grepl("=", a, fixed = TRUE)) {
        key <- sub("^--([^=]+)=.*$", "\\1", a)
        val <- sub("^--[^=]+=", "", a)
      } else {
        key <- sub("^--", "", a)
        if (i == length(args) || grepl("^-", args[[i+1]])) {
          val <- "TRUE"
        } else {
          val <- args[[i+1]]
          i <- i + 1
        }
      }

      # normalize keys to our option names
      key_norm <- gsub("-", "_", tolower(key))
      key2 <- switch(
        key_norm,
        "prefer_data_as_x" = "prefer_data_as_X",
        key_norm
      )

      if (!key2 %in% names(opts) && key2 != "help") {
        stop("Unknown flag: --", key)
      }
      if (key2 == "help") {
        opts$help <- TRUE
      } else {
        opts[[key2]] <- val
      }
      i <- i + 1
      next
    }

    stop("Unrecognized argument: ", a)
  }

  # coerce booleans
  bool_keys <- c("overwrite","inspect_only","quiet","join_layers","prefer_data_as_X",
                 "add_ident_to_obs","coerce_factors","drop_bad_metadata","keep_h5seurat")
  for (bk in bool_keys) {
    if (is.character(opts[[bk]])) opts[[bk]] <- .as_bool(opts[[bk]])
  }

  opts
}

# --------------------------- Seurat I/O helpers -------------------------------

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}

load_seurat_any <- function(path, object_name = NULL, verbose = TRUE) {
  .require_pkg("SeuratObject")

  if (!is.character(path) || length(path) != 1L) stop("`--input` must be a single file path.")
  if (!file.exists(path)) stop("Input file not found: ", path)

  ext <- tolower(tools::file_ext(path))

  if (ext == "rds") {
    obj <- readRDS(path)
    if (!inherits(obj, "Seurat")) stop("RDS did not contain a Seurat object. Class: ", paste(class(obj), collapse=", "))
    if (verbose) message("Loaded Seurat object from RDS: ", path)
    return(obj)
  }

  if (ext %in% c("rdata","rda")) {
    env <- new.env(parent = emptyenv())
    obj_names <- load(path, envir = env)

    if (length(obj_names) == 0) stop("No objects found in: ", path)

    # If user specified object name, use it
    if (!is.null(object_name)) {
      if (!object_name %in% obj_names) {
        stop("Object '", object_name, "' not found in ", path, ". Available: ", paste(obj_names, collapse=", "))
      }
      obj <- get(object_name, envir = env)
      if (!inherits(obj, "Seurat")) {
        stop("Object '", object_name, "' is not a Seurat object. Class: ", paste(class(obj), collapse=", "))
      }
      if (verbose) message("Loaded Seurat object '", object_name, "' from: ", path)
      return(obj)
    }

    # Otherwise: auto-detect Seurat objects
    objs <- mget(obj_names, envir = env, ifnotfound = NA)
    is_seu <- vapply(objs, function(x) inherits(x, "Seurat"), logical(1))
    seu_names <- names(objs)[is_seu]

    if (length(seu_names) == 1) {
      obj <- objs[[seu_names]]
      if (verbose) message("Auto-selected Seurat object '", seu_names, "' from: ", path)
      return(obj)
    }

    if (length(seu_names) == 0) {
      stop("No Seurat object found in ", path, ". Objects present: ", paste(obj_names, collapse=", "))
    }

    stop(
      "Multiple Seurat objects found in ", path, ": ", paste(seu_names, collapse=", "),
      "\nSpecify one with --object <name>."
    )
  }

  stop("Unsupported input extension: .", ext, " (supported: .rds, .RData, .rda)")
}

.get_feature_metadata_df <- function(assay_obj) {
  # Seurat v4 Assay has @meta.features; Seurat v5 Assay5 often uses @meta.data
  if (!isS4(assay_obj)) return(NULL)
  sn <- methods::slotNames(assay_obj)
  if ("meta.features" %in% sn) return(assay_obj@meta.features)
  if ("meta.data" %in% sn) return(assay_obj@meta.data)
  NULL
}

.set_feature_metadata_df <- function(assay_obj, df) {
  if (!isS4(assay_obj)) return(assay_obj)
  sn <- methods::slotNames(assay_obj)
  if ("meta.features" %in% sn) {
    assay_obj@meta.features <- df
    return(assay_obj)
  }
  if ("meta.data" %in% sn) {
    assay_obj@meta.data <- df
    return(assay_obj)
  }
  assay_obj
}

.get_assay_data_compat <- function(seu, assay, name = c("counts","data","scale.data")) {
  .require_pkg("SeuratObject")
  name <- match.arg(name)

  # Seurat v5 prefers layer= ; v4 uses slot=
  mat <- tryCatch(SeuratObject::GetAssayData(seu, assay = assay, layer = name), error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  mat <- tryCatch(SeuratObject::GetAssayData(seu, assay = assay, slot = name), error = function(e) NULL)
  mat
}

.dim_or_na <- function(x) {
  if (is.null(x)) return(c(NA_integer_, NA_integer_))
  d <- dim(x)
  if (is.null(d)) c(NA_integer_, NA_integer_) else d
}

.fmt_dim <- function(d) {
  if (any(is.na(d))) return("NA")
  paste0(d[1], " x ", d[2])
}

# --------------------------- Inspection report --------------------------------

inspect_seurat_object <- function(seu, assay = NULL, quiet = FALSE) {
  .require_pkg("SeuratObject")

  assays <- SeuratObject::Assays(seu)
  def_assay <- SeuratObject::DefaultAssay(seu)

  if (is.null(assay)) {
    assay <- if ("RNA" %in% assays) "RNA" else def_assay
  }
  if (!assay %in% assays) {
    stop("Assay '", assay, "' not found. Available: ", paste(assays, collapse=", "))
  }
  SeuratObject::DefaultAssay(seu) <- assay

  ncells <- ncol(seu)

  # Determine feature count for selected assay using counts->data->scale fallback
  d_counts <- .dim_or_na(.get_assay_data_compat(seu, assay, "counts"))
  d_data   <- .dim_or_na(.get_assay_data_compat(seu, assay, "data"))
  d_scale  <- .dim_or_na(.get_assay_data_compat(seu, assay, "scale.data"))

  nfeat <- if (!any(is.na(d_counts))) d_counts[1] else if (!any(is.na(d_data))) d_data[1] else d_scale[1]

  # Pipeline “presence” flags
  has_counts <- !any(is.na(d_counts)) && all(d_counts > 0)
  has_data   <- !any(is.na(d_data))   && all(d_data > 0)
  has_scale  <- !any(is.na(d_scale))  && all(d_scale > 0)

  # Variable features
  var_feats <- tryCatch(SeuratObject::VariableFeatures(seu, assay = assay), error = function(e) character(0))
  n_varfeat <- length(var_feats)

  # Cell metadata QC columns
  md <- seu[[]]
  md_cols <- colnames(md)
  qc_candidates <- c(
    paste0("nCount_", assay),
    paste0("nFeature_", assay),
    "percent.mt", "percent_mito", "pct_counts_mt",
    "nCount_RNA","nFeature_RNA"
  )
  qc_present <- intersect(qc_candidates, md_cols)

  # Feature metadata stats columns
  feat_md <- .get_feature_metadata_df(seu[[assay]])
  feat_cols <- if (!is.null(feat_md)) colnames(feat_md) else character(0)
  hvg_stat_candidates <- c(
    "vst.mean","vst.variance","vst.variance.standardized",
    "mvp.mean","mvp.dispersion","mvp.dispersion.scaled",
    "mean","variance","dispersion","dispersion.scaled",
    "highly_variable","highly_variable_rank"
  )
  hvg_stats_present <- intersect(hvg_stat_candidates, feat_cols)

  # Reductions
  reductions <- tryCatch(names(seu@reductions), error = function(e) character(0))

  # Commands log (best-effort)
  cmd_names <- tryCatch(names(seu@commands), error = function(e) character(0))
  has_cmd <- length(cmd_names) > 0

  # Heuristics for “was normalization likely run?”
  # (Not perfect; we’re reporting presence/absence signals.)
  norm_signals <- list(
    normalized_matrix_present = has_data,
    scaled_matrix_present = has_scale,
    find_variable_features_present = (n_varfeat > 0) || any(grepl("FindVariableFeatures", cmd_names, ignore.case = TRUE)),
    normalize_data_command_logged = any(grepl("NormalizeData", cmd_names, ignore.case = TRUE)),
    sctransform_command_logged = any(grepl("SCTransform", cmd_names, ignore.case = TRUE)),
    sct_assay_present = ("SCT" %in% assays),
    pca_present = ("pca" %in% tolower(reductions)),
    umap_present = ("umap" %in% tolower(reductions))
  )

  if (!quiet) {
    cat("============================================================\n")
    cat("Seurat inspection report\n")
    cat("------------------------------------------------------------\n")
    cat("Assays            :", paste(assays, collapse = ", "), "\n")
    cat("Default assay     :", def_assay, "\n")
    cat("Selected assay    :", assay, "\n")
    cat("Cells (ncol)      :", ncells, "\n")
    cat("Genes/features    :", nfeat, "\n")
    cat("------------------------------------------------------------\n")

    cat("Matrix dimensions (features x cells) for assay '", assay, "':\n", sep="")
    cat("  counts     :", .fmt_dim(d_counts), "\n")
    cat("  data       :", .fmt_dim(d_data), "\n")
    cat("  scale.data :", .fmt_dim(d_scale), "\n")

    # Seurat v5 layers (if available)
    layers <- tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) character(0))
    if (length(layers) > 0) {
      cat("  layers     :", paste(layers, collapse = ", "), "\n")
    }

    cat("------------------------------------------------------------\n")
    cat("Normalization / statistics signals (presence/absence):\n")
    cat("  - counts present                 :", has_counts, "\n")
    cat("  - normalized 'data' present      :", has_data, "\n")
    cat("  - scaled 'scale.data' present    :", has_scale, "\n")
    cat("  - VariableFeatures count         :", n_varfeat, "\n")
    cat("  - QC columns in meta.data        :", if (length(qc_present)) paste(qc_present, collapse=", ") else "none detected", "\n")
    cat("  - HVG/stat columns in var        :", if (length(hvg_stats_present)) paste(hvg_stats_present, collapse=", ") else "none detected", "\n")
    cat("  - Reductions present             :", if (length(reductions)) paste(reductions, collapse=", ") else "none", "\n")

    if (has_cmd) {
      cat("  - Commands logged (obj@commands) :", length(cmd_names), "\n")
      # Show a compact subset of “pipeline-ish” commands
      pipe_cmds <- cmd_names[grepl("NormalizeData|SCTransform|FindVariableFeatures|ScaleData|RunPCA|RunUMAP|RunTSNE", cmd_names, ignore.case = TRUE)]
      if (length(pipe_cmds)) {
        cat("    pipeline commands detected     :", paste(pipe_cmds, collapse=", "), "\n")
      } else {
        cat("    pipeline commands detected     : none matched common names\n")
      }
    } else {
      cat("  - Commands logged (obj@commands) : none / unavailable\n")
    }

    cat("------------------------------------------------------------\n")
    cat("Interpretation (heuristic):\n")
    if (norm_signals$sctransform_command_logged || norm_signals$sct_assay_present) {
      cat("  * SCTransform signal: YES (SCT assay and/or command log)\n")
    } else {
      cat("  * SCTransform signal: NO\n")
    }
    if (norm_signals$normalize_data_command_logged || has_data) {
      cat("  * NormalizeData signal: YES (data matrix and/or command log)\n")
    } else {
      cat("  * NormalizeData signal: NO (data matrix missing and no NormalizeData command log)\n")
    }
    if (norm_signals$find_variable_features_present) {
      cat("  * HVG/FindVariableFeatures signal: YES (VariableFeatures and/or command log)\n")
    } else {
      cat("  * HVG/FindVariableFeatures signal: NO\n")
    }
    cat("============================================================\n\n")
  }

  invisible(list(
    assay = assay,
    ncells = ncells,
    nfeatures = nfeat,
    dims = list(counts = d_counts, data = d_data, scale.data = d_scale),
    layers = tryCatch(SeuratObject::Layers(seu[[assay]]), error = function(e) character(0)),
    n_variable_features = n_varfeat,
    qc_columns_present = qc_present,
    hvg_stat_columns_present = hvg_stats_present,
    reductions = reductions,
    commands = cmd_names,
    signals = norm_signals
  ))
}

# --------------------------- Export helpers -----------------------------------

.sanitize_df_for_hdf5 <- function(df,
                                  name = "data.frame",
                                  coerce_factors = TRUE,
                                  drop_unsupported = TRUE,
                                  quiet = FALSE) {
  if (is.null(df) || !is.data.frame(df)) return(df)

  out <- df
  bad_cols <- character(0)

  for (col in names(out)) {
    x <- out[[col]]

    # Factor -> character (best for Python/Scanpy categories)
    if (is.factor(x) && coerce_factors) {
      out[[col]] <- as.character(x)
      next
    }

    # Dates/times -> character
    if (inherits(x, c("POSIXct","POSIXlt","Date"))) {
      out[[col]] <- as.character(x)
      next
    }

    # Lists / nested objects -> drop (default) to avoid HDF5 serialization problems
    if (is.list(x) || is.matrix(x) || is.data.frame(x)) {
      bad_cols <- c(bad_cols, col)
      next
    }

    # Non-atomic types -> try stringify else drop
    if (!(is.atomic(x) || is.character(x))) {
      try_chr <- tryCatch(as.character(x), error = function(e) NULL)
      if (is.null(try_chr) || length(try_chr) != nrow(out)) {
        bad_cols <- c(bad_cols, col)
      } else {
        out[[col]] <- try_chr
      }
    }
  }

  if (length(bad_cols) > 0) {
    msg <- paste0("Unsupported ", name, " columns: ", paste(bad_cols, collapse = ", "))
    if (drop_unsupported) {
      if (!quiet) message(msg, " -> dropping for h5ad compatibility.")
      out[bad_cols] <- NULL
    } else {
      stop(msg, " -> set --drop-bad-metadata TRUE or remove them manually.")
    }
  }

  out
}

.maybe_join_layers <- function(seu, assay, quiet = FALSE) {
  .require_pkg("SeuratObject")
  a <- seu[[assay]]

  # JoinLayers exists in SeuratObject v5+; do nothing if unavailable
  if (!"JoinLayers" %in% getNamespaceExports("SeuratObject")) return(seu)

  lay <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  if (length(lay) == 0) return(seu)

  looks_split <- length(lay) > 3 || any(grepl("\\.", lay))
  if (looks_split) {
    if (!quiet) message("Joining layers in assay '", assay, "' (", length(lay), " layers)...")
    a2 <- tryCatch(SeuratObject::JoinLayers(a), error = function(e) a)
    seu[[assay]] <- a2
  }
  seu
}

.drop_scale_data_for_export <- function(seu, assay, quiet = FALSE) {
  a <- seu[[assay]]

  # Seurat v5: delete layer
  lay <- tryCatch(SeuratObject::Layers(a), error = function(e) character(0))
  if ("scale.data" %in% lay) {
    if (!quiet) message("Dropping layer 'scale.data' from assay '", assay, "' for export (so AnnData.X uses data).")
    a$scale.data <- NULL
    seu[[assay]] <- a
    return(seu)
  }

  # Seurat v4: clear slot if present
  if (isS4(a) && "scale.data" %in% methods::slotNames(a)) {
    if (!quiet) message("Clearing slot '@scale.data' from assay '", assay, "' for export.")
    a@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
    seu[[assay]] <- a
  }

  seu
}

seurat_to_h5ad <- function(seu,
                           out_h5ad,
                           assay = NULL,
                           overwrite = FALSE,
                           join_layers = TRUE,
                           prefer_data_as_X = TRUE,
                           add_ident_to_obs = TRUE,
                           coerce_factors = TRUE,
                           drop_bad_metadata = TRUE,
                           keep_h5seurat = FALSE,
                           quiet = FALSE) {
  .require_pkg("SeuratObject")
  .require_pkg("SeuratDisk")

  assays <- SeuratObject::Assays(seu)
  if (is.null(assay)) assay <- if ("RNA" %in% assays) "RNA" else SeuratObject::DefaultAssay(seu)
  if (!assay %in% assays) stop("Assay '", assay, "' not present. Available: ", paste(assays, collapse=", "))

  SeuratObject::DefaultAssay(seu) <- assay

  if (is.null(out_h5ad) || !nzchar(out_h5ad)) stop("Output path --out must be provided or derivable.")
  out_dir <- dirname(out_h5ad)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(out_h5ad) && !overwrite) stop("Output exists: ", out_h5ad, " (set --overwrite TRUE)")
  h5seurat_path <- sub("\\.h5ad$", ".h5Seurat", out_h5ad, ignore.case = TRUE)

  if (file.exists(h5seurat_path) && !overwrite) {
    stop("Intermediate exists: ", h5seurat_path, " (set --overwrite TRUE or change --out)")
  }

  exp <- seu

  if (join_layers) exp <- .maybe_join_layers(exp, assay = assay, quiet = quiet)

  if (add_ident_to_obs) {
    exp$seurat_ident <- as.character(SeuratObject::Idents(exp))
  }

  # Sanitize obs/var metadata for h5ad friendliness
  if (coerce_factors || drop_bad_metadata) {
    exp@meta.data <- .sanitize_df_for_hdf5(
      exp@meta.data,
      name = "meta.data",
      coerce_factors = coerce_factors,
      drop_unsupported = drop_bad_metadata,
      quiet = quiet
    )

    for (a in SeuratObject::Assays(exp)) {
      md <- .get_feature_metadata_df(exp[[a]])
      if (!is.null(md) && is.data.frame(md)) {
        md2 <- .sanitize_df_for_hdf5(
          md,
          name = paste0("feature-metadata(", a, ")"),
          coerce_factors = coerce_factors,
          drop_unsupported = drop_bad_metadata,
          quiet = quiet
        )
        exp[[a]] <- .set_feature_metadata_df(exp[[a]], md2)
      }
    }
  }

  # Ensure AnnData.X ends up as "data" (not scale.data)
  if (prefer_data_as_X) {
    exp <- .drop_scale_data_for_export(exp, assay = assay, quiet = quiet)
  }

  if (!quiet) {
    message("Starting conversion: Seurat -> h5Seurat -> h5ad")
    message("  h5Seurat intermediate: ", h5seurat_path)
    message("  h5ad output          : ", out_h5ad)
  }

  convert_once <- function(obj) {
    SeuratDisk::SaveH5Seurat(obj, filename = h5seurat_path, overwrite = overwrite, verbose = !quiet)
    SeuratDisk::Convert(h5seurat_path, dest = "h5ad", assay = assay, overwrite = overwrite, verbose = !quiet)

    auto_h5ad <- sub("\\.h5Seurat$", ".h5ad", h5seurat_path, ignore.case = TRUE)
    if (file.exists(auto_h5ad) && normalizePath(auto_h5ad) != normalizePath(out_h5ad)) {
      if (file.exists(out_h5ad) && overwrite) file.remove(out_h5ad)
      file.rename(auto_h5ad, out_h5ad)
    }

    if (!file.exists(out_h5ad)) stop("Conversion completed but output not found: ", out_h5ad)
    out_h5ad
  }

  # Try conversion; if it fails due to Assay5 quirks, retry after coercing to v3 Assay
  out <- tryCatch(
    convert_once(exp),
    error = function(e) {
      aobj <- exp[[assay]]
      if (inherits(aobj, "Assay5") || inherits(aobj, "StdAssay")) {
        if (!quiet) message("Retrying conversion after casting assay '", assay, "' to v3 'Assay' class...")
        if (overwrite) {
          if (file.exists(h5seurat_path)) try(file.remove(h5seurat_path), silent = TRUE)
          if (file.exists(out_h5ad)) try(file.remove(out_h5ad), silent = TRUE)
        }
        exp2 <- exp
        exp2[[assay]] <- as(object = aobj, Class = "Assay")
        SeuratObject::DefaultAssay(exp2) <- assay
        convert_once(exp2)
      } else {
        stop(e)
      }
    }
  )

  if (!keep_h5seurat && file.exists(h5seurat_path)) {
    file.remove(h5seurat_path)
  }

  if (!quiet) message("Wrote h5ad: ", out)
  invisible(out)
}

# --------------------------------- MAIN ---------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) { .usage(); quit(status = 1) }

  opts <- .parse_args(args)
  if (!is.null(opts$help) && isTRUE(opts$help)) { .usage(); quit(status = 0) }

  if (is.null(opts$input)) {
    stop("Missing required --input <file.rds|file.RData|file.rda>")
  }

  # Default output path if not specified
  if (is.null(opts$out)) {
    base <- basename(opts$input)
    base <- sub("\\.(rds|rda|rdata)$", "", base, ignore.case = TRUE)
    opts$out <- file.path(dirname(opts$input), paste0(base, ".h5ad"))
  }

  seu <- load_seurat_any(opts$input, object_name = opts$object, verbose = !opts$quiet)

  # Inspection always runs
  rep <- inspect_seurat_object(seu, assay = opts$assay, quiet = opts$quiet)

  # Convert unless inspect-only
  if (!opts$inspect_only) {
    seurat_to_h5ad(
      seu = seu,
      out_h5ad = opts$out,
      assay = rep$assay,
      overwrite = opts$overwrite,
      join_layers = opts$join_layers,
      prefer_data_as_X = opts$prefer_data_as_X,
      add_ident_to_obs = opts$add_ident_to_obs,
      coerce_factors = opts$coerce_factors,
      drop_bad_metadata = opts$drop_bad_metadata,
      keep_h5seurat = opts$keep_h5seurat,
      quiet = opts$quiet
    )
  } else {
    if (!opts$quiet) message("inspect-only=TRUE, skipping conversion.")
  }
}

tryCatch(
  main(),
  error = function(e) {
    message("ERROR: ", conditionMessage(e))
    message("Run with --help to see usage.")
    quit(status = 1)
  }
)
