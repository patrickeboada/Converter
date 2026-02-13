# Simulated Inputs

These mock files are intended to approximate likely Seurat structures for conversion testing.

## 1) Uapinyoying-like integrated object (RDS)
- File: `s5f_gio.tab.uap.integrated.wt.lm_drop.rds`
- Assays: RNA + SCT
- RNA layers: counts, data, scale.data, plus split-like layers (counts.batchA/B, data.batchA/B)
- Meta.data includes: nCount_RNA, nFeature_RNA, percent.mt, pct_counts_mt, condition, batch, sample_id, celltype
- Includes a list-column (`qc_notes`) to test metadata sanitization/drop behavior
- Feature metadata includes VST-like columns: vst.mean, vst.variance, vst.variance.standardized, highly_variable
- Reductions: pca, umap
- Commands slot includes common pipeline command names and SCTransform signal

## 2) Mckellard-like reference object (RData)
- File: `scMuscle_mm10_slim_v1-1.RData`
- Contains Seurat object named `seu` (plus a couple of non-Seurat helper objects)
- Assays: RNA only
- RNA layers: counts + data (no scale.data by default)
- Meta.data includes: nCount_RNA, nFeature_RNA, percent.mt, pct_counts_mt, muscle_region, age_weeks, condition, cell_subclass
- Feature metadata includes MVP-like columns: mean, variance, dispersion, dispersion.scaled, highly_variable_rank
- Reductions: pca, umap
- Commands slot includes NormalizeData / FindVariableFeatures / RunPCA
