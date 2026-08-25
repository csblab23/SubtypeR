# SubtypeR

`SubtypeR` assigns a molecular subtype to every cell in a single-cell RNA-seq
dataset using [UCell](https://github.com/carmonalab/UCell) gene-signature
scoring — per cell.

Given a **raw count matrix**, `SubtypeR`:

1. Runs `Seurat::SCTransform` normalization, retaining **all genes**
   (`return.only.var.genes = FALSE`), so no signature gene is ever silently
   dropped for not being highly variable.
2. Scores a panel of gene signatures per cell with `UCell`.
3. Assigns each cell a subtype using a per-signature quantile threshold
   (default: 75th percentile):
   - passes exactly one signature -> that signature's name
   - passes more than one signature -> `"Mixed"`
   - passes none -> `"Unknown"`

## Two modes

- **Automatic** — use a built-in, curated signature panel. Currently ships
  with a high-grade serous ovarian cancer (HGSOC) panel (Differentiated /
  Immunoreactive / Mesenchymal / Proliferative). See `list_builtin_panels()`.
- **Manual** — supply your own `.gmt` file to run the same pipeline on any
  cancer type or any custom gene-signature set.

## Installation

```r
# install.packages("devtools")
devtools::install_github("csblab23/SubtypeR")
```

## Quick start

```r
library(SubtypeR)

# Automatic mode (built-in HGSOC panel)
res <- RunSubtypeR(raw_counts, mode = "auto", cancer_type = "HGSOC")
table(res$results$Subtype)

# Manual mode (any cancer type / any signatures)
res <- RunSubtypeR(raw_counts, mode = "manual", gmt_file = "my_signatures.gmt")

# QC plots
plot_subtype_scores(res$results)
plot_subtype_composition(res$results)
```

`res` is a list with:
- `seurat_object` — SCTransform-normalized Seurat object with `Subtype` and
  per-signature UCell scores added to `meta.data`
- `results` — a data.frame, one row per cell, with scores + `Subtype`
- `signatures_used` — the gene signatures actually scored, after intersecting
  with genes present in the data

See `vignette("SubtypeR-intro")` for a full walkthrough, including how to
tune the quantile threshold and how to start from an existing Seurat object.

## Citation / acknowledgment

Subtype calling logic adapted from an internal HGSOC single-cell analysis
workflow. Built on [`Seurat`](https://satijalab.org/seurat/) and
[`UCell`](https://github.com/carmonalab/UCell).
