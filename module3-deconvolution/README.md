# Module 1: scRNA-seq QC, clustering, and validated annotation

A Nextflow pipeline that takes one patient's single-cell RNA-seq data
from Wu et al. 2021 (GSE176078), performs QC filtering, sweeps four
clustering resolutions in parallel, and validates the resulting clusters
against the study's own published cell-type labels using multiple
independent lines of evidence.

## Why this exists

This is the foundation for a larger multi-omic project integrating
scRNA-seq with spatial transcriptomics (see the top-level project
README). Before using this patient's cells as a deconvolution reference
for spatial data, we needed confidence that an independently-built
clustering pipeline actually recovers real cell types — not just that it
runs without errors.

## Pipeline

```
SUBSET_PATIENT (params.patient_id)
        │
        ▼
QC_FILTER (thresholds derived from inspecting the actual data distribution)
        │
        ▼
NORMALIZE_PCA  ──► one base object (normalized, PCA, neighbor graph)
        │
        ▼
   [fan-out: Channel.of(params.resolutions)]
        │
        ▼
CLUSTER_AT_RESOLUTION  (4 parallel branches: resolutions 0.5, 0.3, 0.2, 0.1)
        │
        ▼
   [fan-in: .collect()]
        │
        ▼
COMPARE_RESOLUTIONS  (builds ARI-vs-resolution table, selects the best)
        │
        ▼
INVESTIGATE_DISAGREEMENT → FIND_MARKERS → ANNOTATE_AND_PLOT
```

## Usage

```bash
nextflow run main.nf --patient_id CID44971
```

All thresholds and settings are parameterized in `nextflow.config` —
running on a different patient (see Phase 2, `CID4465`, in the top-level
project plan) requires no code changes, only `--patient_id CID4465`.

## Data

Wu et al. 2021 (*Nature Genetics*), GEO accession GSE176078: scRNA-seq
from 26 primary breast tumors. **Raw sequencing reads are not publicly
available** for this dataset — only processed count matrices were
deposited to GEO, consistent with standard practice for human tumor
sequencing data (the same donor-privacy pattern encountered in the
companion methylation project). This pipeline works from the published,
processed matrix, same as every other paper citing this dataset.

CID44971 (triple-negative breast cancer) was chosen as the primary
validation case: it is the study's own featured/representative patient,
giving the richest published context to validate against.

## Results

**QC filtering:** 7,986 → 7,910 cells passed thresholds (`nFeature_RNA`
between 200–6000, `percent.mito` < 20%) derived from inspecting this
patient's own metric distributions, not copied from the original paper.

**Resolution sweep:**

| Resolution | Clusters | ARI vs. published `celltype_major` |
|-----------:|---------:|------------------------------------:|
| 0.5        | 18       | 0.363                                |
| 0.3        | 18       | 0.419                                |
| 0.2        | 17       | 0.422                                |
| **0.1**    | **13**   | **0.559** (selected)                 |

**Validation — why 0.559, not higher, is actually a good result:**
per-cluster purity against the published `celltype_major` labels ranged
from 96.4% to 100% — every cluster was overwhelmingly dominated by a
single published cell type. The moderate ARI is explained by
*resolution*, not misclassification: several published categories were
split across multiple clusters that each independently corresponded to
a real, finer biological subtype. For example, "T-cells" was split into:

- Cluster 0 → 94.5% CD4+ T cells
- Cluster 1 → 85.7% CD8+ T cells
- Cluster 7 → 88.6% Cycling T-cells

This was independently confirmed by marker genes: cluster 7's top
markers (`MKI67`, `RRM2`, `DLGAP5`, histone genes) are canonical,
textbook cell-cycle/proliferation genes — unambiguous biological
confirmation, not just label agreement. The same pattern held for
cancer epithelial (basal vs. cycling), fibroblast (iCAF vs. myCAF), and
endothelial sub-populations.

## Known issues and fixes

**Channel exhaustion bug (queue channel vs. value channel):** the
initial pipeline ran only 1 of 4 intended parallel resolution branches.
`NORMALIZE_PCA.out.base_object` emits a single item; when paired with
the 4-item `resolution_ch` as two separate process inputs, Nextflow
pairs channels *positionally* (like zipping two lists) by default. The
single-item channel was exhausted after pairing with the first
resolution, so the other three branches never launched — silently, with
no error. Fixed by converting the base object into a reusable **value
channel** with `.first()`, which broadcasts it to every item in the
resolution channel instead of being consumed once. This is a
well-documented Nextflow gotcha, and a genuine example of the
distinction between queue channels (consumed as they're read) and value
channels (reusable across multiple pairings).

## Structure

```
module1-scrna/
├── main.nf              <- workflow wiring, including the fan-out/fan-in
├── nextflow.config      <- all parameters (patient ID, thresholds, resolutions)
├── modules/*.nf          <- one process definition per step
├── bin/*.R               <- one standalone R script per process
└── results/              <- pipeline output (01_subset through 08_final)
```

## Requirements

- Nextflow ≥ 23.04
- R with Seurat, Matrix, dplyr, ggplot2, mclust installed and on PATH
  (this module runs R locally rather than in containers —
  `docker.enabled = false` in `nextflow.config`)
