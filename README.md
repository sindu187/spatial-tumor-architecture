# Spatial Architecture of a Triple-Negative Breast Tumor

A six-module Nextflow pipeline integrating single-cell RNA-seq and spatial
transcriptomics to characterize the cellular composition, copy-number
architecture, and intercellular signaling of a triple-negative breast
cancer (TNBC) tumor — with quantitative validation against pathologist
ground truth at every stage.

## Overview

Each Visium spatial spot captures RNA from several co-located cells at
once, blending their identities together. This project uses a patient's
own single-cell RNA-seq data to computationally resolve that mixture,
then builds on the result to ask increasingly specific biological
questions: which cell types are where, which regions carry
copy-number-altered (malignant) signal, whether that signal is spatially
organized into coherent tumor regions, and which cell types are
signaling to each other.

All six modules are genuine Nextflow subworkflows, composed into one
master pipeline with real channel-based data flow between them — not
independent scripts invoked in sequence.

## Data

[Wu et al. 2021, *Nature Genetics*](https://www.nature.com/articles/s41588-021-00911-1)
(GEO: [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)) —
a single study providing **matched** scRNA-seq and Visium spatial
transcriptomics from the same breast tumors, including pathologist-
annotated histology labels per spatial spot. Patient **CID44971**
(triple-negative breast cancer) was used throughout, chosen for its rich
published annotation.

## Pipeline architecture

```
MODULE 1                          MODULE 2
scRNA-seq reference                spatial QC & clustering
(cell-type annotation,             (SCTransform, resolution-swept
 marker validation)                 clustering vs. pathology)
        │                                  │
        │  final_annotated,                │  investigated_spatial
        │  markers, labels                 │  (used by Modules 3 & 4)
        │                                  │
        └──────────────► MODULE 3 ◄────────┘
                    cell-type deconvolution
                 (two independent methods)
                          │              │
             final_deconvolved      (Module 2 output)
                          │              │
                    MODULE 6            MODULE 4
              cell-cell communication    CNV inference +
                 (CellChat)              pathology validation
                                              │
                                    final_cnv_annotated
                                              │
                                         MODULE 5
                                    spatial statistics
                                      (Moran's I)
```

## Key findings

- **Cell-type annotation** (Module 1) independently recovered fine-grained
  T-cell, cancer, and stromal subtypes beyond the published reference
  labels, confirmed via marker gene expression.
- **Spatial clustering** (Module 2) showed strong agreement with
  pathologist-drawn tissue regions, while resolving finer sub-regions
  within some pathology categories.
- **Deconvolution** (Module 3), benchmarked across two independent
  methods, showed strong agreement for abundant cell populations (Cancer
  Basal SC: ρ=0.82) and identified a specific limitation of anchor-based
  transfer for rarer, lower-abundance populations.
- **CNV inference** (Module 4) showed copy-number signal significantly
  elevated in pathologist-annotated malignant regions (p=9.5×10⁻⁹),
  quantitatively validated with ROC-AUC rather than visual inspection
  alone.
- **Spatial statistics** (Module 5) confirmed CNV signal is genuinely
  spatially clustered (Moran's I=0.575, permutation p=0.001) rather than
  randomly distributed.
- **Cell-cell communication** (Module 6) identified CXCL12-CXCR4 —
  a known driver of tumor-associated macrophage recruitment and an active
  TNBC therapeutic target — as a dominant significant signaling axis.

Full results, methods, and discussion are in each module's own README.

## Usage

Run the complete pipeline end to end:

```bash
nextflow run main.nf
```

Or run any module independently:

```bash
cd module1-scrna && nextflow run main.nf --patient_id CID44971
```

## Repository structure

```
.
├── main.nf, nextflow.config        <- master pipeline
├── module1-scrna/                  <- scRNA-seq QC, clustering, annotation
├── module2-spatial/                <- Visium spatial QC & clustering
├── module3-deconvolution/          <- cell-type deconvolution
├── module4-cnv/                    <- spatial CNV inference
├── module5-spatial-stats/          <- spatial statistics (Moran's I)
├── module6-cellchat/               <- cell-cell communication
└── data/raw/                       <- input data (not tracked in git)
```

Each module directory contains its own README with full methods,
results, and a repeatable Nextflow pipeline (`main.nf`, `nextflow.config`,
`modules/`, `bin/`).

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [Docker](https://www.docker.com/) (used for Module 4's CNV inference step)
- R with Seurat, CellChat, spdep, and Bioconductor annotation packages
  (see individual module READMEs for full dependency lists)

## Author

Sinduja Chandrasekaran — [github.com/sindu187](https://github.com/sindu187) · [LinkedIn](https://linkedin.com/in/sinduja.chandrasekaran)
