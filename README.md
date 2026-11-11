# Spatial Architecture of a Triple-Negative Breast Tumor

A six-module Nextflow pipeline that integrates single-cell RNA-seq and Visium spatial transcriptomics from one triple-negative breast cancer (TNBC) patient. It resolves cellular composition, copy-number architecture and cell-cell signaling, and checks the results against pathologist annotations at each stage.

The modules are Nextflow subworkflows composed into one master pipeline. Data passes between them through channels.

## Data

[Wu et al. 2021, *Nature Genetics*](https://www.nature.com/articles/s41588-021-00911-1) (GEO: [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)). It provides matched scRNA-seq and Visium data from the same tumors, with a pathologist label for each spot. Patient **CID44971** is used throughout. Only processed count matrices are public for this dataset, so raw reads are not used.

## Pipeline

```
MODULE 1                          MODULE 2
scRNA-seq reference               spatial QC & clustering
        │                                 │
        └──────────────► MODULE 3 ◄───────┘
                  cell-type deconvolution
                     │              │
                     ▼              ▼
                 MODULE 6       MODULE 4
              cell-cell comm.   CNV inference
               (CellChat)           │
                                    ▼
                                MODULE 5
                            spatial statistics
```

## Modules and key results

| # | Module | What it does | Key result |
|---|--------|--------------|------------|
| 1 | `module1-scrna` | QC, clustering at four resolutions in parallel, marker-based annotation | 7,910 cells kept after QC. Best resolution 0.1 (13 clusters, ARI 0.559). Cluster purity against published labels 96–100%. |
| 2 | `module2-spatial` | Spot QC, SCTransform, resolution sweep scored against pathology | Best resolution 0.5, ARI 0.521 against pathologist regions. |
| 3 | `module3-deconvolution` | Per-spot cell-type composition, using anchor-based label transfer and marker scoring | The two methods agree for abundant types (Cancer Basal SC, Spearman ρ = 0.82) and less for rare types (lowest ρ = 0.38). |
| 4 | `module4-cnv` | inferCNV on spatial data (Docker), validated against pathology | CNV signal is higher in malignant regions (Wilcoxon p = 9.5×10⁻⁹, AUC 0.602). Discrimination is modest because each spot mixes several cell types. |
| 5 | `module5-spatial-stats` | Moran's I on the CNV signal | CNV signal is spatially clustered (Moran's I = 0.575, permutation p = 0.001). |
| 6 | `module6-cellchat` | Spatial-mode CellChat ligand-receptor inference | 1,426 significant interactions. CXCL12–CXCR4 is a dominant axis, with the strongest signal in macrophages. |

## Usage

```bash
# full pipeline
nextflow run main.nf

# single module
cd module1-scrna && nextflow run main.nf --patient_id CID44971
```

Parameters (patient ID, QC thresholds, resolutions) are set in each `nextflow.config`.

## Repository structure

```
.
├── main.nf, nextflow.config   master pipeline
├── module1-scrna/ … module6-cellchat/
│   ├── main.nf, nextflow.config
│   ├── modules/               one process per step
│   └── bin/                   one R script per process
├── results/                   summary tables and figures
└── data/raw/                  input data (not tracked)
```

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 23.04
- [Docker](https://www.docker.com/) (Module 4's inferCNV step only)
- R with Seurat, mclust, spdep, CellChat and Bioconductor annotation packages

## Author

Sinduja Chandrasekaran — [github.com/sindu187](https://github.com/sindu187) · [LinkedIn](https://linkedin.com/in/sinduja.chandrasekaran)
