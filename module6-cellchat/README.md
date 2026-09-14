# Module 6: Cell-cell communication analysis (CellChat)

A Nextflow pipeline inferring ligand-receptor signaling between the
cell-type groups identified in Module 3, using CellChat's spatial mode
(signaling constrained by real physical distance between spots, not
just co-expression regardless of location).

## The question

Modules 1–5 established what cell types are present, where they sit
spatially, which regions carry CNV signal, and that signal is spatially
clustered. This module asks a different question: **are these cell
types actually communicating** with each other via known ligand-receptor
signaling — and specifically, do any pathways relevant to this TNBC
tumor's biology (immune checkpoint signaling, macrophage recruitment,
tumor-stroma crosstalk) show up as significant?

## Pipeline

```
PREPARE_CELLCHAT_INPUT (normalized expression + Module 3's predicted
                         cell-type groups + spatial coordinates + Space
                         Ranger scale factors)
        │
        ▼
RUN_CELLCHAT (spatial-mode CellChat: distance-constrained ligand-
              receptor inference, CellChatDB.human)
        │
        ▼
SUMMARIZE_SIGNALING (full results + curated pathways of interest)
        │
        ▼
VISUALIZE_SIGNALING (aggregated network + pathway-specific plots)
```

## Usage

```bash
nextflow run main.nf --patient_id CID44971
```

## Results

**1,426 significant ligand-receptor interactions identified**, using
the 7 cell-type groups Module 3's deconvolution actually assigned as a
dominant spot identity (CAFs, Cancer Basal SC, Endothelial, Luminal
Progenitors, Macrophage, Myoepithelial, PVL Immature).

**Overall signaling landscape** is dominated by extracellular matrix
pathways (COLLAGEN, LAMININ, FN1, THBS) — consistent with a
desmoplastic tumor stroma, and a reasonable sanity check that the
method recovers biologically sensible structure rather than noise.

**CXCL12-CXCR4 emerged as a dominant, highly significant signaling
axis** (communication probability up to 0.038, most interactions
p = 0), most strongly among macrophages (Macrophage→Macrophage:
prob = 0.038) and between macrophages, CAFs, endothelial cells, and
cancer cells broadly. This is a well-established biological finding —
CXCL12/CXCR4 signaling is known to recruit and retain tumor-associated
macrophages, and is an active therapeutic target in TNBC — recovered
here from an independently built pipeline rather than assumed.

**The PD-L1/PD-1 immune checkpoint axis did not appear among
significant results.** This is consistent with, not contradictory to,
Module 3's finding that anchor-based deconvolution assigned no spot a
dominant T-cell identity in this sample — PDCD1 (PD-1) is a T-cell
receptor, so its absence here reflects the same underlying limitation
already identified and explained in Module 3, rather than a new,
unexplained gap.

## Structure

```
module6-cellchat/
├── main.nf
├── nextflow.config
├── modules/*.nf
├── bin/
│   ├── prepare_cellchat_input.R
│   ├── run_cellchat.R
│   ├── summarize_signaling.R
│   └── visualize_signaling.R
└── results/
```

## Requirements

- R with `CellChat` (jinworks/CellChat), `Seurat`, `jsonlite`
- CellChat has a substantial dependency chain (NMF, circlize,
  ComplexHeatmap, and several CRAN packages with compiled system
  dependencies); a pre-built conda binary is unavailable for Intel
  macOS, so installation here was via `remotes::install_github()` with
  several dependencies resolved individually.
