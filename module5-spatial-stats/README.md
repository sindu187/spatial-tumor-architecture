# Module 4: Spatial CNV inference with quantitative pathology validation

A Nextflow pipeline that infers copy-number variation (CNV) directly
from spatial gene expression (Module 2's data), using pathology-
annotated normal spots as the baseline, then quantitatively tests
whether the inferred signal actually discriminates malignant from
normal tissue , but with a real statistical test and a
discrimination metric.

## The concept

Cancer cells often gain or lose large chromosomal regions. A cell with
an extra copy of a chromosome arm tends to produce proportionally more
RNA from every gene on that arm — not because any single gene became
more "active," but purely from having extra genetic template. Any one
gene's expression is too noisy to detect this reliably, but averaging
expression over many neighboring genes (ordered by real genomic
position) along a chromosome reveals consistent regional shifts that
individual-gene noise would otherwise hide. This is the principle
behind `inferCNV` and every RNA-based CNV inference method.

## Pipeline

```
BUILD_GENE_ORDER_FILE (Bioconductor annotation -> gene chr/start/end)
        │
        ▼
PREPARE_CNV_INPUT (raw counts + pathology-based reference/observation groups)
        │
        ▼
RUN_INFERCNV (Docker: trinityctat/infercnv -- see below)
        │
        ▼
VALIDATE_CNV_AGAINST_PATHOLOGY (Wilcoxon test + ROC-AUC)
        │
        ▼
VISUALIZE_CNV_ON_TISSUE
```

**Design choices:**
- Gene genomic positions are built from official Bioconductor annotation
  packages (`TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db`) rather
  than an unofficial downloaded file — fully reproducible from source,
  no unverifiable external file.
- Pathology categories ("Normal + stroma + lymphocytes", "Stroma",
  "Stroma + adipose tissue") serve as the CNV-neutral reference baseline;
  "DCIS" and "Invasive cancer + lymphocytes" are tested as the malignant
  populations.
- HMM-based CNV state calling is disabled (`HMM = FALSE`) — only the
  underlying denoised expression signal is needed for our own
  validation, and HMM adds substantial runtime without benefit here.
- `RUN_INFERCNV` runs in the official `trinityctat/infercnv` Docker
  image rather than a local install. `inferCNV` depends on JAGS (a
  system-level Bayesian sampling library, separate from R itself) which
  is notoriously difficult to install reliably on macOS. Using the
  tool maintainers' own container sidesteps this entirely, while the
  other four steps (standard R/Bioconductor packages) run locally.

## Usage

```bash
nextflow run main.nf --patient_id CID44971
```

## Results

**CNV signal is real and statistically significant, but discrimination
is modest:**

| Comparison | Wilcoxon p-value | AUC |
|---|---|---:|
| Malignant (DCIS + Invasive) vs. Normal-like | 9.5×10⁻⁹ | 0.602 |
| DCIS alone vs. Normal-like | 4.8×10⁻⁴ | 0.576 |
| Invasive cancer alone vs. Normal-like | 3.3×10⁻⁹ | 0.624 |

The p-values confirm a real, highly significant difference in CNV
signal between malignant and normal-like spots — not noise. The AUCs
(0.58–0.62, versus 0.5 = no discrimination) indicate the effect, while
real, is modest: most individual spots' CNV scores overlap substantially
between categories.

**Splitting DCIS and invasive cancer** showed slightly better
discrimination for invasive carcinoma alone (AUC 0.624 vs. 0.576 for
DCIS), consistent with invasive disease generally carrying greater
genomic instability than pre-invasive DCIS. However, the improvement
over the combined analysis (0.602) was small — histological stage
grouping is a minor contributing factor, not the dominant one.

**The dominant limitation is most likely spot-level cellular admixture** —
directly consistent with Module 3's deconvolution results, which showed
no spot is a pure population (every spot is a mixture of cancer cells,
CAFs, immune cells, etc. in varying proportions). Even spots labeled
"Invasive cancer" by a pathologist looking at tissue morphology contain
real proportions of non-malignant cells (as Module 3 quantified
directly), which dilutes the average CNV signal per spot toward
neutral. This is a fundamental resolution limit of spot-based (as
opposed to single-cell) spatial transcriptomics, not a flaw specific to
this pipeline — and it's a finding that connects directly back to
Module 3's results rather than standing alone.

## Known issues and fixes

- **Missing Bioconductor annotation packages** — `TxDb.Hsapiens.UCSC.hg38.knownGene`,
  `org.Hs.eg.db`, `pROC` needed explicit individual installation; a
  batch `BiocManager::install(c(...))` call silently stopped after the
  first package rather than installing all three, requiring them to be
  installed one at a time with success confirmed before proceeding.
- **`GetAssayData(..., slot = "counts")` is defunct** in current
  SeuratObject (≥5.0.0) — replaced with `layer = "counts"`. A clean
  one-line fix, but a reminder that pinning package versions matters for
  long-term reproducibility as dependencies evolve.

## Diagnostics

`bin/investigate_dcis_vs_invasive.R` is a standalone diagnostic (not
part of the pipeline DAG) used to test whether combining DCIS and
Invasive cancer into one "Malignant" category diluted the discrimination
signal. Kept for transparency of the investigation.

## Limitations and next steps

- **Spot-level admixture is the dominant limitation** (see Results). A
  natural extension would be restricting the CNV validation to only the
  spots Module 3 estimated as >70% cancer-cell composition, to test
  whether discrimination improves substantially when admixture is
  reduced — directly testing the admixture-dilution hypothesis rather
  than just proposing it.
- **HMM-based discrete CNV state calling** (gain/loss/neutral per
  region, rather than a continuous score) was disabled for runtime
  reasons; enabling it would allow identification of *which* specific
  chromosomal regions are altered, not just an aggregate score.

## Structure

```
module4-cnv/
├── main.nf
├── nextflow.config
├── modules/*.nf
├── bin/
│   ├── build_gene_order_file.R
│   ├── prepare_cnv_input.R
│   ├── run_infercnv.R
│   ├── validate_cnv_against_pathology.R
│   ├── visualize_cnv_on_tissue.R
│   └── investigate_dcis_vs_invasive.R   <- diagnostic, not in pipeline DAG
└── results/
```
