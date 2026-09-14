#!/usr/bin/env nextflow
/*
 * MASTER PIPELINE: end-to-end tumor spatial architecture analysis
 * (patient CID44971, triple-negative breast cancer).
 *
 * Composes all 6 modules as genuine subworkflows, with actual channel
 * data flow between them -- Module 3 receives Module 1 and Module 2's
 * outputs as real channels emitted from those subworkflows, not by
 * independently reading hardcoded file paths from disk. This is what
 * distinguishes real subworkflow composition from just running six
 * pipelines in sequence and hoping the file paths line up.
 *
 * Dependency graph:
 *
 *   MODULE1 (raw scRNA-seq)         MODULE2 (raw spatial data)
 *        |                                |
 *        | final_annotated,               | investigated_spatial
 *        | markers_csv,                   |    (used by BOTH
 *        | cluster_labels                 |     Module 3 and 4)
 *        |                                |
 *        +------------> MODULE3 <---------+
 *        |            (deconvolution)     |
 *        |                  |             |
 *        |    final_deconvolved            +------------> MODULE4
 *        |                  |                             (CNV)
 *        |                  v                                |
 *        |              MODULE6                  final_cnv_annotated
 *        |         (cell-cell comm.)                          |
 *        |                                                    v
 *        |                                                MODULE5
 *        |                                          (spatial statistics)
 *
 * Usage:
 *   nextflow run main.nf
 *   (override any parameter, e.g. nextflow run main.nf --patient_id CID4465
 *    to run the Phase 2 independent-application patient)
 */

nextflow.enable.dsl = 2

include { SCRNA_MODULE1          } from './module1-scrna/main.nf'
include { SPATIAL_MODULE2        } from './module2-spatial/main.nf'
include { DECONVOLUTION_MODULE3  } from './module3-deconvolution/main.nf'
include { CNV_MODULE4            } from './module4-cnv/main.nf'
include { SPATIAL_STATS_MODULE5  } from './module5-spatial-stats/main.nf'
include { CELLCHAT_MODULE6       } from './module6-cellchat/main.nf'

workflow {

    // --- Module 1: scRNA-seq reference (no upstream dependency) ---
    scrna_matrix_ch    = Channel.fromPath(params.scrna_matrix_mtx,   checkIfExists: true)
    scrna_genes_ch     = Channel.fromPath(params.scrna_genes_tsv,    checkIfExists: true)
    scrna_barcodes_ch  = Channel.fromPath(params.scrna_barcodes_tsv, checkIfExists: true)
    scrna_metadata_ch  = Channel.fromPath(params.scrna_metadata_csv, checkIfExists: true)
    scrna_resolutions_ch = Channel.of(params.scrna_resolutions)

    SCRNA_MODULE1(
        scrna_matrix_ch, scrna_genes_ch, scrna_barcodes_ch, scrna_metadata_ch,
        params.patient_id, scrna_resolutions_ch
    )

    // --- Module 2: spatial data (no upstream dependency) ---
    spaceranger_ch        = Channel.fromPath(params.spatial_spaceranger_dir, checkIfExists: true)
    spatial_metadata_ch   = Channel.fromPath(params.spatial_metadata_csv,    checkIfExists: true)
    spatial_resolutions_ch = Channel.of(params.spatial_resolutions)

    SPATIAL_MODULE2(
        spaceranger_ch, spatial_metadata_ch, params.patient_id, spatial_resolutions_ch
    )

    // Module 2's output feeds BOTH Module 3 and Module 4 below. Without
    // .first(), this single-item channel would be exhausted by whichever
    // subworkflow consumes it first, leaving the other with nothing --
    // the same queue-channel-exhaustion lesson learned the hard way in
    // Module 1's original resolution-sweep bug.
    spatial_result = SPATIAL_MODULE2.out.investigated_spatial.first()

    // --- Module 3: deconvolution (real inputs from Modules 1 and 2) ---
    DECONVOLUTION_MODULE3(
        SCRNA_MODULE1.out.final_annotated,
        SCRNA_MODULE1.out.markers_csv,
        SCRNA_MODULE1.out.cluster_labels,
        spatial_result,
        params.patient_id
    )

    // --- Module 4: CNV inference (real input from Module 2) ---
    CNV_MODULE4(
        spatial_result,
        params.patient_id
    )

    // --- Module 5: spatial statistics (real input from Module 4) ---
    SPATIAL_STATS_MODULE5(
        CNV_MODULE4.out.final_cnv_annotated,
        params.patient_id
    )

    // --- Module 6: cell-cell communication (real input from Module 3) ---
    scalefactors_ch = Channel.fromPath(params.cellchat_scalefactors_json, checkIfExists: true)

    CELLCHAT_MODULE6(
        DECONVOLUTION_MODULE3.out.final_deconvolved,
        scalefactors_ch,
        params.patient_id
    )
}
