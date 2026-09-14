#!/usr/bin/env nextflow
/*
 * Module 3: cell-type deconvolution of spatial spots, using the
 * Module 1 scRNA-seq reference. Two independent methods are run and
 * compared: anchor-based label transfer and marker-gene module scoring.
 *
 *

 *   nextflow run main.nf --patient_id CID44971
 */

nextflow.enable.dsl = 2

include { ANCHOR_LABEL_TRANSFER          } from './modules/anchor_label_transfer.nf'
include { MARKER_SCORE_DECONVOLUTION     } from './modules/marker_score_deconvolution.nf'
include { COMPARE_DECONVOLUTION_METHODS  } from './modules/compare_deconvolution_methods.nf'
include { VISUALIZE_COMPOSITION_ON_TISSUE} from './modules/visualize_composition_on_tissue.nf'

workflow DECONVOLUTION_MODULE3 {
    take:
    scrna_object    // Module 1's final_annotated.rds
    scrna_markers   // Module 1's all_cluster_markers.csv
    cluster_labels  // Module 1's cluster_labels.csv
    spatial_object  // Module 2's investigated_spatial.rds
    patient_id

    main:
    spatial_val = spatial_object.first()

    ANCHOR_LABEL_TRANSFER(scrna_object, spatial_val)
    MARKER_SCORE_DECONVOLUTION(spatial_val, scrna_markers, cluster_labels)

    COMPARE_DECONVOLUTION_METHODS(
        ANCHOR_LABEL_TRANSFER.out.predictions_csv,
        MARKER_SCORE_DECONVOLUTION.out.marker_scores
    )

    VISUALIZE_COMPOSITION_ON_TISSUE(
        ANCHOR_LABEL_TRANSFER.out.anchor_object,
        patient_id
    )

    emit:
    final_deconvolved = VISUALIZE_COMPOSITION_ON_TISSUE.out.final_object
}

// Default entry point 
workflow {
    scrna_ch   = Channel.fromPath(params.scrna_object,   checkIfExists: true)
    markers_ch = Channel.fromPath(params.scrna_markers,  checkIfExists: true)
    labels_ch  = Channel.fromPath(params.cluster_labels, checkIfExists: true)
    spatial_ch = Channel.fromPath(params.spatial_object, checkIfExists: true)

    DECONVOLUTION_MODULE3(scrna_ch, markers_ch, labels_ch, spatial_ch, params.patient_id)
}
