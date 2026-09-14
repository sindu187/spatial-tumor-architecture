#!/usr/bin/env nextflow
/*
 * Module 1: scRNA-seq QC, resolution-swept clustering, and
 * marker-validated cell-type annotation.
 *

 *
 * 
 *   nextflow run main.nf --patient_id CID44971
 */

nextflow.enable.dsl = 2

include { SUBSET_PATIENT          } from './modules/subset_patient.nf'
include { QC_FILTER                } from './modules/qc_filter.nf'
include { NORMALIZE_PCA            } from './modules/normalize_pca.nf'
include { CLUSTER_AT_RESOLUTION    } from './modules/cluster_at_resolution.nf'
include { COMPARE_RESOLUTIONS      } from './modules/compare_resolutions.nf'
include { INVESTIGATE_DISAGREEMENT } from './modules/investigate_disagreement.nf'
include { FIND_MARKERS             } from './modules/find_markers.nf'
include { ANNOTATE_AND_PLOT        } from './modules/annotate_and_plot.nf'

workflow SCRNA_MODULE1 {
    take:
    matrix_mtx
    genes_tsv
    barcodes_tsv
    metadata_csv
    patient_id
    resolutions

    main:
    SUBSET_PATIENT(matrix_mtx, genes_tsv, barcodes_tsv, metadata_csv, patient_id)
    QC_FILTER(SUBSET_PATIENT.out.raw_object)
    NORMALIZE_PCA(QC_FILTER.out.filtered_object)

    base_object_val = NORMALIZE_PCA.out.base_object.first()
    resolution_ch = resolutions.flatten()

    CLUSTER_AT_RESOLUTION(base_object_val, resolution_ch)

    all_ari_results = CLUSTER_AT_RESOLUTION.out.ari_result.collect()
    all_clustered    = CLUSTER_AT_RESOLUTION.out.clustered_object.collect()

    COMPARE_RESOLUTIONS(all_ari_results, all_clustered)
    INVESTIGATE_DISAGREEMENT(COMPARE_RESOLUTIONS.out.best_object)
    FIND_MARKERS(INVESTIGATE_DISAGREEMENT.out.investigated_object)
    ANNOTATE_AND_PLOT(FIND_MARKERS.out.markers_object, patient_id)

    emit:
    final_annotated = ANNOTATE_AND_PLOT.out.final_object
    markers_csv     = FIND_MARKERS.out.markers_csv
    cluster_labels  = ANNOTATE_AND_PLOT.out.labels_csv
}

// Default entry 
workflow {
    matrix_ch   = Channel.fromPath(params.matrix_mtx,   checkIfExists: true)
    genes_ch    = Channel.fromPath(params.genes_tsv,    checkIfExists: true)
    barcodes_ch = Channel.fromPath(params.barcodes_tsv, checkIfExists: true)
    metadata_ch = Channel.fromPath(params.metadata_csv, checkIfExists: true)
    resolutions_ch = Channel.of(params.resolutions)

    SCRNA_MODULE1(
        matrix_ch, genes_ch, barcodes_ch, metadata_ch,
        params.patient_id, resolutions_ch
    )
}
