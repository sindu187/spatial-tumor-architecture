#!/usr/bin/env nextflow
/*
 * Module 2: Visium spatial QC, SCTransform, resolution-swept clustering,
 * and validation against pathologist-annotation.
 *
 * 
 *

 *   nextflow run main.nf --patient_id CID44971
 */

nextflow.enable.dsl = 2

include { LOAD_SPATIAL                    } from './modules/load_spatial.nf'
include { QC_FILTER_SPATIAL                } from './modules/qc_filter_spatial.nf'
include { SCTRANSFORM_PCA                  } from './modules/sctransform_pca.nf'
include { CLUSTER_AT_RESOLUTION_SPATIAL    } from './modules/cluster_at_resolution_spatial.nf'
include { COMPARE_RESOLUTIONS_SPATIAL      } from './modules/compare_resolutions_spatial.nf'
include { INVESTIGATE_DISAGREEMENT_SPATIAL } from './modules/investigate_disagreement_spatial.nf'
include { VISUALIZE_ON_TISSUE              } from './modules/visualize_on_tissue.nf'

workflow SPATIAL_MODULE2 {
    take:
    spaceranger_dir
    metadata_csv
    patient_id
    resolutions

    main:
    LOAD_SPATIAL(spaceranger_dir, metadata_csv, patient_id)
    QC_FILTER_SPATIAL(LOAD_SPATIAL.out.raw_object)
    SCTRANSFORM_PCA(QC_FILTER_SPATIAL.out.filtered_object)

    base_object_val = SCTRANSFORM_PCA.out.base_object.first()
    resolution_ch = resolutions.flatten()

    CLUSTER_AT_RESOLUTION_SPATIAL(base_object_val, resolution_ch)

    all_ari = CLUSTER_AT_RESOLUTION_SPATIAL.out.ari_result.collect()
    all_clustered = CLUSTER_AT_RESOLUTION_SPATIAL.out.clustered_object.collect()

    COMPARE_RESOLUTIONS_SPATIAL(all_ari, all_clustered)
    INVESTIGATE_DISAGREEMENT_SPATIAL(COMPARE_RESOLUTIONS_SPATIAL.out.best_object)
    VISUALIZE_ON_TISSUE(INVESTIGATE_DISAGREEMENT_SPATIAL.out.investigated_object, patient_id)

    emit:
    investigated_spatial = INVESTIGATE_DISAGREEMENT_SPATIAL.out.investigated_object
}

// Default entry 
workflow {
    spaceranger_ch = Channel.fromPath(params.spaceranger_dir, checkIfExists: true)
    metadata_ch    = Channel.fromPath(params.metadata_csv,    checkIfExists: true)
    resolutions_ch = Channel.of(params.resolutions)

    SPATIAL_MODULE2(spaceranger_ch, metadata_ch, params.patient_id, resolutions_ch)
}
