#!/usr/bin/env nextflow
/*
 * Module 5: formal spatial statistics -- does CNV signal show genuine
 * spatial clustering (coherent tumor regions), or is it scattered
 * randomly, Tested with a permutation-based Moran's I test on an
 * explicit k-nearest-neighbor spatial graph built from real Visium
 * spot coordinates.
 *
 *

 *   nextflow run main.nf --patient_id CID44971
 */

nextflow.enable.dsl = 2

include { BUILD_SPATIAL_NEIGHBORS       } from './modules/build_spatial_neighbors.nf'
include { TEST_SPATIAL_AUTOCORRELATION  } from './modules/test_spatial_autocorrelation.nf'
include { VISUALIZE_SPATIAL_CLUSTERING  } from './modules/visualize_spatial_clustering.nf'

workflow SPATIAL_STATS_MODULE5 {
    take:
    spatial_object  // Module 4's final_cnv_annotated.rds
    patient_id

    main:
    spatial_val = spatial_object.first()

    BUILD_SPATIAL_NEIGHBORS(spatial_val)

    TEST_SPATIAL_AUTOCORRELATION(
        BUILD_SPATIAL_NEIGHBORS.out.neighbors,
        spatial_val
    )

    VISUALIZE_SPATIAL_CLUSTERING(
        TEST_SPATIAL_AUTOCORRELATION.out.plotting_data,
        TEST_SPATIAL_AUTOCORRELATION.out.results,
        spatial_val,
        patient_id
    )

    emit:
    moran_results = TEST_SPATIAL_AUTOCORRELATION.out.results
}

// Default entry 
workflow {
    spatial_ch = Channel.fromPath(params.spatial_object, checkIfExists: true)
    SPATIAL_STATS_MODULE5(spatial_ch, params.patient_id)
}
