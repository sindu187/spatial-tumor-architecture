#!/usr/bin/env nextflow
/*
 * Module 4: spatial CNV inference (inferCNV) with quantitative
 * validation against pathologist-annotated malignant/normal regions.
 *
 * 
 *
 * 
 *   nextflow run main.nf --patient_id CID44971
 */

nextflow.enable.dsl = 2

include { BUILD_GENE_ORDER_FILE          } from './modules/build_gene_order_file.nf'
include { PREPARE_CNV_INPUT               } from './modules/prepare_cnv_input.nf'
include { RUN_INFERCNV                    } from './modules/run_infercnv.nf'
include { VALIDATE_CNV_AGAINST_PATHOLOGY  } from './modules/validate_cnv_against_pathology.nf'
include { VISUALIZE_CNV_ON_TISSUE         } from './modules/visualize_cnv_on_tissue.nf'

workflow CNV_MODULE4 {
    take:
    spatial_object  // Module 2's investigated_spatial.rds
    patient_id

    main:
    spatial_val = spatial_object.first()

    // Explicitly staged so it is visible inside RUN_INFERCNV's Docker
    // container (see note above and in run_infercnv.nf).
    infercnv_script_ch = Channel.fromPath("${moduleDir}/bin/run_infercnv.R", checkIfExists: true)

    BUILD_GENE_ORDER_FILE(spatial_val)
    PREPARE_CNV_INPUT(spatial_val, BUILD_GENE_ORDER_FILE.out.gene_order)

    RUN_INFERCNV(
        PREPARE_CNV_INPUT.out.counts_matrix,
        PREPARE_CNV_INPUT.out.annotation_file,
        PREPARE_CNV_INPUT.out.gene_order_filtered,
        infercnv_script_ch
    )

    VALIDATE_CNV_AGAINST_PATHOLOGY(
        RUN_INFERCNV.out.infercnv_object,
        PREPARE_CNV_INPUT.out.annotation_file
    )

    VISUALIZE_CNV_ON_TISSUE(
        spatial_val,
        VALIDATE_CNV_AGAINST_PATHOLOGY.out.cnv_scores,
        patient_id
    )

    emit:
    final_cnv_annotated = VISUALIZE_CNV_ON_TISSUE.out.final_object
}

// Default entry 
workflow {
    spatial_ch = Channel.fromPath(params.spatial_object, checkIfExists: true)
    CNV_MODULE4(spatial_ch, params.patient_id)
}
