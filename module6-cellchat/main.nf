#!/usr/bin/env nextflow
/*
 * Module 6 (stretch goal): cell-cell communication analysis via
 * CellChat, using Module 3's predicted cell-type groups and real
 * spatial coordinates.
 *

 *

 *   nextflow run main.nf --patient_id CID44971
 */

nextflow.enable.dsl = 2

include { PREPARE_CELLCHAT_INPUT } from './modules/prepare_cellchat_input.nf'
include { RUN_CELLCHAT           } from './modules/run_cellchat.nf'
include { SUMMARIZE_SIGNALING    } from './modules/summarize_signaling.nf'
include { VISUALIZE_SIGNALING    } from './modules/visualize_signaling.nf'

workflow CELLCHAT_MODULE6 {
    take:
    deconvolved_object  // Module 3's final_deconvolved.rds
    scalefactors_json   // raw Space Ranger scale factors (not module-derived)
    patient_id

    main:
    PREPARE_CELLCHAT_INPUT(deconvolved_object, scalefactors_json)
    RUN_CELLCHAT(PREPARE_CELLCHAT_INPUT.out.cellchat_input)
    SUMMARIZE_SIGNALING(RUN_CELLCHAT.out.cellchat_object)
    VISUALIZE_SIGNALING(RUN_CELLCHAT.out.cellchat_object, patient_id)
}

// Default entry 
workflow {
    deconvolved_ch  = Channel.fromPath(params.deconvolved_object, checkIfExists: true)
    scalefactors_ch = Channel.fromPath(params.scalefactors_json,  checkIfExists: true)

    CELLCHAT_MODULE6(deconvolved_ch, scalefactors_ch, params.patient_id)
}
