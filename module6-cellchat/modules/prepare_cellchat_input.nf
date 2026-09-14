process PREPARE_CELLCHAT_INPUT {
    tag "${params.patient_id}"
    publishDir "${params.cellchat_outdir}/01_prepare", mode: 'copy'

    input:
    path deconvolved_object
    path scalefactors_json

    output:
    path "cellchat_input.rds", emit: cellchat_input

    script:
    """
    Rscript ${moduleDir}/../bin/prepare_cellchat_input.R ${deconvolved_object} ${scalefactors_json}
    """
}
