process RUN_CELLCHAT {
    tag "${params.patient_id}"
    publishDir "${params.cellchat_outdir}/02_cellchat", mode: 'copy'

    input:
    path cellchat_input

    output:
    path "cellchat_obj.rds", emit: cellchat_object

    script:
    """
    Rscript ${moduleDir}/../bin/run_cellchat.R ${cellchat_input} ${params.min_cells}
    """
}
