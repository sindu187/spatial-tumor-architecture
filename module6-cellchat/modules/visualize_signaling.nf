process VISUALIZE_SIGNALING {
    tag "${params.patient_id}"
    publishDir "${params.cellchat_outdir}/04_final", mode: 'copy'

    input:
    path cellchat_object
    val  patient_id

    output:
    path "aggregated_communication_network.png"
    path "pathways_of_interest_bubble.png", optional: true

    script:
    """
    Rscript ${moduleDir}/../bin/visualize_signaling.R ${cellchat_object} ${patient_id}
    """
}
