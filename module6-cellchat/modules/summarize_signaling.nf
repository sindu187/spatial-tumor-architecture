process SUMMARIZE_SIGNALING {
    tag "${params.patient_id}"
    publishDir "${params.cellchat_outdir}/03_summary", mode: 'copy'

    input:
    path cellchat_object

    output:
    path "all_significant_interactions.csv"
    path "highlighted_interactions.csv"
    path "pathway_ranking.csv"

    script:
    """
    Rscript ${moduleDir}/../bin/summarize_signaling.R ${cellchat_object}
    """
}
