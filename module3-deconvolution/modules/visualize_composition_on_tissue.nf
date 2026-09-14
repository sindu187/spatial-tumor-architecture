process VISUALIZE_COMPOSITION_ON_TISSUE {
    tag "${params.patient_id}"
    publishDir "${params.deconvolution_outdir}/04_final", mode: 'copy'

    input:
    path anchor_object
    val  patient_id

    output:
    path "predicted_celltype_vs_pathology.png"
    path "celltype_score_heatmaps.png"
    path "final_deconvolved.rds", emit: final_object

    script:
    """
    Rscript ${moduleDir}/../bin/visualize_composition_on_tissue.R ${anchor_object} ${patient_id}
    """
}
