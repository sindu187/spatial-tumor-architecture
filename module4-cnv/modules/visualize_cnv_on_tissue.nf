process VISUALIZE_CNV_ON_TISSUE {
    tag "${params.patient_id}"
    publishDir "${params.cnv_outdir}/05_final", mode: 'copy'

    input:
    path spatial_object
    path cnv_scores
    val  patient_id

    output:
    path "cnv_vs_pathology.png"
    path "final_cnv_annotated.rds", emit: final_object

    script:
    """
    Rscript ${moduleDir}/../bin/visualize_cnv_on_tissue.R ${spatial_object} ${cnv_scores} ${patient_id}
    """
}
