process VISUALIZE_ON_TISSUE {
    tag "${params.patient_id}"
    publishDir "${params.spatial_outdir}/07_final", mode: 'copy'

    input:
    path best_object
    val  patient_id

    output:
    path "clusters_vs_pathology.png"
    path "final_spatial_annotated.rds"

    script:
    """
    Rscript ${moduleDir}/../bin/visualize_on_tissue.R ${best_object} ${patient_id}
    """
}
