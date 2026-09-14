process COMPARE_RESOLUTIONS_SPATIAL {
    tag "${params.patient_id}"
    publishDir "${params.spatial_outdir}/05_compare_resolutions", mode: 'copy'

    input:
    path ari_results
    path clustered_objects

    output:
    path "resolution_comparison_spatial.csv", emit: comparison_table
    path "best_clustered_spatial.rds",        emit: best_object

    script:
    """
    Rscript ${moduleDir}/../bin/compare_resolutions_spatial.R ${ari_results} -- ${clustered_objects}
    """
}
