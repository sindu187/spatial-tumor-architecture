process COMPARE_RESOLUTIONS {
    tag "${params.patient_id}"
    publishDir "${params.scrna_outdir}/05_compare_resolutions", mode: 'copy'

    input:
    path ari_results
    path clustered_objects

    output:
    path "resolution_comparison.csv", emit: comparison_table
    path "best_clustered.rds",        emit: best_object

    script:
    """
    Rscript ${moduleDir}/../bin/compare_resolutions.R ${ari_results} -- ${clustered_objects}
    """
}
