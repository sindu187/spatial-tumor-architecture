process TEST_SPATIAL_AUTOCORRELATION {
    tag "${params.patient_id}"
    publishDir "${params.spatial_stats_outdir}/02_moran_test", mode: 'copy'

    input:
    path neighbors
    path spatial_object

    output:
    path "moran_test_results.csv",        emit: results
    path "moran_data_for_plotting.rds",   emit: plotting_data

    script:
    """
    Rscript ${moduleDir}/../bin/test_spatial_autocorrelation.R ${neighbors} ${spatial_object} ${params.n_permutations}
    """
}
