process VISUALIZE_SPATIAL_CLUSTERING {
    tag "${params.patient_id}"
    publishDir "${params.spatial_stats_outdir}/03_final", mode: 'copy'

    input:
    path plotting_data
    path results_csv
    path spatial_object
    val  patient_id

    output:
    path "moran_scatterplot.png"
    path "cnv_spatial_clustering.png"

    script:
    """
    Rscript ${moduleDir}/../bin/visualize_spatial_clustering.R ${plotting_data} ${results_csv} ${spatial_object} ${patient_id}
    """
}
