process BUILD_SPATIAL_NEIGHBORS {
    tag "${params.patient_id}"
    publishDir "${params.spatial_stats_outdir}/01_neighbors", mode: 'copy'

    input:
    path spatial_object

    output:
    path "spatial_neighbors.rds", emit: neighbors

    script:
    """
    Rscript ${moduleDir}/../bin/build_spatial_neighbors.R ${spatial_object} ${params.k_neighbors}
    """
}
