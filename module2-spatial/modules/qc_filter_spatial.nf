process QC_FILTER_SPATIAL {
    tag "${params.patient_id}"
    publishDir "${params.spatial_outdir}/02_qc_filter", mode: 'copy'

    input:
    path raw_object

    output:
    path "filtered_spatial.rds", emit: filtered_object

    script:
    """
    Rscript ${moduleDir}/../bin/qc_filter_spatial.R ${raw_object} ${params.spatial_min_features} ${params.spatial_max_mito}
    """
}
