process QC_FILTER {
    tag "${params.patient_id}"
    publishDir "${params.scrna_outdir}/02_qc_filter", mode: 'copy'

    input:
    path raw_object

    output:
    path "filtered.rds", emit: filtered_object

    script:
    """
    Rscript ${moduleDir}/../bin/qc_filter.R \
        ${raw_object} ${params.scrna_min_features} ${params.scrna_max_features} ${params.scrna_max_mito}
    """
}
