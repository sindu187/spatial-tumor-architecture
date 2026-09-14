process NORMALIZE_PCA {
    tag "${params.patient_id}"
    publishDir "${params.scrna_outdir}/03_normalize_pca", mode: 'copy'

    input:
    path filtered_object

    output:
    path "base_object.rds", emit: base_object

    script:
    """
    Rscript ${moduleDir}/../bin/normalize_pca.R ${filtered_object} ${params.scrna_npcs}
    """
}
