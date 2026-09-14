process SCTRANSFORM_PCA {
    tag "${params.patient_id}"
    publishDir "${params.spatial_outdir}/03_sctransform_pca", mode: 'copy'

    input:
    path filtered_object

    output:
    path "base_spatial.rds", emit: base_object

    script:
    """
    Rscript ${moduleDir}/../bin/sctransform_pca.R ${filtered_object} ${params.spatial_npcs}
    """
}
