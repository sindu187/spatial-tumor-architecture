process BUILD_GENE_ORDER_FILE {
    tag "${params.patient_id}"
    publishDir "${params.cnv_outdir}/01_gene_order", mode: 'copy'

    input:
    path spatial_object

    output:
    path "gene_order.txt", emit: gene_order

    script:
    """
    Rscript ${moduleDir}/../bin/build_gene_order_file.R ${spatial_object}
    """
}
