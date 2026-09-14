process PREPARE_CNV_INPUT {
    tag "${params.patient_id}"
    publishDir "${params.cnv_outdir}/02_prepare_input", mode: 'copy'

    input:
    path spatial_object
    path gene_order

    output:
    path "counts_matrix.txt",        emit: counts_matrix
    path "gene_order_filtered.txt",  emit: gene_order_filtered
    path "annotation_file.txt",      emit: annotation_file

    script:
    """
    Rscript ${moduleDir}/../bin/prepare_cnv_input.R ${spatial_object} ${gene_order}
    """
}
