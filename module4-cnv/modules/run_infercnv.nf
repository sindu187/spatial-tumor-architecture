process RUN_INFERCNV {
    tag "${params.patient_id}"
    publishDir "${params.cnv_outdir}/03_infercnv", mode: 'copy'
    container 'trinityctat/infercnv:latest'

    input:
    path counts_matrix
    path annotation_file
    path gene_order_filtered
    path script_file  // staged explicitly as a file input so it is
                        // correctly mounted inside the Docker container --
                        // moduleDir-based host paths (used by every other,
                        // non-containerized process) are not automatically
                        // visible across the container boundary.

    output:
    path "final_infercnv_obj.rds", emit: infercnv_object
    path "infercnv_output/*",      optional: true

    script:
    """
    Rscript ${script_file} \
        ${counts_matrix} ${annotation_file} ${gene_order_filtered} \
        "${params.ref_groups}" ${params.cnv_cutoff} ${params.num_threads}
    """
}
