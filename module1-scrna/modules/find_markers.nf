process FIND_MARKERS {
    tag "${params.patient_id}"
    publishDir "${params.scrna_outdir}/07_markers", mode: 'copy'

    input:
    path investigated_object

    output:
    path "all_cluster_markers.csv", emit: markers_csv
    path "with_markers.rds", emit: markers_object

    script:
    """
    Rscript ${moduleDir}/../bin/find_markers.R ${investigated_object}
    """
}
