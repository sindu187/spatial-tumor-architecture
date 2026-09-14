process MARKER_SCORE_DECONVOLUTION {
    tag "${params.patient_id}"
    publishDir "${params.deconvolution_outdir}/02_marker_scoring", mode: 'copy'

    input:
    path spatial_object
    path markers_csv
    path cluster_labels_csv

    output:
    path "marker_scores.csv", emit: marker_scores
    path "spatial_with_marker_scores.rds"

    script:
    """
    Rscript ${moduleDir}/../bin/marker_score_deconvolution.R \
        ${spatial_object} ${markers_csv} ${cluster_labels_csv} ${params.n_marker_genes}
    """
}
