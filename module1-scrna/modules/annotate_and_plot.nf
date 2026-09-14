process ANNOTATE_AND_PLOT {
    tag "${params.patient_id}"
    publishDir "${params.scrna_outdir}/08_final", mode: 'copy'

    input:
    path markers_object
    val  patient_id

    output:
    path "cluster_labels.csv", emit: labels_csv
    path "annotated_umap.png"
    path "final_annotated.rds", emit: final_object

    script:
    """
    Rscript ${moduleDir}/../bin/annotate_and_plot.R ${markers_object} ${patient_id}
    """
}
