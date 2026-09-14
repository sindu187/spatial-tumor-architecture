process ANCHOR_LABEL_TRANSFER {
    tag "${params.patient_id}"
    publishDir "${params.deconvolution_outdir}/01_anchor_transfer", mode: 'copy'

    input:
    path scrna_object
    path spatial_object

    output:
    path "anchor_predictions.csv", emit: predictions_csv
    path "spatial_with_anchor.rds", emit: anchor_object

    script:
    """
    Rscript ${moduleDir}/../bin/anchor_label_transfer.R ${scrna_object} ${spatial_object}
    """
}
