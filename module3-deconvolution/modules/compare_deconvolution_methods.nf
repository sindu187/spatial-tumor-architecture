process COMPARE_DECONVOLUTION_METHODS {
    tag "${params.patient_id}"
    publishDir "${params.deconvolution_outdir}/03_compare_methods", mode: 'copy'

    input:
    path anchor_predictions
    path marker_scores

    output:
    path "method_agreement.csv"

    script:
    """
    Rscript ${moduleDir}/../bin/compare_deconvolution_methods.R ${anchor_predictions} ${marker_scores}
    """
}
