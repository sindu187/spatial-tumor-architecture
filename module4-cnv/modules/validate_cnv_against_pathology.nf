process VALIDATE_CNV_AGAINST_PATHOLOGY {
    tag "${params.patient_id}"
    publishDir "${params.cnv_outdir}/04_validate", mode: 'copy'

    input:
    path infercnv_object
    path annotation_file

    output:
    path "cnv_scores_per_spot.csv",     emit: cnv_scores
    path "cnv_validation_summary.csv"
    path "cnv_score_boxplot.png"

    script:
    """
    Rscript ${moduleDir}/../bin/validate_cnv_against_pathology.R \
        ${infercnv_object} ${annotation_file} "${params.ref_groups}" "${params.malignant_groups}"
    """
}
