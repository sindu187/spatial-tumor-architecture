process INVESTIGATE_DISAGREEMENT {
    tag "${params.patient_id}"
    publishDir "${params.scrna_outdir}/06_investigate", mode: 'copy'

    input:
    path best_object

    output:
    path "cluster_purity.csv"
    path "minor_type_breakdown.csv"
    path "investigated.rds", emit: investigated_object

    script:
    """
    Rscript ${moduleDir}/../bin/investigate_disagreement.R ${best_object}
    """
}
