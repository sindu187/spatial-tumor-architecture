process INVESTIGATE_DISAGREEMENT_SPATIAL {
    tag "${params.patient_id}"
    publishDir "${params.spatial_outdir}/06_investigate", mode: 'copy'

    input:
    path best_object

    output:
    path "cluster_purity_spatial.csv"
    path "full_breakdown_spatial.csv"
    path "investigated_spatial.rds", emit: investigated_object

    script:
    """
    Rscript ${moduleDir}/../bin/investigate_disagreement_spatial.R ${best_object}
    """
}
