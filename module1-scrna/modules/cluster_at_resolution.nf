process CLUSTER_AT_RESOLUTION {
    tag "resolution=${resolution}"
    publishDir "${params.scrna_outdir}/04_resolution_sweep", mode: 'copy'

    input:
    path base_object
    val  resolution

    output:
    path "clustered_res*.rds", emit: clustered_object
    path "ari_res*.csv",       emit: ari_result

    script:
    """
    Rscript ${moduleDir}/../bin/cluster_at_resolution.R ${base_object} ${resolution}
    """
}
