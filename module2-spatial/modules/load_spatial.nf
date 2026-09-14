process LOAD_SPATIAL {
    tag "${params.patient_id}"
    publishDir "${params.spatial_outdir}/01_load", mode: 'copy'

    input:
    path spaceranger_dir
    path metadata_csv
    val  patient_id

    output:
    path "raw_spatial.rds", emit: raw_object

    script:
    """
    Rscript ${moduleDir}/../bin/load_spatial.R ${spaceranger_dir} ${metadata_csv} ${patient_id}
    """
}
