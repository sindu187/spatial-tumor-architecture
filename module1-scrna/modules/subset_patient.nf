process SUBSET_PATIENT {
    tag "${params.patient_id}"
    publishDir "${params.scrna_outdir}/01_subset", mode: 'copy'

    input:
    path matrix_mtx
    path genes_tsv
    path barcodes_tsv
    path metadata_csv
    val  patient_id

    output:
    path "raw_unfiltered.rds", emit: raw_object

    script:
    """
    Rscript ${moduleDir}/../bin/subset_patient.R \
        ${matrix_mtx} ${genes_tsv} ${barcodes_tsv} ${metadata_csv} ${patient_id}
    """
}
