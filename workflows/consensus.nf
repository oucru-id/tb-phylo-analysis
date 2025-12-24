nextflow.enable.dsl = 2

process BUILD_CONSENSUS {
    publishDir "${params.results_dir}/consensus/", mode: 'copy'

    input:
    path fhir_file
    path reference

    output:
    path "*.fasta", emit: fasta

    script:
    """
    python3 $baseDir/scripts/fhir_consensus.py \\
        --input ${fhir_file} \\
        --reference ${reference} \\
        --output ${fhir_file.baseName.replaceAll('.fhir', '').replaceAll('.json', '')}.consensus.fasta
    """
}

workflow CONSENSUS {
    take:
    fhir_files
    reference

    main:
    BUILD_CONSENSUS(fhir_files, reference)
}