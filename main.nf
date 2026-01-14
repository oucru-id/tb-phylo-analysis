#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
    Mycobacterium tuberculosis Federated Phylogeny & Visualization Pipeline
    Developed by SPHERES Lab Team
"""

include { PHYLO_ANALYSIS } from './workflows/phylo.nf'
include { VISUALIZATION }  from './workflows/visualization.nf'
include { FEDERATED_NEXTSTRAIN } from './workflows/federated.nf'
include { VERSIONS }       from './workflows/utils.nf'

process FETCH_FROM_FHIR {
    publishDir "${params.results_dir}/fetched_data", mode: 'copy'

    input:
    val url
    val auth
    val since
    val token_url
    val client_id
    val client_secret
    val scope

    output:
    path "*.json", emit: json_files

    script:
    def date_arg = since ? "--since ${since}" : ""
    """
    python3 $baseDir/scripts/fetch_fhir_data.py \\
        --url "${url}" \\
        --auth "${auth}" \\
        --token-url "${token_url}" \\
        --client-id "${client_id}" \\
        --client-secret "${client_secret}" \\
        --scope "${scope}" \\
        $date_arg
    """
}

workflow {
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true).first()

    if (params.fhir_server_url && params.fhir_server_url != "null") {
        log.info "Using FHIR Server: ${params.fhir_server_url}"
        
        FETCH_FROM_FHIR(
            params.fhir_server_url,
            params.fhir_server_auth,
            params.fetch_since,
            params.fhir_oauth_token_url,
            params.fhir_client_id,
            params.fhir_client_secret,
            params.fhir_scope
        )
        fhir_ch = FETCH_FROM_FHIR.out.json_files.flatten()
    } else {
        log.info "Using Local Directory: ${params.fhir_dir}"
        
        fhir_ch = Channel.fromPath("${params.fhir_dir}/*.json", checkIfExists: true)
    }

    PHYLO_ANALYSIS(fhir_ch, ref_ch)
    VISUALIZATION(PHYLO_ANALYSIS.out.matrix, PHYLO_ANALYSIS.out.metadata, PHYLO_ANALYSIS.out.tree)
    FEDERATED_NEXTSTRAIN(PHYLO_ANALYSIS.out.metadata, PHYLO_ANALYSIS.out.fasta, ref_ch)
    VERSIONS()
}
