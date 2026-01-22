#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
    Mycobacterium tuberculosis Federated Phylogeny & Visualization Pipeline (Local Lab)
    Developed by SPHERES Lab Team
"""

include { PHYLO_ANALYSIS } from './workflows/phylo.nf'
include { VISUALIZATION }  from './workflows/visualization.nf'
include { VERSIONS }       from './workflows/utils.nf'

/* process FETCH_FROM_FHIR {
    publishDir "${params.results_dir}/fetched_data", mode: 'copy'

    input:
    val url
    val auth
    val since

    output:
    path "*.json", emit: json_files

    script:
    def date_arg = since ? "--since ${since}" : ""
    """
    python3 $baseDir/scripts/fetch_fhir_data.py \\
        --url "${url}" \\
        --auth "${auth}" \\
        $date_arg
    """
} */

workflow {
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
    anchor_ch = Channel.fromPath("$baseDir/data/anchor/*.json").collect().ifEmpty([])

    // if (params.fhir_server_url && params.fhir_server_url != "null") {
    //     log.info "Using FHIR Server: ${params.fhir_server_url}"
        
    //     FETCH_FROM_FHIR(params.fhir_server_url, params.fhir_server_auth, params.fetch_since)
    //     fhir_ch = FETCH_FROM_FHIR.out.json_files.flatten()
    // } else {        
        fhir_ch = Channel.fromPath("${params.fhir_dir}/*.json", checkIfExists: true)
    // }

    PHYLO_ANALYSIS(fhir_ch, ref_ch, anchor_ch)
    VISUALIZATION(PHYLO_ANALYSIS.out.matrix, PHYLO_ANALYSIS.out.metadata, PHYLO_ANALYSIS.out.tree)
    VERSIONS()
}
