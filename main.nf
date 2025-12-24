#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
    Mycobacterium tuberculosis (SNP Distance) Phylogeny & Visualization Pipeline (v${params.version})
    Created by: SPHERES Lab Team
"""

include { PHYLO_ANALYSIS } from './workflows/phylo.nf'
include { VISUALIZATION }  from './workflows/visualization.nf'
include { CONSENSUS } from './workflows/consensus.nf'
include { VERSIONS }       from './workflows/utils.nf'

workflow {
    fhir_ch = Channel.fromPath("${params.fhir_dir}/*.json", checkIfExists: true)
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true).first()

    PHYLO_ANALYSIS(fhir_ch, ref_ch)
    VISUALIZATION(PHYLO_ANALYSIS.out.matrix, PHYLO_ANALYSIS.out.metadata, PHYLO_ANALYSIS.out.tree)
    CONSENSUS(fhir_ch, ref_ch)
    VERSIONS()
}
