#!/usr/bin/env nextflow

nextflow.enable.dsl=2

forward_ch= Channel.of(params.forward_primer)
reverse_ch= Channel.of(params.reverse_primer)


process convertPrimer{
    publishDir params.outdir, mode:'copy'

    input:
    val forward_ch
    val reverse_ch
    
    output:
    path '**.fasta'

    shell:
    """
    Primer.sh ${params.forward_primer} ${params.reverse_primer}
    """
    
    
}


//Process Cutprimers

workflow{
    primer_ch= convertPrimer(forward_ch,reverse_ch) | collect | flatten | view
}