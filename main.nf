#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.output_dir= "$projectDir/results"

include { FLUA } from "./modules/fluA.nf"

workflow {
  FLUA()
}
