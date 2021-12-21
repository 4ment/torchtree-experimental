#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

ngen = 10000000
iterations=100000
input_seq="$baseDir/datasets/fluA.HA.1981-1998/fluA.HA.1981-1998.fa"
lsd_dates="$baseDir/datasets/fluA.HA.1981-1998/dates-lsd.txt"
results="${params.output_dir}/fluA.HA.1981-1998/"
beast_results="${results}/beast/"
include { RUN_LSD; RUN_FASTTREE } from "./utils.nf"

process CREATE_BEAST {
  publishDir "${beast_results}/${coalescent}/${clock}/${sitemodel}/${categories}/${susbtmodel}/", mode: 'copy'

  input:
    tuple val(susbtmodel), val(categories), val(sitemodel), val(coalescent), val(clock), val(topology)
  output:
    tuple val(susbtmodel), val(categories), val(sitemodel), val(coalescent), val(clock), val(topology), path("fluA.xml")
  """
  generate-beast.py --input $baseDir/xml/fluA.xml.template \
					--output fluA.xml \
				  --stem fluA \
					--ngen ${ngen} \
					--categories ${categories} \
  				--coalescent ${coalescent} \
  				--clock ${clock} \
					--model ${susbtmodel} \
					--site ${sitemodel} \
          --topology ${topology}
  """
}

process RUN_BEAST {
  publishDir "${beast_results}/${coalescent}/${clock}/${sitemodel}/${categories}/${susbtmodel}/", mode: 'copy'

  input:
    tuple val(susbtmodel), val(categories), val(sitemodel), val(coalescent), val(clock), val(topology), path(xml_file)
  output:
    tuple val(susbtmodel), val(categories), val(sitemodel), val(coalescent), val(clock), val(topology)
    path("fluA.trees")
    path("fluA.log")

  """
  beast ${xml_file}
  """
}

process RUN_TORCHTREE_ROOTED {
  publishDir "${results}/${engine}/${vardist}/${coalescent}/${clock}/${sitemodel}/${categories}/${susbtmodel}/${parameterization}/", mode: 'copy'

  input:
    tuple val(susbtmodel), val(categories), val(sitemodel), val(coalescent), val(clock), val(tree_file), val(vardist), val(parameterization), val(engine)
  output:
    path("torchtree.json")
    path("checkpoint.json")
    path("samples.csv")
    path("samples.trees")
    path("torchtree.log")
    path("torchtree.txt")
  script:
    if (sitemodel == "invariant")
      sitemodel_args = " -I"
    else if (categories == 1)
      sitemodel_args = " "
    else
      sitemodel_args = " -C ${categories}"
    if (engine == "bito")
      engine_arg = " --engine bitorch"
    else
      engine_arg = " "

  """
  torchtree-cli advi -i ${input_seq} \
  					 -t ${tree_file} \
  					 -m ${susbtmodel} \
  					 -C ${categories} \
  					 --tol_rel_obj 0 \
  					 --iter ${params.iterations} \
  					 --clock ${clock} \
  					 --coalescent ${coalescent} \
  					 --heights ${parameterization} \
             -q ${vardist} \
             ${sitemodel_args} \
             ${engine_arg} \
  					 > torchtree.json
  { time \
    torchtree torchtree.json  > torchtree.txt ; } 2> torchtree.log
  """
}

workflow FLUA{
  ch_fasta = Channel.from(input_seq)
  ch_outtree = Channel.from("fluA.HA.1981-1998-fasttree.tree")
  ch_out = Channel.from("${params.output_dir}/fluA.HA.1981-1998")

  RUN_FASTTREE(ch_fasta, ch_outtree, ch_out)

  RUN_LSD(RUN_FASTTREE.out,
      "${lsd_dates}",
      Channel.from(986),
      "${params.output_dir}/fluA.HA.1981-1998")

	subst = Channel.from('JC69', 'HKY', 'GTR')
	categories = Channel.from(1, 4)
  sitemodel = Channel.from('weibull', 'gamma', 'invariant')
  w = subst.combine(categories).combine(Channel.from('gamma'))
  g = subst.combine(Channel.from(4)).combine(Channel.from('weibull'))
  i = subst.combine(Channel.from(1)).combine(Channel.from('invariant'))
  coalescent = Channel.from('constant', 'skyride')
  clock = Channel.from('strict', 'ucln')
  base_ch = w.mix(g,i).combine(coalescent).combine(clock)

  CREATE_BEAST(base_ch.combine(RUN_LSD.out[0])) | RUN_BEAST

  parameterization = Channel.from('ratio')//, 'shift')
  variational = Channel.from('meanfield')//, 'fullrank', 'realnvp')
  engines = Channel.from("", "bito")
  torchtree_ch = base_ch.combine(RUN_LSD.out[0]).combine(variational).combine(parameterization).combine(engines)
  RUN_TORCHTREE_ROOTED(torchtree_ch)
}
