
process RUN_LSD {
  publishDir "${output}", mode: 'copy'

  input:
    path(input_tree)
    path(dates_file)
    val(length)
    val(output)
  output:
    path("lsd.out.date.nexus")
    path "lsd.out.nwk", emit: lsd_tree_newick // branch=subst
    path "lsd.out"

  """
  lsd2 -i ${input_tree} \
       -d ${dates_file} \
       -o lsd.out \
       -s ${length} \
       -l -1 \
       -r a
  """
}

process RUN_FASTTREE {
  publishDir "${output}", mode: 'copy'

  input:
    path(input_file)
    val(output_tree)
    val(output)
  output:
    path(output_tree)

  """
  FastTree -nt ${input_file} > ${output_tree}
  """
}