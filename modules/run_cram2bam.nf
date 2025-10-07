def demux_out = params.output_dir + '/demux_out'


process run_cram2bam {
  cache 'lenient'

  publishDir path: "${demux_out}", pattern: "*.bam", mode: 'copy'

  input:
  tuple val(cram_file_in), val(lane_index), val(pcr7_index), val(pcr5_index)
  val(samplesheet_file_in)
  val(rt_barcode_file_in)
  val(ligation_barcode_file_in)
  val(number_threads)

  output:
  path("*.bam")

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  echo "cram_file: ${cram_file_in}"

  $workflow.projectDir/bin/cram2bam -i ${cram_file_in} \
                                    -s ${samplesheet_file_in} \
                                    -r ${rt_barcode_file_in} \
                                    -l ${ligation_barcode_file_in} \
                                    --lane_index ${lane_index} \
                                    -7 ${pcr7_index} \
                                    -5 ${pcr5_index} \
                                    -t ${number_threads}
  """
}


