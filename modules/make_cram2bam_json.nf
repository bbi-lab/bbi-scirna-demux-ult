process make_cram2bam_json {
  input:
  path(samplesheet_file)
  val(cram_dir)

  output:
  path("cram2bam.json")

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  $workflow.projectDir/bin/make_cram2bam_json.py -i ${samplesheet_file} -d ${cram_dir}
  """
}

