import java.nio.file.Paths;

/*
** Run parameters.
*/


/*
** 
*/
params.bin_dir = workflow.projectDir + '/bin'


/*
** Default barcode file paths.
*/
params.rt_barcode_file_default = "$workflow.projectDir/data/rt.txt"
params.ligation_barcode_file_default = "$workflow.projectDir/data/ligation_megasci.row_sorted.tsv"
params.number_threads_cram2bam = 4

/*
** Set up channels.
*/
samplesheet_file = channel.value(params.samplesheet_json)
ultima_cram_dir = channel.value(params.ultima_cram_dir)
rt_barcode_file_default = channel.value(params.rt_barcode_file_default)
ligation_barcode_file_default = channel.value(params.ligation_barcode_file_default)
number_threads_cram2bam = channel.value(params.number_threads_cram2bam)

/*
** Import modules after defining params.* so that
** the parameters are accessible in the modules.
*/
include {run_check_samplesheet} from './modules/run_check_samplesheet.nf'
include {make_cram2bam_json} from './modules/make_cram2bam_json.nf'
include {run_cram2bam} from './modules/run_cram2bam.nf'


def run_cram2bam_closure = {
  item ->
    def in_file = item['in_file']
    def lane_index = item['lane_index']
    def pcr7_index = item['pcr7_index']
    def pcr5_index = item['pcr5_index']
    [in_file, lane_index, pcr7_index, pcr5_index]
}


/*
** Run pipeline.
*/
workflow {
  run_check_samplesheet(samplesheet_file)
  make_cram2bam_json(samplesheet_file, ultima_cram_dir)
  make_cram2bam_json.out.splitJson().map{run_cram2bam_closure(it)}.set{run_cram2bam_in}
 
run_cram2bam_in.view()

  run_cram2bam(run_cram2bam_in, samplesheet_file, rt_barcode_file_default, ligation_barcode_file_default, number_threads_cram2bam)
}

