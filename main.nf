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
params.p7_barcode_file_default = "$workflow.projectDir/data/p7.txt"
params.p5_barcode_file_default = "$workflow.projectDir/data/p5.txt"
params.rt_barcode_file_default = "$workflow.projectDir/data/rt.txt"
params.ligation_barcode_file_default = "$workflow.projectDir/data/ligation.txt"


/*
** Set up channels.
*/
samplesheet_file = channel.value(params.samplesheet_json)
ultima_cram_dir = channel.value(params.ultima_cram_dir)
p7_barcode_file_default = channel.value(params.p7_barcode_file_default)
p5_barcode_file_default = channel.value(params.p5_barcode_file_default)
rt_barcode_file_default = channel.value(params.rt_barcode_file_default)
ligation_barcode_file_default = channel.value(params.ligation_barcode_file_default)


/*
** Import modules after defining params.* so that
** the parameters are accessible in the modules.
*/
include {run_check_samplesheet} from './modules/run_check_samplesheet.nf'
include { run_rna_rtlig_demux } from './modules/run_rna_rtlig_demux.nf'


/*
** Run pipeline.
*/
workflow {
  run_check_samplesheet(samplesheet_file)
  run_bclconvert.out.flatMap{ make_pairwise_fastq_bclconvert(it) }.set{fastq_pairs}
  run_cram2bam(fastq_pairs, samplesheet_file, rt_barcode_file_default, ligation_barcode_file_default)
}

