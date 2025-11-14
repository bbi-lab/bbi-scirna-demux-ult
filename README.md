# bbi-scirna-demux-ult

## Intro

This *bbi-scirna-demux-ult* pipeline runs a *cram2bam* program to make *.bam* files from Ultima filtered ucram files.

## Installation

Install the following software

- Nextflow: this pipeline uses Nextflow DSL2 so you must install a recent version of Nextflow. I use version 24.10.2 successfully. If you need to run the *bbi-dmux* and *bbi-sci* pipelines too, you will need two different Nextflow version so install the new Nextflow in its own location because the recent versions no longer support DSL1.
- Rust: the *cram2bam* program is written in Rust so you must install the Rust compiler.
- cram2bam: this program is a compiled program written in Rust so you must compile it and copy the executable to the *bbi-scirna-demux-ult/bin* directory.
- python3 interpreter: I use version 3.12.1 successfully.

### Install Nextflow

See the Nextflow installation instructions at

https://www.nextflow.io/docs/latest/install.html

### Install Rust

See the Rust installation instructions at

https://www.rust-lang.org/tools/install

### Build and install *cram2bam*

Run the following commands

```
cd bbi-scirna-demux-ult/src/cram2bam
cargo build --release
cp target/release/cram2bam ../../bin
```

I recommend that you build *cram2bam* on a newer cluster node, for example, s020 on the Shendure cluster.

## Run *bbi-scirna-demux-ult*

I recommend that you use the *run.demux.sh* script in this repository. You must edit the script to use the correct Nextflow program and Nextflow main.nf script.

### Edit the *experiment.config* file.

Edit *experiment.config* to set the following parameters for your run:

```
params.samplesheet_json
params.ultima_cram_dir
params.output_dir
```

*params.samplesheet_json* is the path to your samplesheet JSON file for the run. *params.ultima_cram_dir* is the path to the directory that contains the Ultima filtered CRAM files, and *params.output_dir* is the path to the directory where the processing output is written to.

### Make the samplesheet JSON file.

This is a bit lengthy at this time. The steps are

- run *bbi-scirna-demux/samplesheet/lims2scrunch.py* on a LIMS CSV manifest file to make a CSV samplesheet file where each row describes a sample. Run *lims2scrunch.py --help* command for more information.
- run *bbi-scirna-demux/samplesheet/samplesheet_scrunch.py* on a samplesheet CSV file that is suitable for the *bbi-dmux* pipeline. In the simplest case, the input file has one row per RT well and the output CSV file has one row per sample. *samplesheet_scrunch.py* also adds columns that give the PCR primer wells or columns and rows. Run *samplesheet_scrunch.py --help* command for more information.
- run *bbi-scirna-demux/samplesheet/scirna_samplesheet.py* to convert the scrunched CSV file to a JSON file. *scirna_samplesheet.py* requires a command line parameter that gives the number of lanes used in the sequencing run. Run *scirna_samplesheet.py -d* for detailed documentation. (At this early stage of the program's life, there may be omissions and errors in the documentation.) You may need to edit the scrunched CSV file in a spreadsheet program in order to add columns described in the *scirna_samplesheet.py* documentation.

###

Run the Ultima trimmer on the raw Ultima cram files. The trimmer program runs from a Docker container. See the Ultima trimmer documentation. There is additional documentation in our git repository bbi-ultima/ultima-trimmer. There is a bash script that runs the trimmer in the git repository

  bbi-ultima/ultima-trimmer/run_trimmer.sh

Use the script

  bbi-scirna-demux-ult/scripts/make_cram_symlinks.sh

to make symbolic links to the trimmer cram files. This bbi-scirna-demux-ult pipeline looks for cram files that have the symlink filenames.

### Run *bbi-scirna-demux-ult*

Use the *run.demux.sh* bash script to start the pipeline run.

The output *.bam* files are in the directory demux_out. The *.bam* file contains unaligned read sequences and quality values as well as barcode data in the SAMtags.

