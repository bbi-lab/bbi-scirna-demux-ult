# bbi-scirna-demux-ult

## Intro

This *bbi-scirna-demux-ult* pipeline runs a *cram2bam* program to make *.bam* files from Ultima filtered UCRAM files.

## Summary of Ultima data processing.

1. download the raw Ultima CRAM files
2. run the Ultima trimmer on the raw Ultima CRAM files: see the *bbi-ultima* repo
3. make symbolic links to the trimmed CRAM files: see the *bbi-scirna-demux-ult* repo
4. make a sampleheet JSON file: see thhe *bbi-scirna-demux/samplesheet* repo
5. make an experiment.config file: see the *bbi-scirna-demux-ult* repo
6. run the *bbi-scirna-demux-ult* pipeline: see the *bbi-scirna-demux-ult* repo
7. run the *bbi-scirna-analyze* pipeline: see the *bbi-scirna-analyze* repo

## Installation

Install the following software

- Nextflow: this pipeline uses Nextflow DSL2 so you must install a recent version of Nextflow. I use version 24.10.2 successfully. If you need to run the *bbi-dmux* and *bbi-sci* pipelines too, you will need two different Nextflow version so install the new Nextflow in its own location because the recent versions no longer support DSL1.
- Rust: the *cram2bam* program is written in Rust so you must install the Rust compiler.
- cram2bam: this program is a compiled program written in Rust so it must be compiled and the executable copied to the *bbi-scirna-demux-ult/bin* directory. The *bbi-scirna-demux-ult/bin* directory has a *cram2bam* executable that was made on a shendure GS-IT node so you will not need to compile it unless the executable fails to run on your system.
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

## Make symbolic links to the trimmed CRAM files

After running the Ultima trimmer on the raw Ultima CRAM files, you need to make symbolic links to the trimmed CRAM files because this *bbi-scirna-demux-ult* pipeline expects input CRAM files with standard names. The name format is

```
symlink-<lane_number>-<P5_index_id>.trim.cram
```

For example, *symlink-001-Z0140.trim.cram*. The P5 index id starts with 'Z' followed by a 4 digit integer for the index.

This repository includes the bash script *scripts/make_cram_symlinks.sh* for making these symbolic links. You will need to edit the file to set the shell variables

- lane_id [use a 3-digit integer, starting with 001, for the CRAM files from each wafer]
- in_dir [path to the directory that has the trimmed CRAM files]
- lcram [a command string to select the CRAM files from the wafer]

The *lcram* command string must include the wafer id in order to select CRAM files from the wafer. The script stores the symbolic links in the *in_dir* directory.


### Edit the *experiment.config* file.

Edit *experiment.config* to set the following parameters for your run:

- params.samplesheet_json: the path to your samplesheet JSON file for the run
- params.ultima_cram_dir: the path to the directory that contains symbolic links to the Ultima filtered CRAM files
- params.output_dir: the path to the directory where the data are processed. The output BAM files are written to the directory *$params.output_dir/demux_out*

### Make the samplesheet JSON file.

The file *bbi-scirna-demux/samplesheet/scirna_samplesheet.py* has detailed information for making a samplesheet JSON file. The steps are

- run *bbi-scirna-demux/samplesheet/lims2scrunch.py* on a LIMS CSV manifest file to make a CSV samplesheet file where each row describes a sample. Run *lims2scrunch.py --help* command for more information.
- run *bbi-scirna-demux/samplesheet/samplesheet_scrunch.py* on a samplesheet CSV file that is suitable for the *bbi-dmux* pipeline. In the simplest case, the input file has one row per RT well and the output CSV file has one row per sample. *samplesheet_scrunch.py* also adds columns that give the PCR primer wells or columns and rows. Run *samplesheet_scrunch.py --help* command for more information.
- run *bbi-scirna-demux/samplesheet/scirna_samplesheet.py* to convert the scrunched CSV file to a JSON file. *scirna_samplesheet.py* requires a command line parameter that gives the number of lanes used in the sequencing run. Run *scirna_samplesheet.py -d* for detailed documentation. (At this early stage of the program's life, there may be omissions and errors in the documentation.) You may need to edit the scrunched CSV file in a spreadsheet program in order to add columns described in the *scirna_samplesheet.py* documentation.

### Run *bbi-scirna-demux-ult*

I recommend that you use the *run.demux.sh* script in this repository. You must edit the script to use the correct Nextflow program and Nextflow main.nf script.

### Run *bbi-scirna-analyze

For information about running the *bbi-scirna-analyze* pipeline, see *bbi-scirna-analyze/README.md*.
