#!/bin/bash

#
# Edit the NEXTFLOW path to your Nextflow installation.
#
NEXTFLOW="/net/gs/vol1/home/bge/src/nextflow/nextflow.24.04.2.5914"

#
# Edit the MAIN_NF path to the main.nf file in your
# bbi-scirna-demux installation.
#
MAIN_NF="/net/gs/vol1/home/bge/git/bbi-scirna-analyze/main.nf"


NOW=`date '+%Y%m%d_%H%M%S'`
WORK_DIR="$PWD/work_analyze"
TRACE_FILE="$PWD/trace.analyze.${NOW}.tsv"
CONFIG_FILE="$PWD/experiment.config"

$NEXTFLOW run $MAIN_NF -c $CONFIG_FILE -w $WORK_DIR -with-trace $TRACE_FILE  -resume

