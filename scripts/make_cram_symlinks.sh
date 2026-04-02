#!/bin/bash

#
# Make symlinks to trimmed Ultima CRAM files.
#
# Notes:
#   o  the trimmed Ultima CRAM file name format may varies
#      run-to-run
#   o  the bbi-scirna-demux-ult pipeline expects that input 
#      CRAM files have a standard name
#   o  this script makes symlinks with the standard names
#   o  this script links wafer numbers to the equivalence
#      of lane numbers, which are required by this pipeline.
#   o  the CRAM files from each wafer is assigned a distinct
#      lane number starting with 001 and increasing
#      by one.
#   o  note that each wafer must have distinct p5 indices;
#      that is, no two wafers can have the same p5 indices.
#   o  trimmed CRAM files have names like
# 
#        441449-R100C_GEX_hash_oligo-Z0001-CAGCTCGAATGCGAT.trim.cram
#        ^      ^                    ^
#        wafer  experiment           P5
#        id     id                   index
#
#   o  symlink name format and directory
#       
#        R097/symlink-001-Z0001.trim.cram
#         ^             ^   ^
#        experiment   lane  P5
#        id           id    id
#
#   o  edit and run this script for each wafer in the
#      data set
#   o  set the lane_id, in_dir, lcram, and symlink lines
#      as needed
#   o  edit lcram to select only the CRAM files for one
#      wafer
#

lane_id="001"


in_dir="/net/our_lab/vol1/bbi_raw_data/seq/nobackup/NVUS2024101701-32.trim/R100/raw/441389"
lcram=`ls ${in_dir}/441389-R100E_GEX_hash_oligo-*.trim.cram`

for cram in $lcram
do
  filename=`basename ${cram}`
  p5_id=`echo ${filename} | awk 'BEGIN{FS="-"}{print $3}'`
  ln -s ${cram} ${in_dir}/R100/symlink-${lane_id}-${p5_id}.trim.cram
done

