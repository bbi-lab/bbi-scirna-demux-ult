#!/bin/bash

#
# Make symlinks to trimmed Ultima CRAM files.
#
# Notes:
#   o  the trimmed Ultima CRAM file name format varies from
#      run-to-run
#   o  the bbi-scirna-demux-ult pipeline expects that input
#      CRAM files have a standard name
#   o  this script makes symlinks with the standard names
#   o  this script links wafer numbers to the equivalence
#      of lane numbers, which are required by this pipeline.
#   o  the CRAM files from each wafer is assigned a distinct
#      lane number starting with 001 and increasing
#      by one.
#
#   o  trimmed CRAM files have names like
#
#        R097C_GEX_hash_oligo-436012-Z0001.trim.cram
#        ^                    ^      ^
#        experiment           wafer  P5
#        id                   id     index
#
#   o  symlink name format is
#
#        symlink-001-Z0001.trim.cram
#                ^   ^
#              lane  P5
#              id    id
#
#   o  edit and run this script for each wafer in the
#      data set
#   o  set the lane_id, in_dir, and lcram as needed
#   o  edit lcram to select only the CRAM files for one
#      wafer
#

lane_id="001"

in_dir="/net/bbi/vol2/seq/ug100/NVUS2024101701-09.trim"
lcram=`ls ${in_dir}/R097C_GEX_hash_oligo-436012-*.trim.cram`

for cram in $lcram
do
  filename=`basename ${cram}`
  p5_id=`echo ${filename} | awk 'BEGIN{FS="-"}{print $3}'`
  echo "${in_dir}/symlink-${lane_id}-${p5_id}"
  ln -s ${cram} ${in_dir}/symlink-${lane_id}-${p5_id}
  echo
done
