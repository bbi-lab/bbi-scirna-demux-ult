#!/bin/bash

lane_id="001"

in_dir="/net/bbi/vol2/seq/ug100/NVUS2024101701-02.trim"
lcram=`ls ${in_dir}/*.trim.cram`

for cram in $lcram
do
  p5_id=`echo ${cram} | awk 'BEGIN{FS="-"}{print $4}'`
  echo "${in_dir}/symlink-${lane_id}-${p5_id}"
  ln -s ${cram} ${in_dir}/symlink-${lane_id}-${p5_id}
  echo
done
