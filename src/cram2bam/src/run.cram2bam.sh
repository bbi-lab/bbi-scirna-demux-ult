#!/bin/bash

cram_in="/net/bbi/vol2/seq/ug100/NVUS2024101701-02.trim/symlink-001-Z0097.trim.cram"
samplesheet="/net/bbi/vol2/data/sciRNAseq/ug100/REF3/samplesheet.ref3.json"
rt_barcodes="/net/gs/vol1/home/bge/git/bbi-scirna-demux/data/rt.txt"
lig_barcodes="/net/gs/vol1/home/bge/git/bbi-scirna-demux/data/ligationwellbarcodeMEGA.sorted_by_row.tsv"

../target/release/cram2bam -i ${cram_in} -s ${samplesheet} -r ${rt_barcodes} -l ${lig_barcodes} --lane_index 1 -7 0 -5 97 -t 1
