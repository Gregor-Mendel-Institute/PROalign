#!/bin/bash
set -o errexit
set -o pipefail

T_DIR="test_01"
mkdir -p ${T_DIR}

#

nextflow run ../PROalign.nf -profile test_local \
  --files "test_files/fastq/pe/*_R{1,2}_*.fastq.gz" \
  --output ${T_DIR}/both \
  --seqtype "PE" \
  --FIVEP_UMI true \
  --THREEP_UMI true \
  -w ${T_DIR}/work -resume



if [ ! $(grep "umi_loc=per_read" ${T_DIR}/both/logs_and_QC/fastp/*.log | wc -l ) -gt 0 ]; then
    echo "wrong umi both"
    exit 1
  fi


#  nextflow run ../PROalign.nf -profile test_local \
#    --files "test_files/fastq/pe/*_R{1,2}_*.fastq" \
#    --output ${T_DIR}/five \
#    --seqtype "PE" \
#    --FIVEP_UMI true \
#    --THREEP_UMI false \
#    -w ${T_DIR}/work

#if [ ! $(grep "umi_loc=read2" ${T_DIR}/five/logs_and_QC/fastp/*.log | wc -l ) -gt 0 ]; then
#    echo "wrong umi five"
#    exit 1
#fi



  #nextflow run ../PROalign.nf -profile test_local \
  #    --files "test_files/fastq/pe/*_R{1,2}_*.fastq" \
  #    --output ${T_DIR}/three \
  #    --seqtype "PE" \
  #    --FIVEP_UMI false \
  #    --THREEP_UMI true \
  #    -w ${T_DIR}/work -resume


  #if [ ! $(grep "umi_loc=read1" ${T_DIR}/three/logs_and_QC/fastp/*.log | wc -l ) -gt 0 ]; then
  #      echo "wrong umi three"
  #      exit 1
  #fi


  #nextflow run ../PROalign.nf -profile test_local \
  #    --files "test_files/fastq/pe/*_R{1,2}_*.fastq" \
#      --output ${T_DIR}/none \
#      --seqtype "PE" \
#      --FIVEP_UMI false \
#      --THREEP_UMI false \
#      -w ${T_DIR}/work

#      if [ $(grep "umi" ${T_DIR}/none/logs_and_QC/fastp/*.log | wc -l ) -gt 0 ]; then
#            echo "wrong no umi"
#            exit 1
#      fi
