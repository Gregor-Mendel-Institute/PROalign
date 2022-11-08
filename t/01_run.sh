#!/bin/bash
set -o errexit
set -o pipefail

T_DIR="test_01"

if [ -d "${T_DIR}/both" ]; then rm -Rf ${T_DIR}/both ; fi
mkdir -p ${T_DIR}


nextflow run ../PROalign.nf -profile test_local \
  --files "test_files/fastq/pe/*_R{1,2}_*.fastq.gz" \
  --output ${T_DIR}/both \
  --FIVEP_UMI true \
  --THREEP_UMI true \
  -w ${T_DIR}/work -resume



if [ ! $(grep "umi_loc=per_read" ${T_DIR}/both/logs_and_QC/fastp/*.log | wc -l ) -gt 0 ]; then
    echo "wrong umi both"
    exit 1
  fi


  #### check output ####

  ## right folders
  if [ ! -d ${T_DIR}/both ]; then
      echo "no results"
      exit 1
  fi

  if [ ! $(ls -1q ${T_DIR}/both/ | wc -l) -eq 3 ] ; then # counts,
        echo "wrong number output dirs"
      #exit 1
  fi

  if [ ! -d ${T_DIR}/both/logs_and_QC ]; then
      echo "no logs"
      exit 1
  fi

  if [ ! -d ${T_DIR}/both/bams ]; then
      echo "no bams"
      exit 1
  fi

  if [ ! -d ${T_DIR}/both/bw ]; then
      echo "no bw"
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
