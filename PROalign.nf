#!/usr/bin/env nextflow

/*********************
* parameters section *
**********************/
params.output="$baseDir/results/"

params.THREADS=50 // Threads to use for multithreaded applications
params.UMI_LEN=6  // Length of UMI in basepairs

// UMI Flags (set to Y or N as appropriate)
params.FIVEP_UMI=true // Is there a UMI on the 5' end of the read?
params.THREEP_UMI=true // Is there a UMI on the 3' end of the read?

if (params.FIVEP_UMI & params.THREEP_UMI){
  params.umi=true
  params.umi_loc="per_read"
  }
else if (params.FIVEP_UMI){
  params.umi=true
  params.umi_loc="read2"
  }
else if (params.THREEP_UMI){
  params.umi=true
  params.umi_loc="read1"
}
else {
  params.umi=false
  params.umi_loc=false
}

// Adaptor sequences to clip. Default = Tru-Seq small RNA
params.ADAPTOR_1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
params.ADAPTOR_2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"



// Genomes. Fill in paths.
//GENOME_EXP="/home/jaj256/genome/dm6/dm6Hsp70AaOnly"
//GENOME_SPIKE="/home/jaj256/genome/dm6hg38/dm6hg38" ## USE REPEAT MASKED VERSION!!
//SPIKE_PREFIX="hg38" ## This is the prefix you've used on your spike in chromosomes, ie >hg38chr1
//RDNA="/home/jaj256/genome/dm3hg38/dm3hg38rDNA"

// Mapq value for filtering multimappers
params.MAPQ=10

// channel

files = Channel.fromFilePairs( params.files, size: -1 ).ifEmpty { error "Invalid path: point to path containing fq files using --files flag" }



// unzip

// rename?

//fastqc

//"detecting paired end files"?


//Trimming adapters and filtering rRNA reads

process trim_fastp  {
publishDir "$params.output/logs_and_QC/fastp", mode: 'copy', pattern: '*.html'
publishDir "$params.output/logs_and_QC/fastp", mode: 'copy', pattern: '*.log'

input:
set id, file(reads) from files

output:
set id, file("${id}_R*.fq.gz") into trimmed
file("${id}_fastp.html")
file("${id}_fastp.log")

script:
if (params.umi)
"""
fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${id}_R1.fq.gz \
    -O ${id}_R2.fq.gz \
    --adapter_sequence ${params.ADAPTOR_1} \
    --adapter_sequence_r2 ${params.ADAPTOR_2} \
    --umi \
    --umi_loc=${params.umi_loc} \
    --umi_len=${params.UMI_LEN} \
    --html ${id}_fastp.html \
    -c \
    --overlap_len_require 15 2> ${id}_fastp.log
"""
else
"""
fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o ${id}_R1.fq.gz \
    -O ${id}_R2.fq.gz \
    --adapter_sequence ${params.ADAPTOR_1} \
    --adapter_sequence_r2 ${params.ADAPTOR_2} \
    --html ${id}_fastp.html \
    -c \
    --overlap_len_require 15 2> ${id}_fastp.log
"""
}

//rRNA filter
/*
process rRNA {

  input:
  set id, file(reads) from trimmed

  output:

script:
"""
bowtie2 \
--fast-local \
--un-conc trimmedFastq/${id}.fastq \
--interleaved - \
-x ${RDNA} 2> logs/rRNA/${id}_rRNA_bowtie.log) > /dev/null



"""

}
*/

process align_w_spikein {

 input:
 set id, file(reads) from trimmed

 output:

 script:
 """
 bowtie2 \
    --local \
    --very-sensitive-local \
    --no-unal \
    --no-mixed \
    --no-discordant \
    -x "$GENOME_SPIKE" \
    -1 ${reads[0]} \
    -2  ${reads[1]} \
    2> ${id}_spikeAlign.log) |
    samtools view -hS -f 2 -q ${params.MAPQ} |
    perl -n -e 'print $_ if (/^\@/ || /'${params.SPIKE_PREFIX}'/ ) ' |
    samtools view -b |
    samtools sort -@ -o ${id}.BAM
    samtools index ${id}.BAM
 """

}
//Cleaning up filenames in trimmedFastq (bowtie automatically names PE --un output)


// Aligning to spike in genome to get normalization factors

// Aligning to experimental genome


// deduplicating with UMIs


// Generating infoTable
