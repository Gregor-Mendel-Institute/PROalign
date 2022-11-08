#!/usr/bin/env nextflow

/*********************
* parameters section *
**********************/
params.output="$baseDir/results/"

params.UMI_LEN=6  // Length of UMI in basepairs
// Mapq value for filtering multimappers
params.MAPQ=10
//params.FIVEP_UMI=true // Is there a UMI on the 5' end of the read?
//params.THREEP_UMI=true // Is there a UMI on the 3' end of the read?

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
  exit 1, "Pipeline only makes sense with UMI"
}

// Adaptor sequences to clip. Default = Tru-Seq small RNA
params.ADAPTOR_1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
params.ADAPTOR_2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

params.indexName = "spombe"
params.index = "library://elin.axelsson/index/index_bowtie2_spombe:v2.4.1-release-55"
params.spike_name = "spombe_cerevisiae"
params.index_spike = "library://elin.axelsson/newindex/index_bowtie2_spombe_scerevisiae:v2.4.1-release55"
params.SPIKE_PREFIX = "R64" //This is the prefix you've used on your spike in chromosomes, ie >hg38chr1
params.rdnaName = "spombe_rdna"
params.index_rdna = "library://elin.axelsson/newindex/index_bowtie2_spombe_rdna:v2.4.1_zach"


// channel

files = Channel.fromFilePairs( params.files, size: -1 ).ifEmpty { error "Invalid path: point to path containing fq files using --files flag (Paired-end only)" }


/********************************
*********************************
* START
*********************************
********************************/

/********************************
* Get indices
********************************/

// Spikein
process get_spike_index{

  input:
  params.index_spike

  output:
  file params.spike_name into comb

  script:
  """
  singularity run ${params.index_spike}
  """
}

// experimental
process get_index{

  input:
  params.index

  output:
  file params.indexName into experiment

  script:
  """
  singularity run ${params.index}
  """
}
// rDNA
process get_rdna_index{

  input:
  params.index_rdna

  output:
  file params.rdnaName into rdna

  script:
  """
  singularity run ${params.index_rdna}
  """
}



/********************************
* fastqc
********************************/
files.into{files;qc_files}

process fastqc {
  publishDir "$params.output/logs_and_QC/fastqc", mode: 'copy'

  input:
  set id, file(fq) from qc_files

  output:
  file "*fastqc.*"

  script:
  """
  fastqc $fq
  """

}

/********************************
* preprocess reads
********************************/

//Trimming adapters and umi

process trim_fastp  {
  publishDir "$params.output/logs_and_QC/fastp", mode: 'copy', pattern: '*.html'
  publishDir "$params.output/logs_and_QC/fastp", mode: 'copy', pattern: '*.log'

input:
  set id, file(reads) from files

output:
  set id, file("${id}_R*.fq.gz") into trimmed
  file("${id}_fastp.html")
  file("${id}_fastp.log") into fastp_log

script:
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
}

//rRNA filter

process rRNA {
publishDir "$params.output/logs_and_QC/rRNA", mode: 'copy', pattern: '*.log'
  input:
  set id, file(reads) from trimmed
  file(index) from rdna

  output:
  set id, file("${id}_1.fastq"), file("${id}_2.fastq") into filtered
  file("${id}_rRNA_bowtie.log") into rRNA_bowtie_log

  script:
  """
  bowtie2 \
  -1 ${reads[0]} \
  -2  ${reads[1]} \
  --fast-local \
  --un-conc ${id}.fastq \
  --interleaved - \
  -x ${index}/${index} 2> ${id}_rRNA_bowtie.log
  mv ${id}.1.fastq ${id}_1.fastq
  mv ${id}.2.fastq ${id}_2.fastq
  """

}

/********************************
* alignments
********************************/
filtered.into{filtered_sp;filtered}

// spike in
process align_w_spikein {
publishDir "$params.output/logs_and_QC/align", mode: 'copy', pattern: '*.log'

  input:
  set id, file(read1), file(read2) from filtered_sp
  file(index) from comb

  output:
  set id, file("spike_${id}.sam") into spike
  file("${id}_spikeAlign.log") into spike_log

  script:
  """
  bowtie2 \
      --local \
      --very-sensitive-local \
      --no-unal \
      --no-mixed \
      --no-discordant \
      -x $index/$index \
      -1 ${read1} \
      -2  ${read2} \
      -S spike_${id}.sam 2> ${id}_spikeAlign.log
 """

}

process get_spike_bam {

  input:
  set id, file(reads) from spike
  val params.SPIKE_PREFIX

  output:
  set id, file("${id}_spike.bam"),file("${id}_spike.bam.bai") into spike_bam

  shell:
  '''
  samtools view -hS -f 2 -q !{params.MAPQ} !{reads} |
  perl -n -e 'print $_ if (/^\\@/ || /!{params.SPIKE_PREFIX}/ )' |
  samtools view -b |
  samtools sort -o !{id}_spike.bam
  samtools index !{id}_spike.bam
  '''

}

// Aligning to experimental genome
process align {
publishDir "$params.output/logs_and_QC/align", mode: 'copy', pattern: '*.log'

  input:
  set id, file(read1), file(read2) from filtered
  file(index) from experiment

  output:
  set id, file("${id}.sam") into align
  file("${id}_align.log") into align_log

  script:
  """
  bowtie2 \
  --local \
  --sensitive-local \
  -x $index/$index \
  -1 ${read1} \
  -2 ${read2} \
  -S ${id}.sam 2> ${id}_align.log
  """
}

process filter_align {

  input:
  set id, file(reads) from align

  output:
  set id, file("${id}.bam"),file("${id}.bam.bai") into align_bam

  script:
  """
  samtools view -bS -f 2 -q ${params.MAPQ} ${reads}|
  samtools sort -o ${id}.bam
  samtools index ${id}.bam
  """
}

/********************************
* dedup based on umi
********************************/

allbam=align_bam.combine(Channel.of('main')).concat(spike_bam.combine(Channel.of('spikein')))

process dedup_umi{
  publishDir "$params.output/logs_and_QC/dedup", mode: 'copy', pattern: '*.log'
  publishDir "$params.output/bams/dedup", mode: 'copy', pattern: '*main_deDuped.bam'

  input:
  set id, file(reads), file(index), type from allbam

  output:
  set id, type, file("*_deDuped.bam") into dedup
  file("*.log") into de_log

  script:
  """
  umi_tools dedup \
  -I ${reads} \
  --umi-separator=":" \
  --paired \
  -S ${id}_${type}_deDuped.bam \
  > ${id}_${type}_deDup.log
  """
}
/********************************
* index
********************************/

process index {

  input:
  set id, type, file(reads) from dedup

  output:
  set id, type, file("${reads}"), file("${reads}.bai") into dedup_i

  script:
  """
  samtools index ${reads}
  """
}


/********************************
* Unnormalized bw files
********************************/
dedup_i.into{f_dedup;r_dedup}

process bw_fwd {
publishDir "$params.output/bw/unnorm", mode: 'copy', pattern: '*.bw'

  input:
  set id, type, file(reads), file(bai) from f_dedup

  output:
  set id, file("*_fwd.bw")

  when:
  type=="main"

  script:
  """
  bamCoverage \
  --bam ${reads} \
  --skipNonCoveredRegions \
  --outFileName ${id}_fwd.bw \
  --binSize 1 \
  --normalizeUsing None \
  --Offset 1 \
  --samFlagInclude 82
  """
}

process bw_rev {
publishDir "$params.output/bw/unnorm", mode: 'copy', pattern: '*.bw'

  input:
  set id, type, file(reads), file(bai) from r_dedup

  output:
  set id, file("*_rev.bw")

  when:
  type=="main"

  script:
  """
  bamCoverage \
  --bam ${reads} \
  --skipNonCoveredRegions \
  --outFileName ${id}_rev.bw \
  --binSize 1 \
  --normalizeUsing None \
  --Offset 1 \
  --samFlagInclude 98
  """
}

/********************************
* Summary table
********************************/

process table {
  publishDir "$params.output/logs_and_QC/info", mode: 'copy', pattern: 'infoTable.tsv'

  input:
  file de from de_log.collect()
  file al from align_log.collect()
  file sp from spike_log.collect()
  file fa from fastp_log.collect()
  file rr from rRNA_bowtie_log.collect()

  output:
  file "infoTable.tsv"

  shell:
  '''
  touch infoTable.tsv
  echo -e Name'\t'RawReads'\t'NonDimerReads'\t'%dimer'\t'insertSize'\t'rRNAreads'\t'%rRNA'\t'passedFilters'\t'\
      bowtieConcordant'\t'bowtieMulti'\t'bowtieUnal'\t'bowtieOverallMap%'\t'bowtieConcordant%'\t'\
      bowtieMulti%'\t'bowtieUnal%'\t'uniqueMapped'\t'uniqueMappedNondup'\t'%PCRdups'\t'uniqueMappedSpikein'\t'\
      uniqueMappedSpikeinNondup'\t'spikeInPCRdups% >> infoTable.tsv
  for SAMPLE in $(ls *_align.log | sed 's/_align.log//' )
  do
      echo ${SAMPLE}
      NAME=${SAMPLE}
      RAW_READS=$(cat ${SAMPLE}_fastp.log |
      grep "total reads:" | head -n 1 |
      awk '{print $3}')
      TRIMMED_READS=$(cat ${SAMPLE}_fastp.log |
      grep "total reads:" | tail -n 1 |
      awk '{print $3}')
      PER_DIMER=$(awk -v a="${TRIMMED_READS}" -v b="${RAW_READS}" 'BEGIN { print  (1-a/b)*100 }')%
      INSERT_SIZE=$(cat ${SAMPLE}_fastp.log |
      grep "Insert size peak" |
      awk '{print $8}')
      PASSED_FILTERS=$(cat ${SAMPLE}_align.log |
      grep "reads; of these:$" |
      awk '{print $1}')
      RRNA=$(awk -v a="${TRIMMED_READS}" -v b="${PASSED_FILTERS}" 'BEGIN { print a-b }')
      PER_RRNA=$(awk -v a="${RRNA}" -v b="${RAW_READS}" 'BEGIN { print a/b*100 }')%
      B_CONC=$(cat ${SAMPLE}_align.log |
      grep "aligned concordantly exactly 1 time$" |
      awk '{print $1}')
      B_MULTI=$(cat ${SAMPLE}_align.log |
      grep "aligned concordantly >1 times$" |
      awk '{print $1}')
      B_UNAL=$(cat ${SAMPLE}_align.log |
      grep "aligned concordantly 0 times$" |
      awk '{print $1}')
      B_OAP=$(cat ${SAMPLE}_align.log |
      grep "overall alignment rate$" |
      awk '{print $1}')
      B_CONC_PER=$(awk -v a="${B_CONC}" -v b="${PASSED_FILTERS}" 'BEGIN { print a/b*100}')%
      B_MULTI_PER=$(awk -v a="${B_MULTI}" -v b="${PASSED_FILTERS}" 'BEGIN { print a/b*100 }')%
      B_UNAL_PER=$(awk -v a="${B_UNAL}" -v b="${PASSED_FILTERS}" 'BEGIN { print a/b*100 }')%
      UNIQ_MAPPED=$(cat ${SAMPLE}_main_deDup.log |
      grep "Input Reads:" | awk '{print $10}')
      UNIQ_MAPPED_DEDUP=$(cat ${SAMPLE}_main_deDup.log |
      grep "Number of reads out:" | awk '{print $8}')
      PER_DUPS=$(awk -v a="${UNIQ_MAPPED_DEDUP}" -v b="${UNIQ_MAPPED}" 'BEGIN { print (1-a/b)*100 }')%
      UNIQ_MAPPED_SPIKE=$(cat ${SAMPLE}_spikein_deDup.log |
      grep "Input Reads:" | awk '{print $10}')
      UNIQ_MAPPED_DEDUP_SPIKE=$(cat ${SAMPLE}_spikein_deDup.log |
      grep "Number of reads out:" | awk '{print $8}')
      PER_DUPS_SPIKE=$(awk -v a="${UNIQ_MAPPED_DEDUP_SPIKE}" -v b="${UNIQ_MAPPED_SPIKE}" 'BEGIN { print (1-a/b)*100}')%
      echo ${UNIQ_MAPPED_DEDUP_SPIKE}
      echo ${UNIQ_MAPPED_SPIKE}

      echo -e $NAME'\t'\
      $RAW_READS'\t'\
      $TRIMMED_READS'\t'\
      $PER_DIMER'\t'\
      $INSERT_SIZE'\t'\
      $RRNA'\t'\
      $PER_RRNA'\t'\
      $PASSED_FILTERS'\t'\
      $B_CONC'\t'\
      $B_MULTI'\t'\
      $B_UNAL'\t'\
      $B_OAP'\t'\
      $B_CONC_PER'\t'\
      $B_MULTI_PER'\t'\
      $B_UNAL_PER'\t'\
      $UNIQ_MAPPED'\t'\
      $UNIQ_MAPPED_DEDUP'\t'\
      $PER_DUPS'\t'\
      $UNIQ_MAPPED_SPIKE'\t'\
      $UNIQ_MAPPED_DEDUP_SPIKE'\t'\
      $PER_DUPS_SPIKE  >> infoTable.tsv

      done
  '''
}
