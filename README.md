# PROseq_alignment

Modified from https://github.com/JAJ256/PROseq_alignment.sh

[![DOI](https://zenodo.org/badge/254700530.svg)](https://zenodo.org/badge/latestdoi/254700530)

This is a pipeline script for aligning paired-end PRO-seq data that has cells of a different species spiked in for normalization, and uses some combination of random UMI sequences on the ligation end of either the 5' or 3' adapter, or both.

## Development and maintenance
 Implemented in nextflow by Elin Axelsson-Ekker. She also maintains the pipeline. Please contact her for bug reports, feature requests etc. (elin.axelsson@gmi.oeaw.ac.at).

## Files and Data needed
 - Unaligned sequencing data, as fastq (.fastq/.fq/.fastq.gz or .fq.gz) files. MUST be Paired-End

## Output files
The pipeline results in the following folders with output files:

- logs_and_QC
  - align: alignment logs for spikein and experiment alignment
  - dedup: deduplication logs for spikein and experiment alignment
  - fastp: trimming logs for each sample and read
  - fastqc: fastqc output for each sample and read
  - rRNA: bowtie2 logs for rRNA filtering step
  - info: "infoTable.tsv" summary of output
- bams
  - dedup: bam from alignment with experimental genome (deduplication)
- bw
  - unnorm: un-normalized bw files


## Organism/annotation
At the moment the pipeline works only for s.pombe with spikein from s. cerevisiae. Contact me if you need other organisms.



## Recommended setup
### -1. Add Hinkskalle to singularity (one time only)
See [How to add Hinkskalle to singularity](#hinkskalle)
### 0. Start a new tmux session or attach to an existing (optional but useful)
See e.g. https://tmuxcheatsheet.com/ and section [Short intro to tmux](#tmux)
### <a name="clone"></a>1. Clone the repo
From  e.g. *your user folder in the lab folder* do:
```
git clone https://ngs.vbcf.ac.at/repo/berger_pipelines/PROalign<>.git
```
you might be asked for user and password, use the ones from your **forskalle** account
### 2. Create a folder with the inputs files (or links to the input files)
e.g:
```
mkdir -p fastq
ln -s /paths/to/files/*.fastq fastq/
```
As long as you use links and do not actually copy the files, this folder can be
created in your user folder on the lab folder, e.g. in the  folder that you cloned or in the folder above.
### 3. Then load the following modules:
**NOTE** On CBE, you can skip loading singularity
```
module load nextflow/21.04.1  #or higher if you prefer
module load singularity/3.6.2
```
**NOTE** if you are in a tmux session were this was already done, this step can be skipped!
### 4. If you have not already, then go to the  folder.
```
cd <>
```



## <a name="run"></a>How to run
**Make sure you have:**
 1. completed the setup steps
 2. that you are on the cluster with the correct modules loaded
 3. moved into the PROalign folder

To use the pipeline on data with UMI on both sides:


```
    nextflow run PROalign.nf -profile slurm --files "/path/to/files/<>" \
    --FIVEP_UMI true \
    --THREEP_UMI true \
    -w /path/to/scratch/dir
```

Here the "/path/to/files/" is the path to the folder you created earlier.
The -w sets the path to where the pipeline should be run and where the intermediate
files will be stored. **It should be on the scratch and not in the lab folder!**

**NOTE** that the files path needs to be within citation marks ("") otherwise
 nextflow will only use one bam file from the folder and not all.



## Parameters
### Required

- --files: path to input files, see also examples below.
- --FIVEP_UMI: five prime umi? false or true
- --THREEP_UMI: three prime umi? false or true

**Examples files parameter**

- **Fastq alt1:**
  - --files "path/to/files/*_{1,2}.fastq"
  - --files "path/to/files/*_{1,2}.fq"
- **Fastq alt2: compressed**
  - --files "path/to/files/*_{1,2}.fastq.gz"
  - --files "path/to/files/*_{1,2}.fq.gz"


In the above examples /path/to/files/ is the path to the folder with the data. **Note the citation marks around the files path!!**
Fastq files have to have the extension .fq or .fastq. (.fq.gz or .fastq.gz )

If the input data is **paired end and in fastq format**, then the first and second reads should be indicated in the file names. Usually the name of the files end with _1 and _2 followed by .fq or .fastq (or .fq.gz/.fastq.gz if compressed).

The **files from the NGS facility** are ending with 'R1_001.fastq.gz' and 'R2_001.fastq.gz'. To use those files set the files flag to:

 --files "path/to/files/\*_R{1,2}_001.fastq.gz"



### With defaults:
- --output: the name of the output folder. [ "results" ]
- --UMI_LEN: [ 6 ]
- --MAPQ: [ 10 ]


**More nextflow options:**

 - to resume add -resume.
 - to run in background add -bg

See also https://www.nextflow.io/docs/latest/cli.html for full list of options.



## <a name="hinkskalle"></a>How to add Hinkskalle to singularity
As the pipeline uses containers from the ngs registry (aka "Hinkskalle"), you
need to do the following (on the CBE cluster):
```
    singularity remote add --no-login hinkskalle singularity.ngs.vbcf.ac.at
    singularity remote use hinkskalle
```

**It's enough if you do this once.** The settings will then be saved a folder
(.singularity) in your home directory and you do not need to worry about it
any more. **UNLESS** you later change the remote registry (singularity remote use).
Then you have to re-run the last line before you start the pipeline.

**NOTE:** If you are using another cluster (e.g. outside of VBC),
you need to make sure Singularity (version 3.4 or higher) is i
nstalled and loaded BEFORE you add the hinklist registry.



##  <a name="tmux"></a>Short intro to tmux
tmux is a very convent tool to use on the cluster. It has many benefits but some of the most obvious are:
1. It lets you set up an environment on the cluster that you can move in an out from as much you want.
E.g. if you have a tmux session called 'proseq' you can load the modules needed to run the pipeline
(nextflow and singularity) once and then by re-attaching to this session the modules are already loaded.
2. Scripts that are running in the session will continue to run even if you logout from the cluster.

See also e.g. https://en.wikipedia.org/wiki/Tmux/

To start a new session type:
  ```
  tmux new -s proseq
  ```
This will create a tmux session with the name proseq and attach you to it.
To detach from the session press *Ctrl + b* and then *d*
To attach to an existing session type
```
  tmux a -t proseq
```
This will put you into the existing session called proseq. Note that the session needs to be created before you can re-attach to it.
