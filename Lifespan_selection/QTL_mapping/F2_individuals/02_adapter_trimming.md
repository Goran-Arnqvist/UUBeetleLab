# Adapter trimming

Remove Illumina adapters using [cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html).

## Original scripts
Marek initially shared the script below, which is used after demultiplexing to remove Illumina adapters:

- `A00.Haplo_VX_cutadapt.sge.sh`:
```bash
#!/bin/bash

#$ -pe parallel 10
#$ -l h_vmem=2G
#$ -l 'hostname=!(node524|node525|node514|node526)'
#$ -N cutadaptInfo
#$ -o ./
#$ -j y
#$ -S /bin/bash
#$ -cwd

if [ ! -e "/tmp/mkucka" ];then mkdir /tmp/mkucka;fi

unset PYTHONPATH
export LD_LIBRARY_PATH=/fml/chones/local/lib
export PYTHONPATH=/fml/chones/local/lib/python3.8/site-packages:$PYTHONPATH

file=$1;
echo $file

fbname=$(basename $file _R1_001.fastq.gz)
echo $FBNAME

dir=$(dirname $file)
echo $dir


#### remove DARK-GGGG, trim 5`Tn5ME from Read1; 3`RC-Tn5ME from of Read2, MutTn5ME from 5`Read2, RC-MutTn5ME from 3`Read1, MutTn5ME from 5`Read1
#### then remove reads if one in pair shorter than 30bp
cutadapt --nextseq-trim=10 -g AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT -G ACTTGTGTATAAGAGACAG -a CTGTCTCTTATACACAAGT -g ACTTGTGTATAAGAGACAG \
        --cores=2 -O 15 -m 30 --times=2 \
        -o /tmp/mkucka/$FBNAME\_R1_001.cutadapt.fastq.gz -p /tmp/mkucka/$FBNAME\_R2_001.cutadapt.fastq.gz \
        ./$FBNAME\_R1_001.fastq.gz ./$FBNAME\_R2_001.fastq.gz \
        --too-short-output=/tmp/mkucka/$FBNAME\_R1.tooshort.fastq2.gz --too-short-paired-output=/tmp/mkucka/$FBNAME\_R2.tooshort.fastq2.gz

mv /tmp/mkucka/$FBNAME*.gz ./

```

Then, he sent a message notifying that he found out that his cutadapt trimming was no very efficient at removing stuff (especially it was missing -G AGATGTGTATAAGAGACAG), so he revised it to do much better job at the '5-end of reads. Here is the modified cutadapt trimming script, (you also could just run the second part of the scrip on you existing cutadapted reads to not have to do it again from scratch):

 - `A00.Haplo_VX_cutadapt.sge.sh`:
```bash
file=$1;
echo $file

FBNAME=$(basename $file _R1_001.fastq.gz)
echo $FBNAME

dir=$(dirname $file)
echo $dir

cutadapt --nextseq-trim=20 -A CTGTCTCTTATACACATCT -a CTGTCTCTTATACACATCT -a CTGTCTCTTATACACAAGT -A CTGTCTCTTATACACAAGT \
       --cores=10 -O 15 --times=2 \
       -o $FBNAME\_R1_001.fastq.1.gz -p $FBNAME\_R2_001.fastq.1.gz \
       $FBNAME\_R1_001.fastq.gz $FBNAME\_R2_001.fastq.gz

cutadapt -G ACTTGTGTATAAGAGACAG -g ACTTGTGTATAAGAGACAG -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG \
       --cores=10 -O 19 --times=4 -m 30 \
       -o $FBNAME\_R1_001.cutadapt.fastq.gz -p $FBNAME\_R2_001.cutadapt.fastq.gz \
       $FBNAME\_R1_001.fastq.1.gz $FBNAME\_R2_001.fastq.1.gz \
       --too-short-output=$FBNAME\_R1.tooshort.fastq.gz --too-short-paired-output=$FBNAME\_R2.tooshort.fastq.gz

rm $FBNAME*.fastq.1.gz
```

## Adapted scripts to Uppmax

Make the `file_R1.list` file, containing the full path to the R1 file of all pools (it should be 14 in total):
```bash
find /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads -name *R1_001.fastq.gz -type f > /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/file_R1.list

# Explore file
cat file_R1.list

/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleB_S121_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleI_S128_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleD_S123_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleH_S127_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleN_S133_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleJ_S129_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleM_S132_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleG_S126_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleE_S124_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleF_S125_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleL_S131_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleK_S130_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleC_S122_L004_tag_R1_001.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads/VF-3324-BeetleA_S120_L004_tag_R1_001.fastq.gz

```

Script `02-adapter_trim.sh`:

```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -M snowy
#SBATCH -p core -n 8
#SBATCH -t 1-00:00:00
#SBATCH -J adapter_trim
#SBATCH -e adapter_trim_%J_%A_%a.err
#SBATCH -o adapter_trim_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-14

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads'
FILE_LIST=${OUTPUT_DIR}/file_R1.list

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: ${FILE_PATH}

# Get file name
FBNAME=$(basename $FILE_PATH _R1_001.fastq.gz)
# Get directory path
INPUT_DIR=$(dirname $FILE_PATH)
# Number of cores to use
N_CORES='8'
# Min length
MIN_LENGTH='30'

# Load packages in Uppmax
module load bioinfo-tools samtools/1.14 FastQC/0.11.9

# Go to the directory where the output will be saved
cd $OUTPUT_DIR

## 1. Run cutadapt

# Remove DARK-GGGG, trim 5`Tn5ME from Read1;
# 3`RC-Tn5ME from of Read2, MutTn5ME from 5`Read2, RC-MutTn5ME from 3`Read1, MutTn5ME from 5`Read1,
# then remove reads if one in pair shorter than 30 bp
cutadapt --nextseq-trim=20 -A CTGTCTCTTATACACATCT -a CTGTCTCTTATACACATCT -a CTGTCTCTTATACACAAGT -A CTGTCTCTTATACACAAGT \
--cores=$N_CORES -O 15 --times=2 \
-o $OUTPUT_DIR/$FBNAME\_R1_001.fastq.1.gz -p $OUTPUT_DIR/$FBNAME\_R2_001.fastq.1.gz \
$INPUT_DIR/$FBNAME\_R1_001.fastq.gz $INPUT_DIR/$FBNAME\_R2_001.fastq.gz

cutadapt -G ACTTGTGTATAAGAGACAG -g ACTTGTGTATAAGAGACAG -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG \
--cores=$N_CORES -O 19 --times=4 -m $MIN_LENGTH \
-o $OUTPUT_DIR/$FBNAME\_R1_001.cutadapt.fastq.gz -p $OUTPUT_DIR/$FBNAME\_R2_001.cutadapt.fastq.gz \
$INPUT_DIR/$FBNAME\_R1_001.fastq.1.gz $INPUT_DIR/$FBNAME\_R2_001.fastq.1.gz \
--too-short-output=$OUTPUT_DIR/$FBNAME\_R1.tooshort.fastq.gz --too-short-paired-output=$OUTPUT_DIR/$FBNAME\_R2.tooshort.fastq.gz

# Remove intermediate files
rm $OUTPUT_DIR/$FBNAME*.fastq.1.gz

# cutadapt parameters:
# --nextseq-trim: Quality trimming of reads using two-color chemistry (NextSeq)
# -m length: Discard processed reads that are shorter than LENGTH
# -O: Minimum overlap (default: 3)
# --paired-output, -p: By default, all processed reads, no matter whether they were trimmed or not, are written to
# the output file specified by the -o option (or to standard output if -o was not provided).
# For paired-end reads, the second read in a pair is always written to the file specified by the -p option.
# --times: Number of multiple rounds of adapter removal
# --too-short-output FILE: Instead of discarding the reads that are too short according to -m, write them to FILE (in FASTA/FASTQ format).

## 2. Run fastQC on the trimmed reads file

# Create directory for FastQC results.
if [ -d $OUTPUT_DIR/FastQC_results ]; then echo "Directory /OUTPUT_DIR/FastQC_results exists"; else mkdir $OUTPUT_DIR/FastQC_results; fi

fastqc -t $N_CORES $OUTPUT_DIR/$FBNAME\_R1_001.cutadapt.fastq.gz $OUTPUT_DIR/$FBNAME\_R2_001.cutadapt.fastq.gz -o $OUTPUT_DIR/FastQC_results

```
Submitted batch job 6961978 on cluster snowy
**Runtime: 00-6:00:00** in total (whole array of 14 jobs), but each file separately would run in much less time than that (1:19:09)


Script to get fastQC reports, because I forgot to run it before, `02-fastqc.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 4
#SBATCH -t 5:00:00
#SBATCH -J fastqc
#SBATCH -e fastqc_%J_%A_%a.err
#SBATCH -o fastqc_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-14

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads'
FILE_LIST=${OUTPUT_DIR}/file_R1.list

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: ${FILE_PATH}

# Get file name
FBNAME=$(basename $FILE_PATH _R1_001.fastq.gz)
# Number of cores to use
N_CORES='4'
# Min length
MIN_LENGTH='30'

# Load packages in Uppmax
module load bioinfo-tools FastQC/0.11.9

# Go to the working directory
cd $OUTPUT_DIR

## Run fastQC on the trimmed reads file

# Create directory for FastQC results.
if [ -d $OUTPUT_DIR/FastQC_results ]; then echo "Directory /OUTPUT_DIR/FastQC_results exists"; else mkdir $OUTPUT_DIR/FastQC_results; fi

fastqc -t $N_CORES $OUTPUT_DIR/$FBNAME\_R1_001.cutadapt.fastq.gz $OUTPUT_DIR/$FBNAME\_R2_001.cutadapt.fastq.gz -o $OUTPUT_DIR/FastQC_results

```
Submitted batch job 30013842
**Runtime: 00-00:49:39**

Generate a single FastQC report for all samples using MultiQC:
```bash
# Load program in Uppmax
module load bioinfo-tools MultiQC/1.12

cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/FastQC_results

multiqc .

```
