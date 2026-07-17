# Read mapping and summary statistics

<!-- TOC depthFrom:2 depthTo:3 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Best practice recommendations](#best-practice-recommendations)
- [Overview](#overview)
- [1. Read mapping with BWA](#1-read-mapping-with-bwa)
	- [Original scripts](#original-scripts)
	- [Adapted scripts to Uppmax](#adapted-scripts-to-uppmax)
- [2. Split BAM files into separate individuals](#2-split-bam-files-into-separate-individuals)
	- [Original scripts](#original-scripts)
	- [Adapted scripts to Uppmax](#adapted-scripts-to-uppmax)
	- [2.2 Fix sample labeling](#22-fix-sample-labeling)
	- [2.3 Rename BAM files](#23-rename-bam-files)
- [3. Generate read mapping stats](#3-generate-read-mapping-stats)

<!-- /TOC -->

## Best practice recommendations

From the [heliconius repository](https://github.com/evolgenomics/HeliconiusHaplotagging/tree/main/Demultiplexing):

- Since around mid-2021, we and others have found concrete advantages to using barcode-first read mappers like [EMA](https://github.com/arshajii/ema). Specifically, EMA uses BWA's API to place reads, but does a better job of taking linked-reads (or "read clouds") into account. Please see their repo for details.

For this reason, we recommend substituting EMA for the read mapping step. To do so, use our script `ema_prep.sh` to pre-process and sort the reads prior to mapping.

While we find EMA to be less polished than BWA and it currently involves additional overheads, we feel strongly that the improved read mapping, especially in complex regions, is well worth you trouble. Please consider adopting this recommendation in your pipeline.

- Communication in Slack:

>FC
BUT just for your information - *you’re dealing with brood samples. You don’t need the extra precision of EMA. If you’re more comfortable with bwa, the difference is going to be very marginal. It matters for population samples*.
In crosses, the meiotic resolution is going to limit the information you’ll need linkage information for it. Think of it here as having some fun with a fancy datatype and you can go further, but for your main analysis, I think it’s a complete overkill.
So just *do your bwa, but use the -C switch to retain your BX tag*. bwa is far more polished and handles parallelism much smoother. That’d be my suggestion. For fun, you might well get almost basepair level resolution for your crossover junction, though.

>AF
Question:
by mapping reads with BWA and then phasing haplotypes with STITCH, is it enough to get information about which alleles come from the mother and which ones from the father?

>FC
yup
probably v. straight-forward.

>AF
Neat! But then, what is the benefit of using EMA? My understanding is that EMA uses the barcodes to do the read mapping, while BWA does not (edited)

>FC
yup, EMA shines in complex regions. In your case, the meiotic resolution is the limiting factor because you’ll have large chunks of chromosome in linkage. So unless you’re so very unlucky that your genome is plastered with 80% repetitive element, you’re bound to have single-copy “simple” regions in each 1-10Mbp window to focus on. At those loci where mapping is unique and therefore non-ambiguous, EMA does exactly the same job as bwa [EMA uses bwa under-the-hood]
tl/dr: in your dataset you’ll end up paying the mental overhead of doing EMA without getting much marginal benefit.


## Overview

A two-step process:

1. Map reads of the pools (14 in this case) against the reference genome with BWA, using the `-C` switch to retain the `BX` tag (necessary later). So, a single BAM file per pool is created

2. Split BAM files with an awk script that extracts reads of individuals based on the `BX` tag information and then creates separate BAM files for each individual


## 1. Read mapping with BWA

### Original scripts

- Read placement, using the `-C` switch to include the extra comment tags
```bash
bwa mem -C -t 50 <helera1_demo_dir>/Heliconius_erato_demophoon_v1.fa \ $file ${file/R1_001/R2_001} \ -R "@RG\tID:$fbname\tSM:$fbname\tLB:$fbname\tPL:Illumina.HiSeq3000.2x150" | samtools view -bh - > /tmp/mkucka/$fbname.erato.bam

samtools sort \ -@ 50 -l 9 \ -T /tmp/mkucka/$fbname.tmpsort \ -o /tmp/mkucka/$fbname.erato.sorted.bam \ /tmp/mkucka/$fbname.erato.bam
```

- MarkDuplicates, using BX-aware options
```bash
java -Xmx12g -XX:ParallelGCThreads=64 -jar <picard_dir>/picard.jar MarkDuplicates \ I=$file \ O=/tmp/mkucka/$fbname.pMarkdup.bam \ M=$dir/$fbname.pMarkdup.metrics \ CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT
```

### Adapted scripts to Uppmax

- Create the file list:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files

find /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads -name *R1_001.cutadapt.fastq.gz -type f > /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/file_R1.list
```
Explore file
```bash
cat file_R1.list

/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleK_S130_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleC_S122_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleA_S120_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleG_S126_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleJ_S129_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleD_S123_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleE_S124_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleN_S133_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleM_S132_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleL_S131_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleH_S127_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleI_S128_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleB_S121_L004_tag_R1_001.cutadapt.fastq.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/02-clean-reads/VF-3324-BeetleF_S125_L004_tag_R1_001.cutadapt.fastq.gz
```

Script `03-1-read_mapping.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -M snowy
#SBATCH -p core -n 10
#SBATCH -t 3-00:00:00
#SBATCH -J read_map
#SBATCH -e read_map_%J_%A_%a.err
#SBATCH -o read_map_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-14

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files'
REF_FASTA='/proj/snic2020-6-128/private/a_obtectus_poolseq/00-genome/A.obtectus_v2.0.fasta'
FILE_LIST=$OUTPUT_DIR/file_R1.list
PICARD_JAR='/sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar'

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: $FILE_PATH

# Get file name
FBNAME=$(basename $FILE_PATH _R1_001.cutadapt.fastq.gz)
# Get directory path
INPUT_DIR=$(dirname $FILE_PATH)
# Number of cores to use
N_CORES='10'

# Load packages in Uppmax
module load bioinfo-tools bwa/0.7.17 samtools/1.14 picard/2.23.4

# Go to the working directory
cd $OUTPUT_DIR

#### 1) Read placement
echo -e $(date -u) ": Read mapping began..."

# Map reads to the reference genome using BWA with the -C switch to include the extra comment tags
bwa mem -C -t $N_CORES $REF_FASTA \
$INPUT_DIR/$FBNAME\_R1_001.cutadapt.fastq.gz $INPUT_DIR/$FBNAME\_R2_001.cutadapt.fastq.gz \
-R "@RG\tID:$FBNAME\tSM:$FBNAME\tLB:$FBNAME\tPL:Illumina.NovaseqS4.2x150" | samtools view -bh - > $FBNAME.bam

#### 2) Sort reads
echo -e $(date -u) ": Read sorting began..."

# Create directory for temporal files
if [ -d $FBNAME.tmpsort ]; then echo "tmp/ exists in OUTPUT_DIR"; else mkdir $FBNAME.tmpsort; fi

# Sort reads by coordinate using SAMTOOLS
samtools sort \
-@ $N_CORES -l 9 \
-T $FBNAME.tmpsort \
-o $FBNAME.sorted.bam \
$FBNAME.bam && rm $FBNAME.bam
# --compression-level=COMPRESSION_LEVEL: Compression level to use for sorted BAM, from 0 (known as uncompressed BAM in samtools) to 9

#### 3. Mark duplicates
echo -e $(date -u) ": Mark duplicates began..."

# Mark duplicated reads using BX-aware options
java -Xmx60g -XX:ParallelGCThreads=$N_CORES -jar $PICARD_JAR MarkDuplicates \
I=$FBNAME.sorted.bam \
O=$FBNAME.pMarkdup.bam \
M=$FBNAME.pMarkdup.metrics \
CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT && rm $FBNAME.sorted.bam

# Clean up
rm -R $FBNAME.tmpsort

if [ -d Markdup-metrics ]; then echo "Markdup-metrics"; else mkdir Markdup-metrics; fi
mv $FBNAME.pMarkdup.metrics Markdup-metrics

```
Submitted batch job 6962250 on cluster snowy
**Runtime: 01-08:37:37**

Files generated:
```bash
total 1.6T
-rw-rw-r--  1 angela snic2020-6-128  47G Sep 28 02:22 VF-3324-BeetleA_S120_L004_tag.bam
-rw-rw-r--  1 angela snic2020-6-128 3.5M Sep 28 07:31 VF-3324-BeetleA_S120_L004_tag.pMarkdup.bai
-rw-rw-r--  1 angela snic2020-6-128  31G Sep 28 07:31 VF-3324-BeetleA_S120_L004_tag.pMarkdup.bam
-rw-rw-r--  1 angela snic2020-6-128 4.0K Sep 28 07:31 VF-3324-BeetleA_S120_L004_tag.pMarkdup.metrics
-rw-rw-r--  1 angela snic2020-6-128  29G Sep 28 03:38 VF-3324-BeetleA_S120_L004_tag.sorted.bam
drwxrwsr-x  2 angela snic2020-6-128  28K Sep 28 03:38 VF-3324-BeetleA_S120_L004_tag.tmpsort/
-rw-rw-r--  1 angela snic2020-6-128  53G Sep 29 04:39 VF-3324-BeetleB_S121_L004_tag.bam
-rw-rw-r--  1 angela snic2020-6-128 3.6M Sep 29 10:11 VF-3324-BeetleB_S121_L004_tag.pMarkdup.bai
-rw-rw-r--  1 angela snic2020-6-128  35G Sep 29 10:11 VF-3324-BeetleB_S121_L004_tag.pMarkdup.bam
-rw-rw-r--  1 angela snic2020-6-128 3.9K Sep 29 10:11 VF-3324-BeetleB_S121_L004_tag.pMarkdup.metrics
-rw-rw-r--  1 angela snic2020-6-128  32G Sep 29 06:08 VF-3324-BeetleB_S121_L004_tag.sorted.bam
drwxrwsr-x  2 angela snic2020-6-128  32K Sep 29 06:08 VF-3324-BeetleB_S121_L004_tag.tmpsort/
-rw-rw-r--  1 angela snic2020-6-128  49G Sep 28 03:14 VF-3324-BeetleC_S122_L004_tag.bam
-rw-rw-r--  1 angela snic2020-6-128 3.5M Sep 28 08:36 VF-3324-BeetleC_S122_L004_tag.pMarkdup.bai
-rw-rw-r--  1 angela snic2020-6-128  32G Sep 28 08:36 VF-3324-BeetleC_S122_L004_tag.pMarkdup.bam
-rw-rw-r--  1 angela snic2020-6-128 4.0K Sep 28 08:36 VF-3324-BeetleC_S122_L004_tag.pMarkdup.metrics
-rw-rw-r--  1 angela snic2020-6-128  30G Sep 28 04:43 VF-3324-BeetleC_S122_L004_tag.sorted.bam
```


## 2. Split BAM files into separate individuals

Now that all reads per pool are mapped against the reference genome, it is time to split the BAMs into the individual samples. For this, Frank Chan kindly shared a code to do this.

### Original scripts

What the code does in "English" is:
- Save the header for later binary processing of SAMs
- Read the bulk file and process per-line, and split each line into a separate file based on the sample ID that’s identified on the BX tag itself (with the variant here that I’m saving space by putting them through a gzip binary temporary file)
- Then loop through all the files generated from this while adding the header back in for proper BAM processing with samtools.

```bash
# Save header first
samtools view <input.bam> -H > common.header

# Then process individual files depending on the barcode segment type

#if using C-segment:
samtools view <input.bam> | awk '{match($0, "BX:Z:"); sample_id=substr($0,RSTART+8,3)"-"substr($0,RSTART+18,4); print $0"" | "gzip -c > "sample_id".sam.gz"}'
#if using D-segment:
samtools view <input.bam> | awk '{match($0, "BX:Z:"); sample_id=substr($0,RSTART+14,3)"-"substr($0,RSTART+18,4); print $0"" | "gzip -c > "sample_id".sam.gz"}'

# Then reformat these into BAM files
#test for 5 lines
for sample in `ls *.sam.gz | head -5`; do stem=`basename $sample`; cat common.header <(gzip -d $sample -c) | samtools view - -O BAM -o $stem.bam; done
# If it works, then change the bash loop to *.sam.gz to run everything:
for sample in *.sam.gz; do stem=`basename $sample`; cat common.header <(gzip -d $sample -c) | samtools view - -O BAM -o $stem.bam; done

# Tidy up the file names
for i in *.sam.gz.bam; do mv $i ${i/.sam.gz/}; done

# Once you've checked them all and are happy with the results, clean up:
rm *.sam.gz

```

### Adapted scripts to Uppmax

- Create a `bam.list` file:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data

find /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files -name *_L004_tag.pMarkdup.bam -type f > /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/pool_bam.list
```
Explore file:
```bash
cat /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/pool_bam.list

/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleD_S123_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleF_S125_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleB_S121_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleE_S124_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleJ_S129_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleK_S130_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleA_S120_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleG_S126_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleL_S131_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleC_S122_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleH_S127_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleM_S132_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleI_S128_L004_tag.pMarkdup.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/VF-3324-BeetleN_S133_L004_tag.pMarkdup.bam
```

Script `03-2-split_BAMs.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 2
#SBATCH -t 1-00:00:00
#SBATCH -J split_BAMs
#SBATCH -e split_BAMs_%J_%A_%a.err
#SBATCH -o split_BAMs_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-14

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample'
FILE_LIST=$OUTPUT_DIR/pool_bam.list

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: $FILE_PATH

# Get file name
FBNAME=$(basename $FILE_PATH _L004_tag.pMarkdup.bam)
# Get directory path
INPUT_DIR=$(dirname $FILE_PATH)

# Load programs in Uppmax
module load bioinfo-tools samtools/1.14

# Make a directory for the plate
mkdir $OUTPUT_DIR/$FBNAME

# Go to the working directory
cd $OUTPUT_DIR/$FBNAME

#### Extract reads per individual using the X tags and then create separate BAM files

# Save header first
samtools view $INPUT_DIR/$FBNAME\_L004_tag.pMarkdup.bam -H > $FBNAME.common.header

# Process individual files depending on the barcode segment, D-segment in the beetles data set
samtools view $INPUT_DIR/$FBNAME\_L004_tag.pMarkdup.bam | awk '{match($0, "BX:Z:"); sample_id=substr($0,RSTART+14,3)"-"substr($0,RSTART+18,4); print $0"" | "gzip -c > "sample_id".sam.gz"}'

# Reformat these mapped reads into BAM files
for sample in *.sam.gz; do stem=`basename $sample`; cat $FBNAME.common.header <(gzip -d $sample -c) | samtools view - -O BAM -o $stem.bam; done

# Tidy up the file names
for i in *.sam.gz.bam; do mv $i ${i/.sam.gz/}; done

# Create an index for the final bam file.
for i in *.bam; do samtools index -@ 2 $i; done

# Remove intermediate files
rm *.sam.gz
rm $FBNAME.common.header

# Move all bam files to the working directory
mv *.bam* ../

# Delete temporary directory
#cd ..
#rm -R $FBNAME

```
Submitted batch job 30089569 > job 1, run for testing
**Runtime: 00-06:38:46**
Submitted batch job 30096109 > jobs 2-14
**Runtime: 00-10:21:49**

All runs ended successfully.

I noticed that some of the BAM files were "empty" or failed:
```bash
VF-3324-BeetleD_S123$ for i in *.bam; do samtools index -@ 2 $i; done
[E::hts_hopen] Failed to open file D01-N704.bam
[E::hts_open_format] Failed to open file "D01-N704.bam" : Exec format error
samtools index: failed to open "D01-N704.bam": Exec format error
```
Thus, I saved in a file a list of all the samples per plate for which an index file was and not generated:
```bash
ll */ > bam_bai_files_generated.txt
```

Count the number of BAM files, it should be equal to the number of samples (96 samples per plate x 14 plates = 1344s):
```bash
ls D*.bam | wc -l
1358
```
Hmmm. A total of 97 BAM files were generated per plate with codes from D00 to D96, but there should be 96 (as per 96-well plates). Frank clarified that `D00` corresponds to ambiguous reads, equivalent to barcode assignments that are not clearly assigned to a single D01-D96 barcode. So, he recommends to leave those the D00 files out (I moved them to a separate dir):
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample

mkdir ambiguous-D00

ls *bam | wc -l
#1358

for file in `ls *D00-*`; do mv $file ambiguous-D00/; done

ls *bam | wc -l

#1344
```

### 2.2 Fix sample labeling

I noticed that there are 5 samples that have the same ID, which means they potentially were sequenced twice:
```bash
plate_barcode	sample
D35-N711	AB7_26
D23-N714	AB7_26
D85-N709	BA12_2
D92-N710	BA12_2
D19-N704	BA7_27
D36-N707	BA7_27
D44-N706	DC12_3
D19-N710	DC12_3
D75-N705	DC2_13
D18-N711	DC2_13
```
I asked Johanna and Göran about how to treat these samples (merge them, or fix typos?). Johanna revised her files and identified some typos (mislabeled samples) and a dubious sample (DC 12.3, repeated twice):
```bash
Row 2 = AB 7.26
Row 3 = AB 17.26
Row 4 = BA 12.2
Row 5 = BA 2.12
Row 6 = BA 5.27
Row 7 = BA 7.27
Row 8 = DC 12.3
Row 9 = DC 12.3
Row 10 = CD 2.13
Row 11 = DC 2.13
```
Göran indicated that four of five sorted out then - due to lab book to Excel file typos. Good! For the duplicate of DC 12.3, I checked and Rebecca did not score two individuals of DC 12.3. This means that this is most likely due to a "typo" by Johanna when writing down the id's in the handwritten notes from the vial labels. So, we should ditch these two samples as we cannot be absolutely sure what phenotype to associate with them (this would be true even if we would find e.g. an individual with a similar label present in the phenotype data but absent in the seq samples).



### 2.3 Rename BAM files

I made a back-up copy of all BAM files with the original naming in the directory `original_bams_TO_DELETE`:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample

ls *.bam | wc -l
1344

mkdir original_bams_TO_DELETE
cp *.bam* original_bams_TO_DELETE

cd original_bams_TO_DELETE
ls *.bam | wc -l
1344

cd ..
```
**NOTE:** The back up copy of the BAM files should be deleted once the renaming has been proven to be succesful.

Removed samples `D44-N706.bam`and `D19-N710.bam` as their labeling is ambigous (they are called the same `DC12-3`, so they cannot tell appart - phebotype):
```bash
mv D44-N706.bam* ambiguous-D00/
mv D19-N710.bam* ambiguous-D00/

ls *.bam | wc -l
1342
```
These files were removed from the renaming file. The `rename_bam_files_key.txt` looks like:
```bash
head rename_bam_files_key.txt
D35-N711	AB7_26.D35_N711
D23-N714	AB17_26.D23_N714
D85-N709	BA12_2.D85_N709
D92-N710	BA2_12.D92_N710
D19-N704	BA5_27.D19_N704
D36-N707	BA7_27.D36_N707
D75-N705	CD2_13.D75_N705
D18-N711	DC2_13.D18_N711
D01-N701	BA6_8.D01_N701
D02-N701	BA6_14.D02_N701

cat rename_bam_files_key.txt | wc -l
1342
```

Test run of the bash commands to rename files:
```bash
cp rename_bam_files_key.txt test_rename/
cd test_rename

while IFS='' read -r line || [[ -n "$line" ]]; do
	NR=$(echo "$line" | cut -f 1)
	NAME=$(echo "$line" | cut -f 2)
	if [ -f "${NR}.bam" ]; then
		mv "${NR}.bam" "${NAME}.bam"
		mv "${NR}.bam.bai" "${NAME}.bam.bai"
	fi
done < rename_bam_files_key.txt
```

As the script worked fine, I run it for all the files, but first created the script `rename_bam_files_with_key.sh` with nano:
```bash
#!/bin/bash

# Rename files using a text file with two columns separated by a TAB as reference
while IFS='' read -r line || [[ -n "$line" ]]; do
	NR=$(echo "$line" | cut -f 1)
	NAME=$(echo "$line" | cut -f 2)
	if [ -f "${NR}.bam" ]; then
		mv "${NR}.bam" "${NAME}.bam"
		mv "${NR}.bam.bai" "${NAME}.bam.bai"
	fi
done < rename_bam_files_key.txt
```

Then, run it with:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample
bash rename_bam_files_with_key.sh
```
One error message:
```bash
mv: cannot stat ‘D01-N704.bam.bai’: No such file or directory
```
I confirmed the bam file exists but not the index file:
```bash
ll *D01_N704*
-rw-rw-r-- 1 angela snic2020-6-128 1.3M Sep 29 20:30 DC12_2.D01_N704.bam
```
Very tiny file, maybe this explains why there is not index file, I will continue with the read map stats.


## 3. Generate read mapping stats

- Create a `bam_sample.list` file, while excluding the bam files from other subdirectories:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample
ls *.bam | wc -l
1342

cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data

# Create a text file listing all the sample bam files, while excluding bam files within subdirectories
find /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample -name *.bam \
! -path "/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/ambiguous-D00/*" \
! -path "/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/test_rename/*" \
! -path "/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/original_bams_TO_DELETE/*" \
-type f > /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/bam_sample.list
```

Explore the generated file:
```bash
head /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/bam_sample.list
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA7_33.D05_N701.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/AB7_18.D29_N709.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA12_28.D02_N703.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA8_7.D19_N712.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC10_23.D82_N714.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/AB7_19.D40_N705.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/CD2_18.D06_N708.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/CD11_26.D53_N706.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA12_15.D61_N712.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/AB1_29.D47_N705.bam

# Verify the number of rows is correct
cat /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/bam_sample.list | wc -l
1342  # good! This is the expected number of files for (14 plates x 96) - 2 = 1342
```

Script `03-3_qualimap_1-1000.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 4
#SBATCH -t 20:00:00
#SBATCH -J qualimap
#SBATCH -e qualimap_%J_%A_%a.err
#SBATCH -o qualimap_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-1000

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample'
FILE_LIST=${OUTPUT_DIR}/bam_sample.list

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)
# Get file name
FBNAME=$(basename $FILE_PATH .bam)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: ${FILE_PATH}

#### Generate mapping quality summary statistics

# Setup packages in Uppmax
module load bioinfo-tools QualiMap/2.2.1
unset DISPLAY  # Turn display off to avoid problems with X11 in Uppmax

if [ -d $OUTPUT_DIR/Qualimap_stats ]; then echo "Qualimap_stats/ exists"; else mkdir $OUTPUT_DIR/Qualimap_stats; fi

# Go to the working directory
cd $OUTPUT_DIR

qualimap --java-mem-size=24G bamqc -bam $FILE_PATH -ip -outdir $OUTPUT_DIR/Qualimap_stats/$FBNAME -outformat html

if [ -d $OUTPUT_DIR/Samtools_stats ]; then echo "Samtools_stats/ exists"; else mkdir $OUTPUT_DIR/Samtools_stats; fi

module load samtools/1.14

# Generate coverage statistics with samtools
samtools stats -@ 4 $FILE_PATH > $OUTPUT_DIR/Samtools_stats/$FBNAME.stat

```
Submitted batch job 30260684
**Runtime: 00-04:03:00**


# 2022-Oct-7

The job array `#SBATCH -a 1001-1342` did not work in slurm. Thus, I created a new file listing the bam files from 1001 to 1342:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample

bam_sample.list

sed -n '1001,1342 p' bam_sample.list > bam_sample_1001-1342.list
#sed -n '1,1000 p' bam_sample.list > bam_sample_1-1000.list

# Verify if the number of lines in the file is the one expected: 342
wc -l bam_sample_1001-1342.list
342

# Check which file is in position 1001
head -n 1001 bam_sample.list
...
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA6_24.D59_N701.bam

# Verify this file is the first one in the subsetted file
head bam_sample_1001-1342.list
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA6_24.D59_N701.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC14_34.D59_N707.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/AB10_33.D39_N702.bam
```

Script `03-3_qualimap_1001-1342.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 4
#SBATCH -t 20:00:00
#SBATCH -J qualimap
#SBATCH -e qualimap_%J_%A_%a.err
#SBATCH -o qualimap_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-342

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample'
FILE_LIST=${OUTPUT_DIR}/bam_sample_1001-1342.list

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)
# Get file name
FBNAME=$(basename $FILE_PATH .bam)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: ${FILE_PATH}

#### Generate mapping quality summary statistics

# Setup packages in Uppmax
module load bioinfo-tools QualiMap/2.2.1
unset DISPLAY  # Turn display off to avoid problems with X11 in Uppmax

if [ -d $OUTPUT_DIR/Qualimap_stats ]; then echo "Qualimap_stats/ exists"; else mkdir $OUTPUT_DIR/Qualimap_stats; fi

# Go to the working directory
cd $OUTPUT_DIR

qualimap --java-mem-size=24G bamqc -bam $FILE_PATH -ip -outdir $OUTPUT_DIR/Qualimap_stats/$FBNAME -outformat html

if [ -d $OUTPUT_DIR/Samtools_stats ]; then echo "Samtools_stats/ exists"; else mkdir $OUTPUT_DIR/Samtools_stats; fi

module load samtools/1.14

# Generate coverage statistics with samtools
samtools stats -@ 4 $FILE_PATH > $OUTPUT_DIR/Samtools_stats/$FBNAME.stat

```
**Runtime: 00-00:14:00**
Submitted batch job 30316710

Examine error messages:
```bash
head qualimap_*.err
==> qualimap_30260684_30260684_1000.err <==

==> qualimap_30260858_30260684_1.err <==
-------------------------------------------------------------------------------
The following dependent module(s) are not currently loaded: java/sun_jdk1.8.0_92 (required by: R/3.3.2)
-------------------------------------------------------------------------------


==> qualimap_30317604_30316710_46.err <==
Failed to run bamqc
java.lang.RuntimeException: The BAM file is empty or corrupt
	at org.bioinfo.ngs.qc.qualimap.process.BamStatsAnalysis.run(BamStatsAnalysis.java:520)
	at org.bioinfo.ngs.qc.qualimap.main.BamQcTool.execute(BamQcTool.java:242)
	at org.bioinfo.ngs.qc.qualimap.main.NgsSmartTool.run(NgsSmartTool.java:190)
	at org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain.main(NgsSmartMain.java:111)

[E::hts_hopen] Failed to open file /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC12_2.D01_N704.bam
[E::hts_open_format] Failed to open file "/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC12_2.D01_N704.bam" : Exec format error
samtools stats: failed to open "/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC12_2.D01_N704.bam": Exec format error
```


Generate a single Qualimap report for all samples using MultiQC:
```bash
# Load program in Uppmax
module load bioinfo-tools MultiQC/1.12

cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/Qualimap_stats

multiqc .

```
Stdout:
```
|           multiqc | MultiQC Version v1.13 now available!
|           multiqc | Search path : /crex/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/Qualimap_stats
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 67050/67050  
|          qualimap | Found 1341 BamQC reports
|           multiqc | Compressing plot data
|           multiqc | Report      : multiqc_report.html
|           multiqc | Data        : multiqc_data
|           multiqc | MultiQC complete
|           multiqc | 4 flat-image plots used in the report due to large sample numbers
|           multiqc | To force interactive plots, use the '--interactive' flag. See the documentation.
```
