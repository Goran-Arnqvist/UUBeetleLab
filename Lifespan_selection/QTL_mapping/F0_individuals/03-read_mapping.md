# Mapping clean reads to the Pacbio bean beetle genome assembly

# 2022-05-20

## Create a copy of genome files into the new parentals/ directory

I used these commands:
```bash
cd <REDACTED_PATH>

ln -s <REDACTED_PATH> .
```

## Setup input files and scripts for read mapping

I will use the read mapper [bwa](http://bio-bwa.sourceforge.net) to align clean reads against the reference genome. For this, I created various input files required by the script.

- Path to R1 files:
```bash
ls -lh <REDACTED_PATH> | awk '{print $9}' > <REDACTED_PATH>

cat <REDACTED_PATH>

<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
```
- Path to R2 files:
```bash
ls -lh <REDACTED_PATH> | awk '{print $9}' > <REDACTED_PATH>

cat <REDACTED_PATH>

<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
```
- Sample names:
```bash
ls -lh <REDACTED_PATH> | awk '{print $9}' | sed -e 's+/+\t+g' | cut -f10 | sed 's/_R.*$//' | uniq > <REDACTED_PATH>

cat <REDACTED_PATH>

Sample_VA-3193-EL1-F
Sample_VA-3193-EL1-M
Sample_VA-3193-EL4-F
Sample_VA-3193-EL4-M
Sample_VA-3193-LE1-F
Sample_VA-3193-LE1-M
Sample_VA-3193-LE2-F
Sample_VA-3193-LE2-M
```

Script `03-1-read_mapping.sh` that launches 16 separate jobs (one for each sample):

```bash
#!/bin/bash
#SBATCH -A <REDACTED_PROJECT_CODE>
#SBATCH -M <REDACTED_HPC_SITE>
#SBATCH -p core -n 8
#SBATCH -t 4-00:00:00
#SBATCH -J map_reads
#SBATCH -e map_reads_%J_%A_%a.err
#SBATCH -o map_reads_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=<REDACTED_EMAIL>
#SBATCH -a 1-8

# Load required software.
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.4

# Set environment variables.
# Files.
WORK_DIR='<REDACTED_PATH>'
R1_FILE='./R1_files_bwa.txt'
R2_FILE='./R2_files_bwa.txt'
ID_FILE='./sample_IDs_bwa.txt'
REF_FILE='<REDACTED_PATH>'
# Programs.
PICARD_JAR='<REDACTED_PATH>'

# Directory where the BAM files will be stored.
cd $WORK_DIR

# For a given job in the array, set the correspondent sample files.
R1_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $R1_FILE)
R2_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $R2_FILE)
ID_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_FILE)

# Print current file info to stdout for future reference.
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the R1 target: ${R1_target}
echo -e This is the R2 target: ${R2_target}
echo -e This is the ID target: ${ID_target}
echo -e This is the REF genome: ${REF_FILE}

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

echo -e $(date -u) ": Read mapping began..."

# Map reads to reference genome using several threads [-t 8] and mark split alignment [-M]; then sort reads and generate the bam file index.
bwa mem -M -t 8 -R '@RG\tID:${ID_target}\tSM:${ID_target}' ${REF_FILE} ${R1_target} ${R2_target} | samtools view -@ 8 -b -S - > $SNIC_TMP/${ID_target}.bam
samtools sort -@ 8 -T $SNIC_TMP/tmp -o $SNIC_TMP/${ID_target}.sort.bam $SNIC_TMP/${ID_target}.bam && rm $SNIC_TMP/${ID_target}.bam
samtools index -@ 8 $SNIC_TMP/${ID_target}.sort.bam

echo -e $(date -u) ": Read mapping and bam file indexing ended..."

# Mark duplicate reads.
if [ -d $WORK_DIR/MarkDup_metrics ]; then echo "MarkDup_metrics/ exists"; else mkdir $WORK_DIR/MarkDup_metrics; fi
echo -e $(date -u) ": Mark duplicates began..."

java -Xmx60G -jar ${PICARD_JAR} MarkDuplicates \
I=$SNIC_TMP/${ID_target}.sort.bam O=$SNIC_TMP/${ID_target}.sort.MarkDup.bam M=$WORK_DIR/MarkDup_metrics/${ID_target}.MarkDup.txt \
ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=$SNIC_TMP/tmp && rm $SNIC_TMP/${ID_target}.sort.bam

echo -e $(date -u) ": Mark duplicates ended..."
echo -e $(date -u) ": Add read groups began..."

# Add read groups to the bam file.
java -Xmx60G -jar ${PICARD_JAR} AddOrReplaceReadGroups \
I=$SNIC_TMP/${ID_target}.sort.MarkDup.bam O=$SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam \
RGID=${ID_target} RGLB=${ID_target} RGPL=illumina RGPU=${ID_target} RGSM=${ID_target} && rm $SNIC_TMP/${ID_target}.sort.MarkDup.bam

# Create an index for the final bam file.
samtools index -@ 8 $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam

echo -e $(date -u) ": Add read groups and bam index creation ended..."

# Copy final BAM file and its index to the Results directory.
cp $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam $WORK_DIR
cp $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam.bai $WORK_DIR

# Obtain mapping quality summary statistics.
module load QualiMap/2.2.1
unset DISPLAY  # Turn display off to avoid problems with X11 in <REDACTED_HPC_SITE>.

if [ -d $WORK_DIR/Qualimap_results ]; then echo "Qualimap_results/ exists"; else mkdir $WORK_DIR/Qualimap_results; fi

echo -e $(date -u) ": Qualimap began..."

qualimap --java-mem-size=60G bamqc -bam $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam -ip -outdir $WORK_DIR/Qualimap_results/${ID_target} -outformat html

echo -e $(date -u) ": Qualimap ended..."

```
Submitted batch job <REDACTED_JOB_ID> on cluster <REDACTED_HPC_SITE>
**Runtime: 01-03:56:00**??


# 2022-05-22

## Generate a Multiqc report of the qualimap reports
```
cd <REDACTED_PATH>

module load bioinfo-tools MultiQC/1.12
multiqc .

/// MultiQC 🔍 | v1.12

|           multiqc | Search path : <REDACTED_PATH>
|         searching | ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 400/400  
|          qualimap | Found 8 BamQC reports
|           multiqc | Compressing plot data
|           multiqc | Report      : multiqc_report.html
|           multiqc | Data        : multiqc_data
|           multiqc | MultiQC complete

```
