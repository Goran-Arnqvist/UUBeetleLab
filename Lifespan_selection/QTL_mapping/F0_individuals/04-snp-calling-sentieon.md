# SNP calling using GATK implemented with the Sentieon pipelines

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [About sentieon pipelines](#about-sentieon-pipelines)
- [Instructions from Mafalda](#instructions-from-mafalda)
- [Run pipeline on F0 beetles](#run-pipeline-on-f0-beetles)

<!-- /TOC -->

## About sentieon pipelines

Sentieon® develops and supplies a suite of bioinformatics secondary analysis tools that process genomics data with high computing efficiency, fast turnaround time, exceptional accuracy, and 100% consistency. More information in this website: https://support.sentieon.com/manual/, https://support.sentieon.com/quick_start/

https://github.com/Sentieon/sentieon-scripts/blob/master/example_pipelines/germline/joint-calling.sh
https://support.sentieon.com/manual/examples/examples/#dna-pipeline-example-script
https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller

# 2022-Oct-26

## Run Sentieon pipeline on F0 beetles

- Upload licenses to Uppmax

- Upload the latest version of the software package for Linuz shared by Don Freed, don.freed@sentieon.com (https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-202112.05.tar.gz)

- Initiate license in Uppmax (rackham2):
```bash
# Create an empty file for the logs
touch gatk_OTT_parentals_Oct-26-2022.log

# Set environment variables
# Make sure to replace /home/ for /domus/h1/ when setting the path in Uppmax
SENTIEON_DIR='/domus/h1/angela/sentieon-pipelines/sentieon-genomics-202112.05'
SENTIEON_LICENSE='/domus/h1/angela/sentieon-pipelines/licenses/Uppsala_University_node-51.lic'
export SENTIEON_LICENSE='/domus/h1/angela/sentieon-pipelines/licenses/Uppsala_University_node-51.lic'
#export SENTIEON_LICENSE=10.0.10.51:8990
LOG_FILE='/domus/h1/angela/sentieon-pipelines/logs/gatk_OTT_parentals_Oct-26-2022.log'

# Initialize the server (--stop to make a running server stops)
$SENTIEON_DIR/bin/sentieon licsrvr --start --log $LOG_FILE $SENTIEON_LICENSE

# If you use the `--log`, you can `cat` the file to make sure the license has started. It should say something like:
cat $LOG_FILE
# 2022/10/26 20:13:58 Start serving on 10.0.10.51:8990 (sentieon-genomics-202112.05)

2022/10/26 20:13:58 Start serving on 10.0.10.51:8990 (sentieon-genomics-201911)
2022/10/27 17:46:24 cmdline: /domus/h1/angela/sentieon-pipelines/sentieon-genomics-202112.05/libexec/licsrvr --start --log /domus/h1/angela/sentieon-pipelines/logs/gatk_OTT_parentals_Oct-26-2022.log /domus/h1/angela/sentieon-pipelines/licenses/Uppsala_University_node-51.lic
2022/10/27 17:46:24 Start serving on 10.0.10.51:8990 (sentieon-genomics-202112.05)
2022/10/27 17:46:25 Failed to bind to address 10.0.10.51: Address already in use
2022/10/27 17:46:35 Stop serving
2022/10/27 17:51:10 Stop serving
2022/10/27 17:52:22 cmdline: /domus/h1/angela/sentieon-pipelines/sentieon-genomics-202112.05/libexec/licsrvr --start --log /domus/h1/angela/sentieon-pipelines/logs/gatk_OTT_parentals_Oct-26-2022.log /domus/h1/angela/sentieon-pipelines/licenses/Uppsala_University_node-51.lic
2022/10/27 17:52:22 Received update from the master server
2022/10/27 17:52:22 Start serving on 10.0.10.51:8990 (sentieon-genomics-202112.05)
```

- Run HaplotypeCaller

Make sure there is a file listing the BAM files of interest:
```bash
cat /proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/parentals_bam.list

/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-EL4-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-EL1-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-EL4-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-LE1-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-LE2-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-LE2-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-LE1-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/Sample_VA-3193-EL1-F.sort.MarkDup.RG.bam
```
Script `04-sentieon_01_HaplotypeCaller.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 5
#SBATCH -t 2-00:00:00
#SBATCH -J sentieon
#SBATCH -e sentieon_%J_%A_%a.err
#SBATCH -o sentieon_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-8

# Load packages in Uppmax
module load python/2.7.15

# *******************************************
# Script to perform DNA seq variant calling
# adapted from senteion original script
# *******************************************

# Update with the fullpath location of your sample fastq
set -x
BAM_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files'
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files'

# Update with the location of the reference data files
REF_FASTA="/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/00-genome/A.obtectus_v2.0.fasta"

# Set SENTIEON_LICENSE if it is not set in the environment
export SENTIEON_LICENSE='/domus/h1/angela/sentieon-pipelines/licenses/Uppsala_University_node-51.lic'
echo $SENTIEON_LICENSE

# Update with the location of the Sentieon software package
SENTIEON_INSTALL_DIR='/domus/h1/angela/sentieon-pipelines/sentieon-genomics-202112.05'

# Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=$SNIC_TMP

# It is important to assign meaningful names in actual cases.
# It is particularly important to assign different read group names.
BAM_LIST='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/03-bam-files/parentals_bam.list'
BAM_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $BAM_LIST)
#BAM=$(ls *_2pass_star_resultsAligned.sortedByCoord.out.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
SAMPBASE=$(basename $BAM_target)
SAMPLE=${SAMPBASE%.sort.MarkDup.RG.bam*}
SAMPLE_DIR=$(dirname $BAM_target)

echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the SAMPBASE: $SAMPBASE
echo -e This is the SAMPLE: $SAMPLE
echo -e This is the SAMPLE_DIR: $SAMPLE_DIR
#SAMPLE="sprat"

# Other settings
N_THREADS='5' # number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
workdir=$OUTPUT_DIR
#mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir

# Sentieon proprietary compression
bam_option="--bam_compression 1"

# Run GATK 4.1
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REF_FASTA -t $N_THREADS -i $BAM_target --algo Haplotyper --genotype_model multinomial --emit_mode gvcf --emit_conf 30 --call_conf 30 ${SAMPLE}.g.vcf.gz
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REF_FASTA -t $N_THREADS -i $BAM_DIR/${SAMPLE}.RG.bam --algo Haplotyper --genotype_model multinomial --emit_mode gvcf --emit_conf 30 --call_conf 30 ${SAMPLE}.g.vcf.gz

```
Submitted batch job 30696808
**Runtime: 00-08:53:00**

- Run Genotyper:

Make a file listing the gvcf files:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files
cd ..

find /proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files -name *.g.vcf.gz -type f > /proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/gvcf_list.txt

cat gvcf_list.txt

/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-LE1-M.g.vcf.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-EL4-F.g.vcf.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-EL1-F.g.vcf.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-LE2-M.g.vcf.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-EL1-M.g.vcf.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-LE1-F.g.vcf.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-EL4-M.g.vcf.gz
/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/Sample_VA-3193-LE2-F.g.vcf.gz

wc -l gvcf_list.txt
8
```
Script `04-sentieon_02_Genotyper.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p node
#SBATCH -t 1-00:00:00
#SBATCH -J genotyper
#SBATCH -e genotyper_%A_%a.err
#SBATCH -o genotyper_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com

# Load packages in Uppmax
module load python/2.7.15

# Set environmental variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files'
REF_FASTA='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/00-genome/A.obtectus_v2.0.fasta'
GVCF_LIST='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files/gvcf_list.txt'
OUTPUT_PREFFIX='OTT_8ind_F0_sentieon'

# sentieon options
export SENTIEON_LICENSE='/domus/h1/angela/sentieon-pipelines/licenses/Uppsala_University_node-51.lic'
SENTIEON_INSTALL_DIR='/domus/h1/angela/sentieon-pipelines/sentieon-genomics-202112.05'
SENTIEON_TMPDIR=$SNIC_TMP

# Go to the work directory
cd $OUTPUT_DIR

# Run sentienion
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REF_FASTA -t 20 --algo GVCFtyper --emit_mode CONFIDENT $OUTPUT_PREFFIX.vcf.gz $(<$GVCF_LIST)

```
Submitted batch job 30704989
**Runtime: 00-01:06:57**


- Run filtering

Script `04-sentieon_03_Filtering.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 5
#SBATCH -t 2-00:00:00
#SBATCH -J filter_vcf
#SBATCH -e filter_vcf_%J_%A_%a.err
#SBATCH -o filter_vcf_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com

# Load packages in Uppmax
module load bioinfo-tools samtools vcftools bcftools picard/2.20.4
module load GATK/4.1.4.1

# Set environmental variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/04-vcf-files'
REF_FASTA='/proj/snic2020-6-128/private/a_obtectus_QTLmap/parentals/data/00-genome/A.obtectus_v2.0.fasta'
OUTPUT_PREFFIX='OTT_8ind_F0_sentieon'

# Go to the working directory
cd $OUTPUT_DIR

# Run GATK
gatk --java-options "-Xmx30g" VariantFiltration  \
    -R $REF_FASTA \
    -V $OUTPUT_PREFFIX.vcf.gz \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0))" \
    --filter-name "RPRS8" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('QD') && QD < 2.0))" \
    --filter-name "QD2" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('FS') && FS > 60.0))" \
    --filter-name "FS60" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('SOR') && SOR > 3.0))" \
    --filter-name "SOR3" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQ') && MQ < 40.0))" \
    --filter-name "MQ40" \
    --filter-expression "(vc.isSNP() && (vc.hasAttribute('MQRankSum') && MQRankSum < -12.5))" \
    --filter-name "MQ12.5" \
    -G-filter "vc.isSNP() && DP < 3" \
    -G-filter-name "gtDP3" \
    -G-filter "vc.isSNP() && GQ < 20" \
    -G-filter-name "gtGQ20" \
    -O $OUTPUT_PREFFIX.SV.VF.F2.vcf.gz

# set filtered GT to no call
gatk --java-options "-Xmx30g" SelectVariants  \
    -R $REF_FASTA \
    -V $OUTPUT_PREFFIX.SV.VF.F2.vcf.gz \
    --set-filtered-gt-to-nocall \
    -O $OUTPUT_PREFFIX.SV.VF.F2.setGT.vcf.gz

# Filter on sample specific Max Depth on bcftools. I have a final filter for INFO/DP < 3 because the
# command above does not filter ref calls I think.
bcftools filter -S . -e "FMT/DP[0] > 30.9" $OUTPUT_PREFFIX.SV.VF.F2.setGT.vcf.gz | bcftools filter -S . -e "FMT/DP[1] > 36.3" | bcftools filter -S . -e "FMT/DP[2] > 36" | bcftools filter -S . -e "FMT/DP[3] > 36.6" | bcftools filter -S . -e "FMT/DP[4] > 34.2" | bcftools filter -S . -e "FMT/DP[5] > 34.5" | bcftools filter -S . -e "FMT/DP[6] > 30.3" | bcftools filter -S . -e "FMT/DP[7] > 33" | bcftools filter -S . -e "FMT/DP[8] > 39" | bcftools filter -S . -e "FMT/DP[9] > 34.2" | bcftools filter -S . -e "FMT/DP[10] > 38.4" | bcftools filter -S . -e "FMT/DP[11] > 32.4" | bcftools filter -S . -e "FMT/DP[12] > 30.9" | bcftools filter -S . -e "FMT/DP[13] > 32.1" | bcftools filter -S . -e "FMT/DP[14] > 30.6" | bcftools filter -S . -e "FMT/DP[15] > 39.9" | bcftools filter -S . -e "FMT/DP[16] > 33.9" | bcftools filter -S . -e "FMT/DP[17] > 33" | bcftools filter -S . -e "FMT/DP[18] > 30.9" | bcftools filter -S . -e "FMT/DP[19] > 35.1" | bcftools filter -S . -e "FMT/DP < 3" | bcftools view -f PASS -e 'ALT="*" | TYPE~"indel" | ref="N"' -O z -o herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.vcf.gz

# Create an index for the VCF file
tabix herring_sentieon_91ind_190521.SV.VF.F2.maxDPtriple.setGT.inv.vcf.gz

```
