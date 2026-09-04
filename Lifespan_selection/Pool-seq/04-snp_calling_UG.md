# SNP calling of the bean beetle dataset with GATK-UnifiedGenotyper

<!-- MDTOC maxdepth:6 firsth1:2 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Perform SNP calling](#perform-snp-calling)   
   - [Input files](#input-files)   
- [Concatenate VCFs](#concatenate-vcfs)   

<!-- /MDTOC -->


# 2021-03-30

## Perform SNP calling
Perform SNP calling for each chromosome separately and for all unplaced scaffolds, for all the 16 pools available.

### Input files
- List of BAM files:
```bash
cp <REDACTED_PATH> <REDACTED_PATH>

cat <REDACTED_PATH>

<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>
<REDACTED_PATH>

```

- List of chromosomes:
```bash
cat <REDACTED_PATH> | grep ">chr" | sort -V | sed -e 's/>//g' > <REDACTED_PATH>

cat <REDACTED_PATH>
chr_1
chr_2
chr_3
chr_4
chr_5
chr_6
chr_7
chr_8
chr_9
chr_10
```
Create script `04-1-SNPcalling-UG-chr.sh`. This script is a job array that scatter the run by chromosome (26 in the case of herring):
```bash
#!/bin/bash
#SBATCH -A <REDACTED_PROJECT_CODE>
#SBATCH -M <REDACTED_HPC_SITE>
#SBATCH -p core -n 4
#SBATCH -t 7-00:00:00
#SBATCH -J SNPcall_UG
#SBATCH -e SNPcall_UG_%J_%A_%a.err
#SBATCH -o SNPcall_UG_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=<REDACTED_EMAIL>
#SBATCH -a 1-10
# Note that the number of jobs to split in the array should be the total number of chromosomes

# Load required software.
module load bioinfo-tools
module load GATK/3.8-0

# Set environment variables.
# Files.
WORK_DIR='<REDACTED_PATH>'
REF_FILE='<REDACTED_PATH>'
OUTFILE_PREFIX='OTT_16_pools'
BAM_LIST='<REDACTED_PATH>'
CHR_LIST='<REDACTED_PATH>'
# Programs.
GATK_JAR='<REDACTED_PATH>'

echo This is array job: $SLURM_ARRAY_TASK_ID
CHR_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $CHR_LIST)
echo -e This is the BAM target: ${CHR_target}

# Go to the directory where the VCF file will be stored.
cd $WORK_DIR

# Run GATK-UnifiedGenotyper to generate a single VCF file for all pool samples.
java -Xmx30G -jar $GATK_JAR \
-T UnifiedGenotyper \
-R $REF_FILE \
-I $BAM_LIST \
-L ${CHR_target} \
-o $OUTFILE_PREFIX.${CHR_target}.vcf

```
Submitted batch job <REDACTED_JOB_ID> on cluster <REDACTED_HPC_SITE>
**Runtime: 00-21:09:00**

*For unplaced scaffolds:*
Same BAM files list.

- List of unplaced scaffolds:
```bash
cat <REDACTED_PATH> | grep ">scaffold" | sort -V | sed -e 's/>//g' > <REDACTED_PATH>

head <REDACTED_PATH>

scaffold_11
scaffold_12
scaffold_13
scaffold_14
scaffold_15
scaffold_16
scaffold_17
scaffold_18
scaffold_19
scaffold_20

tail <REDACTED_PATH>
scaffold_3793
scaffold_3794
scaffold_3795
scaffold_3796
scaffold_3797
scaffold_3798
scaffold_3799
scaffold_3800
scaffold_3801
scaffold_3802
```

Create script `04-1-SNPcalling-UG-scaff.sh`:
```bash
#!/bin/bash
#SBATCH -A <REDACTED_PROJECT_CODE>
#SBATCH -M <REDACTED_HPC_SITE>
#SBATCH -p core -n 4
#SBATCH -t 10-00:00:00
#SBATCH -J SNPcall_UG
#SBATCH -e SNPcall_UG_%J_%A_%a.err
#SBATCH -o SNPcall_UG_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=<REDACTED_EMAIL>

# Load required software.
module load bioinfo-tools
module load GATK/3.8-0

# Set environment variables.
# Files.
WORK_DIR='<REDACTED_PATH>'
REF_FILE='<REDACTED_PATH>'
OUTFILE_PREFIX='OTT_16_pools'
BAM_LIST='<REDACTED_PATH>'
SCAFFOLD_LIST='<REDACTED_PATH>'
# Programs.
GATK_JAR='<REDACTED_PATH>'

# Go to the directory where the VCF file will be stored.
cd $WORK_DIR

# Run GATK-UnifiedGenotyper to generate a single VCF file for all pool samples. Note that GATK requires files suffix to be .list or .interval for the -I flag.
java -Xmx30G -jar $GATK_JAR \
-T UnifiedGenotyper \
-R $REF_FILE \
-I $BAM_LIST \
-L ${SCAFFOLD_LIST} \
-o $OUTFILE_PREFIX.scaff.vcf

```
Submitted batch job <REDACTED_JOB_ID> on cluster <REDACTED_HPC_SITE>
**Runtime: 00-16:04:00**


All the jobs were completed successfully:
```bash
==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>.err <==
INFO  01:58:38,912 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  01:58:38,912 MicroScheduler -   -> 14846447 reads (5.72% of total) failing BadMateFilter
INFO  01:58:38,913 MicroScheduler -   -> 11064408 reads (4.26% of total) failing DuplicateReadFilter
INFO  01:58:38,913 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  01:58:38,913 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  01:58:38,913 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  01:58:38,913 MicroScheduler -   -> 348901 reads (0.13% of total) failing NotPrimaryAlignmentFilter
INFO  01:58:38,913 MicroScheduler -   -> 537143 reads (0.21% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_1.err <==
INFO  07:03:47,050 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  07:03:47,050 MicroScheduler -   -> 30038924 reads (5.89% of total) failing BadMateFilter
INFO  07:03:47,050 MicroScheduler -   -> 15337275 reads (3.01% of total) failing DuplicateReadFilter
INFO  07:03:47,051 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  07:03:47,051 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  07:03:47,051 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  07:03:47,051 MicroScheduler -   -> 2456698 reads (0.48% of total) failing NotPrimaryAlignmentFilter
INFO  07:03:47,051 MicroScheduler -   -> 1565054 reads (0.31% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_2.err <==
INFO  07:01:06,128 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  07:01:06,128 MicroScheduler -   -> 29236292 reads (5.42% of total) failing BadMateFilter
INFO  07:01:06,129 MicroScheduler -   -> 44554272 reads (8.26% of total) failing DuplicateReadFilter
INFO  07:01:06,129 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  07:01:06,129 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  07:01:06,129 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  07:01:06,129 MicroScheduler -   -> 2444212 reads (0.45% of total) failing NotPrimaryAlignmentFilter
INFO  07:01:06,129 MicroScheduler -   -> 1565541 reads (0.29% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_3.err <==
INFO  06:35:59,279 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  06:35:59,280 MicroScheduler -   -> 29027249 reads (5.58% of total) failing BadMateFilter
INFO  06:35:59,280 MicroScheduler -   -> 31340028 reads (6.03% of total) failing DuplicateReadFilter
INFO  06:35:59,280 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  06:35:59,280 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  06:35:59,280 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  06:35:59,280 MicroScheduler -   -> 2734612 reads (0.53% of total) failing NotPrimaryAlignmentFilter
INFO  06:35:59,280 MicroScheduler -   -> 1665058 reads (0.32% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_4.err <==
INFO  06:41:40,529 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  06:41:40,529 MicroScheduler -   -> 27321868 reads (5.20% of total) failing BadMateFilter
INFO  06:41:40,529 MicroScheduler -   -> 26579031 reads (5.06% of total) failing DuplicateReadFilter
INFO  06:41:40,529 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  06:41:40,530 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  06:41:40,530 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  06:41:40,530 MicroScheduler -   -> 2647100 reads (0.50% of total) failing NotPrimaryAlignmentFilter
INFO  06:41:40,530 MicroScheduler -   -> 1573134 reads (0.30% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_5.err <==
INFO  06:54:19,930 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  06:54:19,930 MicroScheduler -   -> 27902754 reads (5.52% of total) failing BadMateFilter
INFO  06:54:19,930 MicroScheduler -   -> 27289354 reads (5.40% of total) failing DuplicateReadFilter
INFO  06:54:19,930 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  06:54:19,930 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  06:54:19,931 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  06:54:19,931 MicroScheduler -   -> 2594775 reads (0.51% of total) failing NotPrimaryAlignmentFilter
INFO  06:54:19,931 MicroScheduler -   -> 1607760 reads (0.32% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_6.err <==
INFO  05:19:43,858 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  05:19:43,858 MicroScheduler -   -> 24465686 reads (4.62% of total) failing BadMateFilter
INFO  05:19:43,858 MicroScheduler -   -> 91795235 reads (17.33% of total) failing DuplicateReadFilter
INFO  05:19:43,858 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  05:19:43,859 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  05:19:43,859 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  05:19:43,859 MicroScheduler -   -> 2218047 reads (0.42% of total) failing NotPrimaryAlignmentFilter
INFO  05:19:43,859 MicroScheduler -   -> 1378487 reads (0.26% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_7.err <==
INFO  04:51:20,211 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  04:51:20,211 MicroScheduler -   -> 20119768 reads (4.78% of total) failing BadMateFilter
INFO  04:51:20,211 MicroScheduler -   -> 22282339 reads (5.29% of total) failing DuplicateReadFilter
INFO  04:51:20,211 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  04:51:20,212 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  04:51:20,212 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  04:51:20,212 MicroScheduler -   -> 1879919 reads (0.45% of total) failing NotPrimaryAlignmentFilter
INFO  04:51:20,212 MicroScheduler -   -> 1231772 reads (0.29% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_8.err <==
INFO  04:03:28,108 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  04:03:28,108 MicroScheduler -   -> 23334714 reads (6.44% of total) failing BadMateFilter
INFO  04:03:28,108 MicroScheduler -   -> 14528142 reads (4.01% of total) failing DuplicateReadFilter
INFO  04:03:28,109 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  04:03:28,109 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  04:03:28,109 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  04:03:28,109 MicroScheduler -   -> 1973683 reads (0.55% of total) failing NotPrimaryAlignmentFilter
INFO  04:03:28,109 MicroScheduler -   -> 1200117 reads (0.33% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_9.err <==
INFO  03:19:51,936 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  03:19:51,936 MicroScheduler -   -> 18274636 reads (5.72% of total) failing BadMateFilter
INFO  03:19:51,936 MicroScheduler -   -> 24600957 reads (7.70% of total) failing DuplicateReadFilter
INFO  03:19:51,936 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  03:19:51,937 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  03:19:51,937 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  03:19:51,937 MicroScheduler -   -> 1369406 reads (0.43% of total) failing NotPrimaryAlignmentFilter
INFO  03:19:51,937 MicroScheduler -   -> 953196 reads (0.30% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------

==> SNPcall_UG_<REDACTED_JOB_ID>_<REDACTED_JOB_ID>_10.err <==
INFO  00:42:28,636 MicroScheduler -   -> 0 reads (0.00% of total) failing BadCigarFilter
INFO  00:42:28,636 MicroScheduler -   -> 9978395 reads (5.80% of total) failing BadMateFilter
INFO  00:42:28,636 MicroScheduler -   -> 2470076 reads (1.44% of total) failing DuplicateReadFilter
INFO  00:42:28,636 MicroScheduler -   -> 0 reads (0.00% of total) failing FailsVendorQualityCheckFilter
INFO  00:42:28,637 MicroScheduler -   -> 0 reads (0.00% of total) failing MalformedReadFilter
INFO  00:42:28,637 MicroScheduler -   -> 0 reads (0.00% of total) failing MappingQualityUnavailableFilter
INFO  00:42:28,637 MicroScheduler -   -> 509686 reads (0.30% of total) failing NotPrimaryAlignmentFilter
INFO  00:42:28,637 MicroScheduler -   -> 499305 reads (0.29% of total) failing UnmappedReadFilter
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------
```


# 2021-03-31

## Concatenate VCFs
As I run UnifiedGenotyper per chromosome, now I will generate a single VCF containing all the variants and chromosomes.

With this commands I created the lines listing all the chr + scaffolds in ascending order:
```bash
cd <REDACTED_PATH>
ls | grep ".vcf$" | awk '{print $9}' | sort -V | awk '{print "-V "$1" \\"}'

# Note that I manually put the scaff file at the end, after the chromosomes
-V OTT_16_pools.chr_1.vcf \
-V OTT_16_pools.chr_2.vcf \
-V OTT_16_pools.chr_3.vcf \
-V OTT_16_pools.chr_4.vcf \
-V OTT_16_pools.chr_5.vcf \
-V OTT_16_pools.chr_6.vcf \
-V OTT_16_pools.chr_7.vcf \
-V OTT_16_pools.chr_8.vcf \
-V OTT_16_pools.chr_9.vcf \
-V OTT_16_pools.chr_10.vcf \
-V OTT_16_pools.scaff.vcf \
```

Create script `04-2-concatenateVCFs-UG.sh`:
```bash
#!/bin/bash
#SBATCH -A <REDACTED_PROJECT_CODE>
#SBATCH -M <REDACTED_HPC_SITE>
#SBATCH -p core -n 4
#SBATCH -t 20:00:00
#SBATCH -J catVar_UG
#SBATCH -e catVar_UG_%J_%A_%a.err
#SBATCH -o catVar_UG_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=<REDACTED_EMAIL>

# Load required software.
module load bioinfo-tools
module load GATK/3.8-0

# Set environment variables.
# Files.
WORK_DIR='<REDACTED_PATH>'
REF_FILE='<REDACTED_PATH>'
OUTFILE_PREFIX='OTT_16_pools'
# Programs.
GATK_JAR='<REDACTED_PATH>'

# Go to the directory where the VCF file will be stored.
cd $WORK_DIR

# Run CatVariants.
java -Xmx20G -cp $GATK_JAR org.broadinstitute.gatk.tools.CatVariants \
-R $REF_FILE \
-V OTT_16_pools.chr_1.vcf \
-V OTT_16_pools.chr_2.vcf \
-V OTT_16_pools.chr_3.vcf \
-V OTT_16_pools.chr_4.vcf \
-V OTT_16_pools.chr_5.vcf \
-V OTT_16_pools.chr_6.vcf \
-V OTT_16_pools.chr_7.vcf \
-V OTT_16_pools.chr_8.vcf \
-V OTT_16_pools.chr_9.vcf \
-V OTT_16_pools.chr_10.vcf \
-V OTT_16_pools.scaff.vcf \
-assumeSorted \
-out $OUTFILE_PREFIX.UG.rawVar.vcf

```
Submitted batch job <REDACTED_JOB_ID> on cluster <REDACTED_HPC_SITE>
**Runtime: 00-00:9:22**

Concatenation of the separate VCF files per chr and scaffold group was completed successfully:
```bash
cat catVar_UG_*.err

You may find GATK files in <REDACTED_PATH>
The hg38bundle is pre-release.
INFO  10:30:15,662 HelpFormatter - -------------------------------------------------------
INFO  10:30:15,665 HelpFormatter - Program Name: org.broadinstitute.gatk.tools.CatVariants
INFO  10:30:15,675 HelpFormatter - Program Args: -R <REDACTED_PATH> -V OTT_16_pools.chr_1.vcf -V OTT_16_pools.chr_2.vcf -V OTT_16_pools.chr_3.vcf -V OTT_16_pools.chr_4.vcf -V OTT_16_pools.chr_5.vcf -V OTT_16_pools.chr_6.vcf -V OTT_16_pools.chr_7.vcf -V OTT_16_pools.chr_8.vcf -V OTT_16_pools.chr_9.vcf -V OTT_16_pools.chr_10.vcf -V OTT_16_pools.scaff.vcf -assumeSorted -out OTT_16_pools.UG.rawVar.vcf
INFO  10:30:15,684 HelpFormatter - Executing as <REDACTED_REMOTE_TARGET> on Linux 3.10.0-1160.15.2.el7.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.8.0_92-b14.
INFO  10:30:15,685 HelpFormatter - Date/Time: 2021/03/31 10:30:15
INFO  10:30:15,685 HelpFormatter - -------------------------------------------------------
INFO  10:30:15,685 HelpFormatter - -------------------------------------------------------
------------------------------------------------------------------------------------------
Done. ------------------------------------------------------------------------------------------
```
