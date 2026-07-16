# Apply Neff correction the allele depth (AD)

<!-- MDTOC maxdepth:6 firsth1:2 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Apply Neff correction to allele depth (AD)](#apply-neff-correction-to-allele-depth-ad)   

<!-- /MDTOC -->
# 2020-08-24

## Apply Neff correction to allele depth (AD)
I will use the custom python3 script `apply_Neff_correction_AD_v2.py` (below), which applies the Neff correction to the AD files, it receives the pool sizes (not hard-coded anymore) from the command line, and it handles better missing data (0,0) than previous versions:
```python
#!/usr/bin/python3

"""
Date: May 19, 2020
Created by: Angela P. Fuentes-Pardo, e-mail: apfuentesp@gmail.com, Uppsala University.

Script that applies the Neff correction described by Bergland, A. O., Behrman, E. L., O’Brien, K. R., Schmidt, P. S., & Petrov, D. A. (2014).
Genomic Evidence of Rapid and Stable Adaptive Oscillations over Seasonal Time Scales in Drosophila. PLoS Genetics, 10(11), e1004775. https://doi.org/10.1371/journal.pgen.1004775
The input file should correspond to a text file containing the raw allele depth (AD) extracted from a VCF file generated with GATK SNP caller.

** NEW in this version: It can handle NA and 0,0 sites.

Input file format:

CHROM	POS	REF	ALT	1a_West_Ireland_2016.AD	1b_Southwest_Ireland_2016.AD	2_West-Southwest_Ireland_2017.AD	3_Southern_NorthSea_2016.AD	4_Southern_NorthSea_2017.AD	5a_Northern_Portugal_2016.AD	5b_Southern_Portugal_2016.AD	6a_Northern_Portugal_2017.AD	6b_Southern_Portugal_2017.AD	7_NorthAfrica_Mauritania_2016.AD	8_Northern_SpanishShelf_2016.AD	9_Mediterranean_AlboranSea_2018.AD
scaffold_103_arrow_ctg1	2588	T	C	43,0	51,0	49,0	63,0	67,0	52,0	46,5	35,7	37,0	28,5	50,6	59,0
scaffold_103_arrow_ctg1	2685	T	A	56,0	52,7	75,9	78,0	91,0	53,0	52,0	43,0	37,0	56,0	88,12	70,0
scaffold_103_arrow_ctg1	2741	T	G	33,7	35,8	60,9	62,15	64,6	58,16	44,0	49,8	60,18	42,0	62,22	58,15
scaffold_103_arrow_ctg1	28152	T	G	46,17	54,19	63,16	56,15	61,26	48,15	46,20	49,7	42,23	50,25	51,22	58,18
scaffold_103_arrow_ctg1	28739	T	C	39,10	54,13	61,11	49,6	70,17	50,16	31,11	30,0	43,17	77,28	85,16	64,0
scaffold_103_arrow_ctg1	28776	A	T	38,8	65,0	62,7	53,0	72,3	53,4	45,0	25,5	60,0	74,16	63,0	59,0
scaffold_103_arrow_ctg1	28780	A	T	38,8	65,0	63,7	53,0	74,0	64,0	45,0	26,5	60,0	72,15	63,0	59,0
scaffold_103_arrow_ctg1	56087	C	T	0,0	56,13	53,0	68,0	76,0	51,0	48,0	51,0	42,0	114,0	78,0	86,0
scaffold_103_arrow_ctg1	106404	A	C	NA	NA	NA	56,0	73,0	51,0	48,5	41,0	54,4	86,0	68,0	84,0

Output file format:

CHROM	POS	REF	ALT	1a_West_Ireland_2016.AD	1b_Southwest_Ireland_2016.AD	2_West-Southwest_Ireland_2017.AD	3_Southern_NorthSea_2016.AD	4_Southern_NorthSea_2017.AD	5a_Northern_Portugal_2016.AD	5b_Southern_Portugal_2016.AD	6a_Northern_Portugal_2017.AD	6b_Southern_Portugal_2017.AD	7_NorthAfrica_Mauritania_2016.AD	8_Northern_SpanishShelf_2016.AD	9_Mediterranean_AlboranSea_2018.AD
scaffold_103_arrow_ctg1	2588	T	C	30,0	32,0	35,0	47,0	45,0	36,0	24,3	24,5	26,0	21,4	38,5	36,0
scaffold_103_arrow_ctg1	2685	T	A	35,0	30,5	44,6	55,0	55,0	37,0	27,0	29,0	26,0	37,0	57,8	40,0
scaffold_103_arrow_ctg1	2741	T	G	23,5	23,6	38,6	43,11	42,4	36,10	25,0	30,5	33,10	30,0	42,16	32,9
scaffold_103_arrow_ctg1	28152	T	G	27,11	29,11	38,10	40,11	37,16	32,10	21,10	30,5	24,14	30,15	36,16	32,10
scaffold_103_arrow_ctg1	28739	T	C	25,7	30,8	38,7	37,5	42,11	32,11	17,7	22,0	25,11	39,15	55,11	38,0
scaffold_103_arrow_ctg1	28776	A	T	25,6	37,0	39,5	41,0	46,2	36,3	25,0	18,4	36,0	41,9	47,0	36,0
scaffold_103_arrow_ctg1	28780	A	T	25,6	37,0	39,5	41,0	48,0	42,0	25,0	19,4	36,0	40,9	47,0	36,0
scaffold_103_arrow_ctg1	56087	C	T	0,0	31,8	37,0	50,0	49,0	36,0	26,0	33,0	29,0	56,0	55,0	45,0
scaffold_103_arrow_ctg1	106404	A	C	NA	NA	NA	43,0	47,0	36,0	25,3	28,0	33,3	49,0	50,0	45,0

"""

import argparse
import math

def main():

    parser = argparse.ArgumentParser(description="This script applies the neff (effective number of chromosomes) correction described by Bergland et al. (2014) - https://doi.org/10.1371/journal.pgen.1004775 - based on a CSV file with read counts extracted from a VCf file.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="INFILE", type=str,
                        required=True,
                        help="Input CSV file with read counts.")

    parser.add_argument("-o", "--output", metavar="OUTPUT FILE", type=str,
                        required=True,
                        help="Output file with allele frequencies.")

    parser.add_argument("-s", "--poolsizes", metavar="POOL SIZES", type=str,
                        required=True,
                        help="List of pool sizes. Pool size = 2 x number of individuals in the pool (for diploid species).")


    args = parser.parse_args()
    poolSizes = args.poolsizes
    poolSizes_split = poolSizes.split(",")
    #print(poolSizes_split)

    line_count = 0

    outfile = open(args.output, "w")

    # Read through the CSV file line-by-line.
    with open(args.input, "r") as csv_file:
        for line in csv_file:

            # Remove line terminator/newline character.
            line = line.rstrip("\r\n")

            # Split line on tab characters.
            line_split = line.split("\t")

            # If first line, save header and continues to the next line.
            if line_count == 0:
                headerline = line_split
                #print(headerline)

                print('\t'.join(headerline), file=outfile)
                line_count += 1
                continue

            # For each SNP line, do...
            # Print to the output file the columns CHROM, POS, REF, ALT of the current locus.
            print('\t'.join(line_split[0:4]), end="", file=outfile)

            k = 4 # counter for the pool_size list.

            # Loop through each individual sample (i.e. population) = AD (Allele Depth):
            for i in range(4, len(line_split)):

                if(line_split[i] != "NA"):
                    #print("** k = " + str(k) + ", i = " + str(i) + ", i-k = " + str(i-k))

                    # Extract the REF and ALT allele counts (read depth) of a given population by spliting the string in the comma sign.
                    pop_info = line_split[i].split(",")
                    #print(pop_info)

                    Ref_counts = int(pop_info[0])
                    Alt_counts = int(pop_info[1])

                    # Calculate Nc.
                    total_depth = Ref_counts + Alt_counts
                    # List of strings with pool sizes of 12 HOM populations in this order: 1a_West_Ireland_2016.AD  1b_Southwest_Ireland_2016.AD    2_West-Southwest_Ireland_2017.AD    3_Southern_NorthSea_2016.AD 4_Southern_NorthSea_2017.AD 5a_Northern_Portugal_2016.AD    5b_Southern_Portugal_2016.AD    6a_Northern_Portugal_2017.AD    6b_Southern_Portugal_2017.AD    7_NorthAfrica_Mauritania_2016.AD    8_Northern_SpanishShelf_2016.AD 9_Mediterranean_AlboranSea_2018.AD.
                    #list_pool_size = ['100','90','124','192','140','128','60','96','96','114','192','98']
                    list_pool_size = poolSizes_split
                    #print(list_pool_size)
                    list_pool_size = [int(i) for i in list_pool_size]  # Transform strings into integers.
                    pool_size = list_pool_size[i-k]  # Pool size of the current population.
                    #print(pool_size)

                    # Nc calculation from Kolaczkowski et al. 2011, Feder et al. (2012), Bergland et al. (2014), Machado et al. (2016).
                    Nc = math.floor(((total_depth * pool_size)-1)/(total_depth + pool_size))

                    # If REF is the major allele, apply the correction to it.
                    if(Ref_counts > Alt_counts):
                        Crr_Ref_counts = math.floor((Ref_counts/total_depth) * Nc)
                        Crr_Alt_counts = Nc - Crr_Ref_counts

                    # If ALT is the major allele, apply the correction to it.
                    if(Ref_counts < Alt_counts):
                        Crr_Alt_counts = math.floor((Alt_counts/total_depth) * Nc)
                        Crr_Ref_counts = Nc - Crr_Alt_counts

                    # If REF and ALT have the same coverage, set both depths to the same coverage even if the sum is below Nc.
                    if(Ref_counts == Alt_counts and total_depth != 0):
                        Crr_Ref_counts = math.floor((Ref_counts/total_depth) * Nc)
                        Crr_Alt_counts = Crr_Ref_counts

                    # If REF and ALT are zero, set both depths to zero.
                    if(Ref_counts == Alt_counts and total_depth == 0):
                        Crr_Ref_counts = Ref_counts
                        Crr_Alt_counts = Alt_counts

                    #print("Ref_counts: " + str(Ref_counts) + ", Alt_counts: " + str(Alt_counts) + ", Total depth: " + str(total_depth) +
                    # ", pool_size: " + str(pool_size) + ", Nc: " + str(Nc) + ", Crr_Ref_counts: " + str(Crr_Ref_counts) + ", Crr_Alt_counts: " + str(Crr_Alt_counts))

                    # Print the corrected counts to the output file.
                    print("\t" + str(Crr_Ref_counts) + "," + str(Crr_Alt_counts), end="", file=outfile)
                    #print("\t" + str(Crr_Ref_counts) + "," + str(Crr_Alt_counts) + "\t" + str(Crr_Ref_counts + Crr_Alt_counts), end="", file=outfile)

                    #k += 1

                else:
                    # Print NA to the output file.
                    print("\t" + "NA", end="", file=outfile)

            # Make newline
            print("", file=outfile)

            line_count += 1

    outfile.close()


if __name__ == "__main__":
    main()

```

File listing the input VCF files, `VCF_files_calcAF.txt`:
```bash
cd /proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/06-apply-Neff-calc-AF

ls /proj/snic2020-6-128/private/obtectus_poolseq/04-vcf-files/OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DP*.vcf.gz | awk '{print $9}' > /proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/06-apply-Neff-calc-AF/VCF_files_calcAF.txt

cat /proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/06-apply-Neff-calc-AF/VCF_files_calcAF.txt

/proj/snic2020-6-128/private/obtectus_poolseq/04-vcf-files/OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMean1SD.vcf.gz
/proj/snic2020-6-128/private/obtectus_poolseq/04-vcf-files/OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMode1half.vcf.gz
/proj/snic2020-6-128/private/obtectus_poolseq/04-vcf-files/OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.vcf.gz
```

I used the script `06-1-2-3-apply-Neff-calc-AF.sh` to generate all the AF+ files:
```bash
#!/bin/bash
#SBATCH -A snic2020-15-137
#SBATCH -M snowy
#SBATCH -p core -n 4
#SBATCH -t 20:00:00
#SBATCH -J calcAF
#SBATCH -e calcAF_%J_%A_%a.err
#SBATCH -o calcAF_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-3

# Load required programs.
module load bioinfo-tools
module load GATK/3.8-0
module load python3/3.7.2
module load bcftools/1.10

# Set path to required files and directories as environment variables.
echo This is array job: $SLURM_ARRAY_TASK_ID
VCF_FILE='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/06-apply-Neff-calc-AF/VCF_files_calcAF.txt'
INPUT_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $VCF_FILE)
echo -e This is the VCF target/INPUT_PATH: ${INPUT_PATH}

# Set environment variables.
# Files.
WORK_DIR='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/06-apply-Neff-calc-AF'
REF_PATH='/proj/snic2020-6-128/private/obtectus_poolseq/00-genome/A.obtectus_v2.0.fasta'
REF_FILE=$(basename "${REF_PATH}")
#INPUT_PATH='/home/afuentes/project/afuentes/04-VCF-files/NWAtlantic_herring_15_pools.SNPs.hf.DP20-300.GQ10.mono.miss20.maf0.05.vcf.gz'
INPUT_FILE=$(basename "${INPUT_PATH}")
OUTFILE_PREFIX=$(basename "${INPUT_PATH}" .vcf.gz)
# Programs. (Copied the ones developed for the HOM, notebook date 2020-05-20)
GATK_JAR='/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar'
NEFF_CODE='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/code/utility-code/apply_Neff_correction_AD_v2.py'
AF_CODE='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/code/utility-code/calculate_allelefreq_from_AD_v2.py'
TPMX_CODE='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/code/utility-code/transpose_allelefreq_matrix.py'

# Copy required files to the hard drive of the node in Uppmax.
cp $REF_PATH* $SNIC_TMP  # With *, copy the reference and associated files.
cp ${REF_PATH/.fasta/.dict} $SNIC_TMP  # Copy *.dict file of the reference. Use ${parameter/pattern/string} to replace the first occurrence of a pattern with a given string. Note fasta file can ends with .fa or .fasta, change accordingly
cp $INPUT_PATH* $SNIC_TMP  # With *, the VCF file and its index are copied.

# Go to the SNIC_TMP (working) directory.
cd $SNIC_TMP

# Create directory for temporal files.
if [ -d $SNIC_TMP/tmp ]; then echo "tmp/ exists in SNIC_TMP"; else mkdir $SNIC_TMP/tmp; fi

# 6-1) Extract allele depth (AD) from the VCF file.
java -Djava.io.tmpdir=$SNIC_TMP/tmp -XX:ParallelGCThreads=1 -Dsamjdk.use_async_io=true -Dsamjdk.buffer_size=4194304 -Xmx20G -jar $GATK_JAR \
-T VariantsToTable \
-R $REF_FILE \
-V $INPUT_FILE \
-F CHROM -F POS -F REF -F ALT --genotypeFields AD \
-o $OUTFILE_PREFIX.AD.txt

# 6-2) Apply Neff correction to the raw read counts per allele (AD).
python3 $NEFF_CODE \
-i $OUTFILE_PREFIX.AD.txt \
-o $OUTFILE_PREFIX.AD.Neff.txt \
-s 96,96,96,96,96,96,96,96,96,96,96,96,96,96,96,96 # Pool sizes of the 16 pool samples (2 x number of ind. in a pool = 2 x 48 = 96).

# 6-3) Calculate population-level allele frequencies for each SNP.
python3 $AF_CODE \
-i $OUTFILE_PREFIX.AD.Neff.txt \
-o $OUTFILE_PREFIX.AD.Neff.AF.txt

# Transpose AF matrix using a custom python script.
python3 $TPMX_CODE \
-i $OUTFILE_PREFIX.AD.Neff.AF.txt \
-o $OUTFILE_PREFIX.AD.Neff.AF.tp

# Copy result files to working directory.
cp $SNIC_TMP/$OUTFILE_PREFIX.AD.txt $WORK_DIR
cp $SNIC_TMP/$OUTFILE_PREFIX.AD.Neff.txt $WORK_DIR
cp $SNIC_TMP/$OUTFILE_PREFIX.AD.Neff.AF.txt $WORK_DIR
cp $SNIC_TMP/$OUTFILE_PREFIX.AD.Neff.AF.tp $WORK_DIR

```
Submitted batch job 4350989 on cluster snowy
**Duration: 00-00:56:00**
