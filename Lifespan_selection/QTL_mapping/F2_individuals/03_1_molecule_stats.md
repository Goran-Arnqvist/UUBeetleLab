# Estimate molecule stats

<!-- TOC depthFrom:2 depthTo:3 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Original code](#original-code)
- [Adapted code to Uppmax](#adapted-code-to-uppmax)
- [Merge output files](#merge-output-files)
- [Plotting in R](#plotting-in-r)
- [Molecule size](#molecule-size)
- [Reads per molecule](#reads-per-molecule)
- [Code shared by Ryan Franckowiak](#code-shared-by-ryan-franckowiak)
- [paste.awk function (paste code below into text file)](#pasteawk-function-paste-code-below-into-text-file)

<!-- /TOC -->

In the haplotagging data analysis pipeline is recommended to examine `reads per molecule` and `molecule size` summary statistics, for a quick test and info on the quality of the data and getting quick stats on linked read, like pcr-duplication (.metrics file) or N50 of reads-prer-DNA-molecule and n50-DNA-molecule-size.

I generated these statistics using a code shared by Marek Kucka, and I created the R scripts for plotting.

## Original code
Code shared by Marek Kucka, called `MarkDupl-dedup-bxsort-linked-reads-histo.sh`:
```bash
file=$1;
echo $file

fbname=$(basename $file .bam)
echo $fbname

dir=$(dirname $file)
echo $dir

#marks duplicates
java -Xmx8g -XX:ParallelGCThreads=10 -jar /fml/chones/local/picard-2.18.25/picard.jar MarkDuplicates \
       I=$file \
       TMP_DIR=./ \
       O=$dir/$fbname.pMarkdup.bam \
       M=$dir/$fbname.pMark.metrics \
CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT \

#de-duplicates the bam
samtools view -@ 10 -h $fbname.pMarkdup.bam -F 1024 -O bam -o $fbname.dedup.bam

#BX-tag sorts the deduplicated bam
samtools sort -@ 10 -t BX $fbname.dedup.bam -T ./$fbname.tmpsort -O BAM -o $fbname.dedup.bxsorted.bam

#make linked-read molecules off of reads with same BX-tag within 50kbp from each other
bed_write.full.pl $fbname.dedup.bxsorted.bam

#removes 00 containing molecules
awk '!/A00C|C00B|B00D|D00/' $fbname.dedup.bxsorted.linked_reads.full.bed > $fbname.dedup.bxsorted.linked_reads.full.no00.bed

#makes histogram of reads per molecule, e.g. how many DNA molecules have 50 reads or 10 reads etc
cut -f10 $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash -s groupby 1 count 1 | sort -k1nr > $fbname.reads.per.mol.log

#writes histogram of molecule size in 1kb bins, sorted from longest to shortest molecule
awk '{ print $3-$2 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash bin:1000 1 | datamash -s groupby 1 count 1 | sort -k1nr > $fbname.1kb-bin.molecule.histogram

reads=$(cut -f10 $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')
size=$(awk '{ print $3-$2 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')

#calculates N50_reads_per_molecule
cut -f10 $fbname.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr | awk -v val="$reads" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' |tail -1 > $fbname.N50_reads_per_mol.info

#calculates N50_molecule_length
awk '{ print $3-$2 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr |  awk -v val="$size" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' | tail -1 > $fbname.N50_molecule_size.info

#writes all size-sorted DNA molecules; showing Mol-size, BX-BC, number of reads/mol, and chr location
awk '{ print $3-$2"\t"$4"\t"$10"\t"$1"\t"$2"\t"$3 }' $fbname.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr > $fbname_Mol-size_BC_ReadsPerMol_location.log

rm $fbname.dedup.bam
rm $fbname.dedup.bxsorted.bam
```

## Adapted code to Uppmax

Script `03_2_molecule_stats.sh`:

- Make sure the `pool_bam.list` file lists the files you want:
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

wc -l /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/bam_sample.list
1342
```
**NOTE:** The original script shared by Marek uses the perl script called `bed_write.full.pl`. In the github page there is a script called `bed_write.pl`. Marek clarified they are the same, just the name is different, so I used the one available in the repository.

Also, I asked people in Uppmax to install the program `datamash`, which is required for Marek's script.

Given that Uppmax only accepts 1000 jobs at a time, I split the files into two, one goes from 1-1000, and the other one from 1001 to 1342.

Script `03-1-molecule_stats_1-1000.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 2
#SBATCH -t 0-10:00:00
#SBATCH -J molecStats
#SBATCH -e molecStats_%J_%A_%a.err
#SBATCH -o molecStats_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-1000

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample'
#REF_FASTA='/proj/snic2020-6-128/private/a_obtectus_poolseq/00-genome/A.obtectus_v2.0.fasta'
FILE_LIST=$OUTPUT_DIR/bam_sample.list
PICARD_JAR='/sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar'
PERL_SCRIPT_PATH='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/code'

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: $FILE_PATH

# Get file name
FBNAME=$(basename $FILE_PATH .bam)
# Get directory path
INPUT_DIR=$(dirname $FILE_PATH)
# Number of cores to use
N_CORES='2'

# Load packages in Uppmax
module load bioinfo-tools samtools/1.14 picard/2.23.4 perl/5.26.2 datamash/1.8

#### Continuation from the 03-1-read_mapping.sh script ####

# create directory for temporal files
if [ -d $OUTPUT_DIR/Molecule-stats ]; then echo "OUTPUT_DIR/Molecule-stats/ exists in OUTPUT_DIR"; else mkdir $OUTPUT_DIR/Molecule-stats; fi

# Go to the working directory
cd $OUTPUT_DIR/Molecule-stats

# Mark duplicated reads using BX-aware options
java -Xmx12g -XX:ParallelGCThreads=$N_CORES -jar $PICARD_JAR MarkDuplicates \
I=$FILE_PATH \
O=$FBNAME.pMarkdup.bam \
M=$FBNAME.pMarkdup.metrics \
CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT
#CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT && rm $FBNAME.sorted.bam

# create directory for temporal files
if [ -d $FBNAME.tmpsort ]; then echo "tmp/ exists in OUTPUT_DIR"; else mkdir $FBNAME.tmpsort; fi

# De-duplicates the bam (-F 1024 = read is PCR or optical duplicate, 0x400)
samtools view -@ $N_CORES -h $FBNAME.pMarkdup.bam -F 1024 -O bam -o $FBNAME.dedup.bam
#samtools view -@ $N_CORES -h $FILE_PATH -F 1024 -O bam -o $FBNAME.dedup.bam

# BX-tag sorts the deduplicated bam
samtools sort -@ $N_CORES -t BX $FBNAME.dedup.bam -T ./$FBNAME.tmpsort -O bam -o $FBNAME.dedup.bxsorted.bam

# Make linked-read molecules off of reads with same BX-tag within 50 kbp from each other
$PERL_SCRIPT_PATH/bed_write.pl $FBNAME.dedup.bxsorted.bam
#$SCRIPT_PATH/bed_write.full.pl $FBNAME.dedup.bxsorted.bam

# Remove 00 containing molecules
awk '!/A00C|C00B|B00D|D00/' $FBNAME.dedup.bxsorted.linked_reads.full.bed > $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed

# Make histogram of reads per molecule, e.g. how many DNA molecules have 50 reads or 10 reads etc
cut -f10 $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash -s groupby 1 count 1 | sort -k1nr > $FBNAME.reads.per.mol.log

# Write histogram of molecule size in 1 kb bins, sorted from longest to shortest molecule
awk '{ print $3-$2 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash bin:1000 1 | datamash -s groupby 1 count 1 | sort -k1nr > $FBNAME.1kb-bin.molecule.histogram

reads=$(cut -f10 $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')
size=$(awk '{ print $3-$2 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')

# Calculate N50_reads_per_molecule
cut -f10 $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr | awk -v val="$reads" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' | tail -1 > $FBNAME.N50_reads_per_mol.info

# Calculate N50_molecule_length
awk '{ print $3-$2 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr |  awk -v val="$size" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' | tail -1 > $FBNAME.N50_molecule_size.info

# Write all size-sorted DNA molecules; showing Mol-size, BX-BC, number of reads/mol, and chr location
awk '{ print $3-$2"\t"$4"\t"$10"\t"$1"\t"$2"\t"$3 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr > $FBNAME\_Mol-size_BC_ReadsPerMol_location.log

# Clean up
#rm $FBNAME.pMarkdup.bam
#rm $FBNAME.dedup.bam
#rm $FBNAME.dedup.bxsorted.bam
rm -R $FBNAME.tmpsort

```
Submitted batch job 30845657
**Runtime: 00-00:13:00** 3 min per file on average??

Make sure there is a file listing the bam files from 1001 to 1342:
```bash
head /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/bam_sample_1001-1342.list
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA6_24.D59_N701.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC14_34.D59_N707.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/AB10_33.D39_N702.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC2_15.D64_N712.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC15_31.D70_N705.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/AB1_30.D86_N702.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/BA2_3.D47_N710.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/AB6_30.D12_N702.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC1_31.D95_N701.bam
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/DC16_31.D60_N706.bam

wc -l /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/bam_sample_1001-1342.list
342
```

Script `03-1-molecule_stats_1001-1342.sh`:

```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -p core -n 2
#SBATCH -t 0-10:00:00
#SBATCH -J molecStats
#SBATCH -e molecStats_%J_%A_%a.err
#SBATCH -o molecStats_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-342

# Set environment variables
OUTPUT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample'
#REF_FASTA='/proj/snic2020-6-128/private/a_obtectus_poolseq/00-genome/A.obtectus_v2.0.fasta'
FILE_LIST=$OUTPUT_DIR/bam_sample_1001-1342.list
PICARD_JAR='/sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar'
PERL_SCRIPT_PATH='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/code'

# For a given job in the array, set the correspondent sample files
FILE_PATH=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE_LIST)

# Print current file info to stdout for future reference
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the FILE: $FILE_PATH

# Get file name
FBNAME=$(basename $FILE_PATH .bam)
# Get directory path
INPUT_DIR=$(dirname $FILE_PATH)
# Number of cores to use
N_CORES='2'

# Load packages in Uppmax
module load bioinfo-tools samtools/1.14 picard/2.23.4 perl/5.26.2 datamash/1.8

#### Continuation from the 03-1-read_mapping.sh script ####

# create directory for temporal files
if [ -d $OUTPUT_DIR/Molecule-stats ]; then echo "OUTPUT_DIR/Molecule-stats/ exists in OUTPUT_DIR"; else mkdir $OUTPUT_DIR/Molecule-stats; fi

# Go to the working directory
cd $OUTPUT_DIR/Molecule-stats

# Mark duplicated reads using BX-aware options
java -Xmx12g -XX:ParallelGCThreads=$N_CORES -jar $PICARD_JAR MarkDuplicates \
I=$FILE_PATH \
O=$FBNAME.pMarkdup.bam \
M=$FBNAME.pMarkdup.metrics \
CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT
#CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT && rm $FBNAME.sorted.bam

# create directory for temporal files
if [ -d $FBNAME.tmpsort ]; then echo "tmp/ exists in OUTPUT_DIR"; else mkdir $FBNAME.tmpsort; fi

# De-duplicates the bam (-F 1024 = read is PCR or optical duplicate, 0x400)
samtools view -@ $N_CORES -h $FBNAME.pMarkdup.bam -F 1024 -O bam -o $FBNAME.dedup.bam
#samtools view -@ $N_CORES -h $FILE_PATH -F 1024 -O bam -o $FBNAME.dedup.bam

# BX-tag sorts the deduplicated bam
samtools sort -@ $N_CORES -t BX $FBNAME.dedup.bam -T ./$FBNAME.tmpsort -O bam -o $FBNAME.dedup.bxsorted.bam

# Make linked-read molecules off of reads with same BX-tag within 50 kbp from each other
$PERL_SCRIPT_PATH/bed_write.pl $FBNAME.dedup.bxsorted.bam
#$SCRIPT_PATH/bed_write.full.pl $FBNAME.dedup.bxsorted.bam

# Remove 00 containing molecules
awk '!/A00C|C00B|B00D|D00/' $FBNAME.dedup.bxsorted.linked_reads.full.bed > $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed

# Make histogram of reads per molecule, e.g. how many DNA molecules have 50 reads or 10 reads etc
cut -f10 $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash -s groupby 1 count 1 | sort -k1nr > $FBNAME.reads.per.mol.log

# Write histogram of molecule size in 1 kb bins, sorted from longest to shortest molecule
awk '{ print $3-$2 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash bin:1000 1 | datamash -s groupby 1 count 1 | sort -k1nr > $FBNAME.1kb-bin.molecule.histogram

reads=$(cut -f10 $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')
size=$(awk '{ print $3-$2 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | datamash sum 1 | awk '{print $1/2}')

# Calculate N50_reads_per_molecule
cut -f10 $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr | awk -v val="$reads" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' | tail -1 > $FBNAME.N50_reads_per_mol.info

# Calculate N50_molecule_length
awk '{ print $3-$2 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr |  awk -v val="$size" '($1+prev)>val{exit} ($1+prev)<=val; {prev+=$1}' | tail -1 > $FBNAME.N50_molecule_size.info

# Write all size-sorted DNA molecules; showing Mol-size, BX-BC, number of reads/mol, and chr location
awk '{ print $3-$2"\t"$4"\t"$10"\t"$1"\t"$2"\t"$3 }' $FBNAME.dedup.bxsorted.linked_reads.full.no00.bed | sort -k1nr > $FBNAME\_Mol-size_BC_ReadsPerMol_location.log

# Clean up
#rm $FBNAME.pMarkdup.bam
#rm $FBNAME.dedup.bam
#rm $FBNAME.dedup.bxsorted.bam
rm -R $FBNAME.tmpsort

```
Submitted batch job 30845699
**Runtime: 00-00:10:00** 3 min per file on average??

These output files were generated for each sample:
```bash
/Molecule-stats$ head DC9_8.D68_N704*

==> DC9_8.D68_N704.dedup.bxsorted.linked_reads.full.bed <==
chr_1	1506787	1507168	A00C00B01D68-N704	60	-	1506787	1507168	0,0,0	2	117,150	0,231
chr_1	5091829	5091980	A00C00B01D68-N704	0	+	5091829	5091980	0,0,0	1	151	0
chr_1	17060148	17060296	A00C00B01D68-N704	0	-	17060148	17060296	0,0,0	1	148	0
chr_1	39194650	39194888	A00C00B01D68-N704	0	-	39194650	39194888	0,0,0	2	117,151	0,87
chr_1	43302732	43303111	A00C00B01D68-N704	57.5	+	43302732	43303111	0,0,0	2	151,89	0,290
chr_1	67556445	67556596	A00C00B01D68-N704	0	+	67556445	67556596	0,0,0	1	151	0
chr_1	86439832	86439983	A00C00B01D68-N704	25.5	-	86439832	86439983	0,0,0	2	151,151	0,0
chr_1	95236028	95236200	A00C00B01D68-N704	0	-	95236028	95236200	0,0,0	2	116,151	0,21
chr_1	118875650	118875744	A00C00B01D68-N704	0	+	118875650	118875744	0,0,0	2	95,94	0,0
chr_2	5611966	5612317	A00C00B01D68-N704	0	-	5611966	5612317	0,0,0	2	116,150	0,201

==> DC9_8.D68_N704.dedup.bxsorted.linked_reads.full.no00.bed <==
chr_7	62566904	62567358	A01C01B04D68-N704	60	+	62566904	62567358	0,3,0	2	151,119	0,335
chr_4	79043621	79043793	A01C01B07D68-N704	60	+	79043621	79043793	0,3,0	2	151,119	0,53
chr_5	32984258	32984683	A01C01B11D68-N704	26	-	32984258	32984683	0,3,0	2	117,151	0,274
chr_2	97652868	97653107	A01C01B15D68-N704	60	+	97652868	97653107	0,3,0	2	151,119	0,120
chr_5	39569290	39569669	A01C01B17D68-N704	59.5	+	39569290	39569669	0,3,0	2	151,118	0,261
chr_4	41171335	41171447	A01C01B19D68-N704	60	+	41171335	41171447	0,3,0	2	112,112	0,0
chr_8	60875302	60875484	A01C01B21D68-N704	60	-	60875302	60875484	0,3,0	2	119,150	0,32
chr_2	80092250	80092518	A01C01B23D68-N704	53	-	80092250	80092518	0,3,0	2	119,151	0,117
chr_1	17156268	17156535	A01C01B30D68-N704	0	-	17156268	17156535	0,3,0	2	119,151	0,116
chr_2	6879964	6880231	A01C01B30D68-N704	0	-	6879964	6880231	0,3,0	2	119,151	0,116

==> DC9_8.D68_N704.reads.per.mol.log <==
54	1
46	1
45	1
44	2
42	2
40	1
39	1
38	4
37	2
36	10

==> DC9_8.D68_N704.1kb-bin.molecule.histogram <==
217000	1
191000	1
179000	1
176000	1
169000	1
152000	1
141000	1
139000	3
138000	1
137000	1

==> DC9_8.D68_N704.N50_molecule_size.info <==
16506

==> DC9_8.D68_N704.N50_reads_per_mol.info <==
4

==> DC9_8.D68_N704_Mol-size_BC_ReadsPerMol_location.log <==
217314	A66C46B59D68-N704	32	chr_7	63642853	63860167
191157	A46C51B46D68-N704	36	chr_7	94105212	94296369
179730	A87C55B08D68-N704	22	scaffold_20	228250	407980
176106	A64C39B69D68-N704	18	chr_4	93269380	93445486
169431	A43C50B22D68-N704	18	chr_7	94031877	94201308
152101	A29C89B96D68-N704	18	chr_5	50426859	50578960
141994	A74C91B14D68-N704	19	chr_7	94033864	94175858
139887	A37C59B30D68-N704	10	scaffold_20	79348	219235
139664	A29C96B35D68-N704	12	chr_7	94053177	94192841
139351	A71C03B09D68-N704	18	chr_5	50401532	50540883

```

## Merge output files

- To merge the N50 one-liner files generated for each individual, I used:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/Molecule-stats

cat *.N50_molecule_size.info > summary_N50_molecule_size.info
cat *.N50_reads_per_mol.info > summary_N50_reads_per_mol.info
```

- Merge the `*reads.per.mol.log` and `*1kb-bin.molecule.histogram` files:
```bash
awk '{print FILENAME"\t"$0}' *.reads.per.mol.log | sed 's/.reads.per.mol.log//g' > merged_reads.per.mol.log.txt

awk '{print FILENAME"\t"$0}' *.1kb-bin.molecule.histogram | sed 's/.1kb-bin.molecule.histogram//g' > merged_1kb-bin.molecule.histogram.txt
```


## Plotting in R

Script `03-molecule-stats.R`
```R
# Molecule stats haplotagging data

# Set environmental variables
work_dir <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping"
rpm_file <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping/Molecule-stats/merged_reads.per.mol.log.txt"
#rpm_file <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping/Molecule-stats/DC9_8.D68_N704.reads.per.mol.log"
hist_file <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping/Molecule-stats/merged_1kb-bin.molecule.histogram.txt"
#hist_file <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping/Molecule-stats/DC9_8.D68_N704.1kb-bin.molecule.histogram"
#log_file <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping/Molecule-stats/DC9_8.D68_N704_Mol-size_BC_ReadsPerMol_location.log"
N50_molsize_file <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping/Molecule-stats/summary_N50_molecule_size.info"
N50_readsmol_file <- "~/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/analysis/03_read_mapping/Molecule-stats/summary_N50_reads_per_mol.info"

# Set working directory
setwd(work_dir)

# Load packages
library(data.table)
library(ggplot2)
require(scales)

# Load data
rpm_df <- fread(rpm_file, data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)
hist_df <- fread(hist_file, data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)
#log_df <- fread(log_file, data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)
N50_molsize <- fread(N50_molsize_file, data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)
N50_readsmol <- fread(N50_readsmol_file, data.table = FALSE, header = FALSE, stringsAsFactors = FALSE)


# histogram of reads per molecule, e.g. how many DNA molecules have 50 reads or 10 reads, etc
head(rpm_df)

p1 <- ggplot(rpm_df) +
  geom_line(aes(x = V2, y = V3, color = V1), alpha = 0.5) + #geom_histogram
  scale_y_continuous(trans = 'log2') +
  #xlab("reads per molecule") + ylab("count") +
  xlab("reads per molecule") + ylab("log2(count)") +
  theme_bw() +
  theme(text = element_text(size = 18),
        legend.position = "none")

# Save plot
#pdf("plot_n1342_reads_per_molecule.pdf")
pdf("plot_n1342_reads_per_molecule_log2.pdf")
p1
dev.off()

# histogram of molecule size in 1 kb bins
head(hist_df)
require(scales)

p2 <- ggplot(hist_df) +
  geom_line(aes(x = V2, y = V3, color = V1), alpha = 0.5) + #geom_histogram
  scale_y_continuous(trans = 'log2') +
  #xlab("molecule size") + ylab("count") +
  xlab("molecule size") + ylab("log2(count)") +
  theme_bw() +
  theme(text = element_text(size = 18),
        legend.position = "none") +
  scale_x_continuous(labels = comma)

# Save plot
#pdf("plot_n1342_molecule_size.pdf")
pdf("plot_n1342_molecule_size_log2.pdf")
p2
dev.off()


# N50 stats ---------------------------------------------------------------

## Molecule size
sum(is.na(N50_molsize))
N50_molsize_median <- median(N50_molsize$V1)
N50_molsize_median
#[1] 11899

# Plot
p3 <- ggplot(N50_molsize) +
  #geom_histogram(aes(x = V1)) +
  geom_density(aes(x = V1), color = "purple", size = 1) +
  xlab("N50 - molecule size (bp)") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  scale_x_continuous(labels = comma)

# Save plot
pdf("plot_N50_molecule_size.pdf")
p3
dev.off()

## Reads per molecule
sum(is.na(N50_readsmol))
N50_readsmol_median <- median(N50_readsmol$V1)
N50_readsmol_median
#[1] 2

# Plot
p4 <- ggplot(N50_readsmol) +
  #geom_histogram(aes(x = V1)) +
  geom_density(aes(x = V1), color = "blue", size = 1) +
  xlab("N50 - reads per molecule") +
  theme_bw() +
  theme(text = element_text(size = 18)) +
  scale_x_continuous(labels = comma)

# Save plot
pdf("plot_N50_reads_per_molecule.pdf")
p4
dev.off()

```

## Code shared by Ryan Franckowiak

Ryan Franckowiak, Cornell University, shared the code below, which I did not use at the end but was helpful to figure out which input files I should use to generate the plots.

The script is called `Molecule_Scripts.rtf`:
```bash
######################################
# Step 00: Sample file (sample_list.tsv)
######################################
Sample file format:

ShadHap1        L002    C01     SC001   SC      Illumina        NovaSeq6000
ShadHap1        L002    C02     SC002   SC      Illumina        NovaSeq6000
ShadHap1        L002    C03     SC003   SC      Illumina        NovaSeq6000
ShadHap1        L002    C04     SC004   SC      Illumina        NovaSeq6000
ShadHap1        L002    C05     SC005   SC      Illumina        NovaSeq6000
ShadHap1        L002    C06     SC006   SC      Illumina        NovaSeq6000
ShadHap1        L002    C07     SC007   SC      Illumina        NovaSeq6000
ShadHap1        L002    C08     SC008   SC      Illumina        NovaSeq6000
ShadHap1        L002    C09     CO065   CO      Illumina        NovaSeq6000
ShadHap1        L002    C10     CO067   CO      Illumina        NovaSeq6000

######################################
# Step 01: BX sorting
######################################
#!/bin/bash

# Software path
export PATH=/programs/samtools-1.11/bin:$PATH

# Global variables
INFO_FILES="02-info_files"
GENOME="03-genome"
INPUTFOLDER="12-rm_duplicates"
OUTFOLDER="13-linked_reads"
LOGFOLDER="99-log_files"

 N_CORE_MAX=36
 COUNT=0
 THREADS=1

while read i; do
  LIBRARY=$(echo $i | cut -d ' ' -f1)
  LANE=$(echo $i | cut -d ' ' -f2)
  BARCODE=$(echo $i | cut -d ' ' -f3)
  BAMBASE=$LIBRARY'_'$LANE'_'$BARCODE
  # Run samtools
samtools sort -@ $THREADS -t BX ${INPUTFOLDER}/${BAMBASE}'.dedup.bam' -O BAM -o   ${OUTFOLDER}/${BAMBASE}'.dedup.bxsorted.bam' --reference ${GENOME}'/fAloSap1.pri.cur.20210421.fasta' &
         COUNT=$(( COUNT + 1 ))
       if [ $COUNT == $N_CORE_MAX ]; then
          wait
    COUNT=0
  Fi
done<${INFO_FILES}'/sample_list.tsv'

######################################
# Step 02:
######################################
#!/bin/bash

#Global variables
INFO_FILES="02-info_files"
BAMDIR="13-linked_reads"
LOGFOLDER="99-log_files"

N_CORE_MAX=24
COUNT=0
THREADS=1

while read i; do
  LIBRARY=$(echo $i | cut -d ' ' -f1)
  LANE=$(echo $i | cut -d ' ' -f2)
  BARCODE=$(echo $i | cut -d ' ' -f3)
  BAMBASE=$LIBRARY'_'$LANE'_'$BARCODE

 perl ${INFO_FILES}'/bed_write.full.pl' ${BAMDIR}'/'${BAMBASE}'.dedup.bxsorted.bam'  &
         COUNT=$(( COUNT + 1 ))
       if [ $COUNT == $N_CORE_MAX ]; then
          wait
    COUNT=0
  fi
 done<${INFO_FILES}'/sample_list.tsv'

######################################
# Step 03: Removing reads w/ double zero '00' barcodes
######################################
#!/bin/bash

# Global variables
INFO_FILES="02-info_files"
BAMDIR="13-linked_reads"
LOGFOLDER="99-log_files"

# Number of cores
N_CORE_MAX=36
COUNT=0
THREADS=1

while read i; do
  LIBRARY=$(echo $i | cut -d ' ' -f1)
  LANE=$(echo $i | cut -d ' ' -f2)
  BARCODE=$(echo $i | cut -d ' ' -f3)
  BAMBASE=$LIBRARY'_'$LANE'_'$BARCODE

  awk '!($4~/A00|B00|C00|D00/)' ${BAMDIR}'/'${BAMBASE}'.dedup.bxsorted.linked_reads.full.bed' > ${BAMDIR}'/'${BAMBASE}'.dedup.bxsorted.linked_reads.full.no00.bed' &
         COUNT=$(( COUNT + 1 ))
       if [ $COUNT == $N_CORE_MAX ]; then
          wait
    COUNT=0
  fi
 done<${INFO_FILES}'/sample_list.tsv'

######################################
# Step 04:
######################################
#1/bin/bash

# Global variables
INFODIR="02-info_files"
BEDDIR="13-linked_reads"
OUTDIR="13-reads_per_molecule"
LOGDIR="99-log_files"

# Number of cores
N_CORE_MAX=36
COUNT=0
THREADS=1

while read i; do
  LIBRARY=$(echo $i | cut -d ' ' -f1)
  LANE=$(echo $i | cut -d ' ' -f2)
  BARCODE=$(echo $i | cut -d ' ' -f3)
  BEDBASE=$LIBRARY'_'$LANE'_'$BARCODE

  cut -f10 ${BEDDIR}'/'${BEDBASE}'.dedup.bxsorted.linked_reads.full.no00.bed' |
        datamash -s groupby 1 count 1 | sort -k 1 -n > ${BEDDIR}'/'${OUTDIR}'/'${BEDBASE}'.reads.per.mol.log' &
         COUNT=$(( COUNT + 1 ))
       if [ $COUNT == $N_CORE_MAX ]; then
          wait
    COUNT=0
  fi
 done<${INFODIR}'/sample_list.tsv'


######################################
# Step 05:
######################################
#1/bin/bash

# Merge reads_per_molecule files

awk -f ../../02-info_files/paste.awk `cat 13-summary/list_reads_per_molecule.logs` >> 13-summary/paste_reads_per_molecule.log

## paste.awk function (paste code below into text file)

#!/usr/bin/awk -f
BEGIN {
    RS = "(\r\n|\n\r|\r|\n)"
    FS = " *\t *"
    SUBSEP = ":"
}
FNR==1 {
    ++file
}
NF>=2 {
    if ($1 in keynum)
        key = keynum[$1]
    else {
        key = ++keys
        keynum[$1] = key
        keystr[key] = $1
    }
    printf "key = %s, file = %s, value = %s\n", key, file, $2 >/dev/stderr
    value[key,file] = $2
}
END {
    files = file
    for (key = 1; key <= keys; key++) {
        printf "%s", keystr[key]
        for (file = 1; file <= files; file++)
            printf "\t%s", value[key,file]
        printf "\n"
    }
}

######################################
# Step 06: molecule sizes
######################################
#!/bin/bash

# Global variables
INFODIR="02-info_files"
BEDDIR="13-linked_reads"
OUTDIR="13-molecule_size"
LOGDIR="99-log_files"

# Number of cores
N_CORE_MAX=24
COUNT=0
THREADS=1

while read i; do
  LIBRARY=$(echo $i | cut -d ' ' -f1)
  LANE=$(echo $i | cut -d ' ' -f2)
  BARCODE=$(echo $i | cut -d ' ' -f3)
  BEDBASE=$LIBRARY'_'$LANE'_'$BARCODE
  # Calculate molecule size
  awk '{ print $1"\t"$2"\t"$3"\t"$3-$2"\t"$4"\t"$10 }' ${BEDDIR}'/'${BEDBASE}'.dedup.bxsorted.linked_reads.full.no00.bed' |
        sort -k 4 -n > ${BEDDIR}'/'${OUTDIR}'/'${BEDBASE}'_molecule.sizes.sorted.all.bed' &
         COUNT=$(( COUNT + 1 ))
       if [ $COUNT == $N_CORE_MAX ]; then
          wait
    COUNT=0
  fi
 done<${INFODIR}'/sample_list.tsv'

######################################
# Step 07: molecule histogram
######################################
#!/bin/bash

# Global variables
INFODIR="02-info_files"
BEDDIR="13-linked_reads"
OUTDIR="13-molecule_size"
LOGDIR="99-log_files"

# Number of cores
N_CORE_MAX=12
COUNT=0
THREADS=1

while read i; do
  LIBRARY=$(echo $i | cut -d ' ' -f1)
  LANE=$(echo $i | cut -d ' ' -f2)
  BARCODE=$(echo $i | cut -d ' ' -f3)
  BEDBASE=$LIBRARY'_'$LANE'_'$BARCODE
  # Generate histogram bins

  cut -f4 ${BEDDIR}'/'${OUTDIR}'/'${BEDBASE}'_molecule.sizes.sorted.all.bed' | datamash bin:1000 1 | datamash -s groupby 1 count 1 |
        sort -k 1 -n > ${BEDDIR}'/'${OUTDIR}'/13-summary/'${BEDBASE}'_1kb-bin.molecule.histogram' &
         COUNT=$(( COUNT + 1 ))
       if [ $COUNT == $N_CORE_MAX ]; then
          wait
    COUNT=0
  fi
 done<${INFODIR}'/sample_list.tsv'

######################################
# Step 08: barcode (BX-tag)
######################################
#!/bin/bash

# Set timestamp
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Copy script
SCRIPT=$0
NAME=$(basename $0)
SCRIPT_FOLDER="98-scripts_log"

cp $SCRIPT $SCRIPT_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
INFODIR="02-info_files"
BEDDIR="13-linked_reads"
DATADIR="13-molecule_size"
OUTDIR="13-summary"

ls 13-linked_reads/13-molecule_size/*_1kb-bin.molecule.histogram > 13-linked_reads/13-molecule_size/13-summary/list_molecule_size.logs
awk -f 02-info_files/paste.awk `cat 13-linked_reads/13-molecule_size/13-summary/list_molecule_size.logs` >> 13-linked_reads/13-molecule_size/13-summary/paste_molecule_size.log


######################################
# Step 09: Plot in R
######################################
#Load the libraries
library(data.table)
library(formattable)
library(tidyverse)
library(ggplot2)

# Set working directory (i.e., path relative to script)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dirpath=dirname(dirname(getwd()))

# Set directories for testing
datadir <- file.path(dirpath, "13-linked_reads/13-reads_per_molecule/13-summary/")

#######################
# Primary assembly report
reads_per_mol <- fread(file=paste0(datadir, "paste_reads_per_molecule.log"),
                         data.table=FALSE, header = FALSE, stringsAsFactors = FALSE)

sample_list <- fread(file=paste0(datadir, "list_reads_per_molecule.logs"),
                       data.table=FALSE, header = FALSE, stringsAsFactors = FALSE)

# Remove string from sample names
sample_list$V1 <- gsub('.reads.per.mol.log', '', sample_list$V1)

#
sample_list<-rbind(data.frame(V1 = 'reads'), sample_list)

#
colnames(reads_per_mol)<- sample_list$V1

#
SH1<-reads_per_mol %>% dplyr::select(contains("ShadHap1")) %>%
  cbind(reads = reads_per_mol$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

SH2<-reads_per_mol %>% dplyr::select(contains("ShadHap2")) %>%
  cbind(reads = reads_per_mol$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

SH3<-reads_per_mol %>% dplyr::select(contains("ShadHap3")) %>%
  cbind(reads = reads_per_mol$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

SH4<-reads_per_mol %>% dplyr::select(contains("ShadHap4")) %>%
  cbind(reads = reads_per_mol$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

#
shadhap1 <-SH1 %>% mutate(counts=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,counts)
shadhap2 <-SH2 %>% mutate(counts=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,counts)
shadhap3 <-SH3 %>% mutate(counts=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,counts)
shadhap4 <-SH4 %>% mutate(counts=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,counts)

# Join molecule per read files
library(dplyr)
shap_join <- shadhap1 %>% full_join(shadhap2,by ="reads") %>%
             full_join(shadhap3,by ="reads") %>%
             full_join(shadhap4,by ="reads")

# Rename columns
shap_join<- shap_join %>% dplyr::rename(shadhap1 = counts.x, shadhap2= counts.y,
                                 shadhap3 = counts.x.x, shadhap4= counts.y.y )

# Data preparation
library("tidyverse")
shap_wrgl <- shap_join %>%
  dplyr::select(reads, shadhap1, shadhap2, shadhap3, shadhap4) %>%
  gather(key = "variable", value = "value", -reads)
head(shap_wrgl)

# Visualization

p<- ggplot(shap_wrgl, aes(x = reads, y = value)) +
  geom_line(aes(color = variable),size=0.4) +
  scale_color_manual(values=c("darkred", "steelblue","green","orange")) +
  scale_y_continuous(trans='log2') +
  labs(x ="reads-per-molecule", y = "count") +
  theme_bw()

# Save image to directory
p + ggsave(filename = paste0(datadir, "./reads-per-molecule.png"))


#################################

# Set directories for testing
datadir <- file.path(dirpath, "13-linked_reads/13-molecule_size/13-summary/")

# Primary assembly report
molecule_size <- fread(file=paste0(datadir, "paste_molecule_size.log"),
                       data.table=FALSE, header = FALSE, stringsAsFactors = FALSE)

file_list <- fread(file=paste0(datadir, "list_molecule_size.logs"),
                     data.table=FALSE, header = FALSE, stringsAsFactors = FALSE)

# Remove string from samle names
file_list$V1 <- gsub('_1kb-bin.molecule.histogram', '', file_list$V1)

#
file_list<-rbind(data.frame(V1 = 'reads'), file_list)

#
colnames(molecule_size)<- file_list$V1

#
library(tidyselect)

SH1<-molecule_size %>% dplyr::select(contains("ShadHap1")) %>%
  cbind(reads = molecule_size$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

SH2<-molecule_size %>% dplyr::select(contains("ShadHap2")) %>%
  cbind(reads = molecule_size$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

SH3<-molecule_size %>% dplyr::select(contains("ShadHap3")) %>%
  cbind(reads = molecule_size$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

SH4<-molecule_size %>% dplyr::select(contains("ShadHap4")) %>%
  cbind(reads = molecule_size$reads) %>%
  dplyr::select(reads, everything()) %>% replace(is.na(.), 0)

#
shadhap1 <-SH1 %>% mutate(size=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,size)
shadhap2 <-SH2 %>% mutate(size=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,size)
shadhap3 <-SH3 %>% mutate(size=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,size)
shadhap4 <-SH4 %>% mutate(size=rowSums(dplyr::select(.,-reads))) %>% dplyr::select(reads,size)

# Join molecule per read files
library(dplyr)
shap_join <- shadhap1 %>% full_join(shadhap2,by ="reads") %>%
  full_join(shadhap3,by ="reads") %>%
  full_join(shadhap4,by ="reads")

# Rename columns
shap_join<- shap_join %>% dplyr::rename(shadhap1 = size.x, shadhap2= size.y,
                                 shadhap3 = size.x.x, shadhap4= size.y.y )

# Data preparation
library("tidyverse")
shap_wrgl <- shap_join %>%
  dplyr::select(reads, shadhap1, shadhap2, shadhap3, shadhap4) %>%
  gather(key = "variable", value = "value", -reads)

head(shap_wrgl)

# Visualization
p<- ggplot(shap_wrgl, aes(x = reads, y = value)) +
  geom_line(aes(color = variable),size=0.4) +
  scale_color_manual(values=c("darkred", "steelblue","green","orange")) +
  scale_y_continuous(trans='log2') +
  labs(x ="molecule size", y = "count") +
  theme_bw()

# Save image to directory
p + ggsave(filename = paste0(datadir, "./molecule_size.png"))
```
