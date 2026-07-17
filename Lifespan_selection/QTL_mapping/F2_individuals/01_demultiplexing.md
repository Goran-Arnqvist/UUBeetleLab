# Read demultiplexing

Recommended readings: 

https://www.biostars.org/p/344768/
https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl_convert/bcl-convert-v3-7-5-software-guide-1000000163594-00.pdf
http://www.bea.ki.se/documents/bcl2fastq.pdf
https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/direct-demultiplexing
https://ipyrad.readthedocs.io/en/latest/5-demultiplexing.html


## Haplotagging background
Code and barcode files related to processing haplotagging data can be found [here](https://github.com/evolgenomics/haplotagging/blob/master/README.md).

## Software dependencies
- bcl2fastq v2.18.0 or above.
- bwa v0.6 or above, or now EMA
- libgzstream

## Strategy
Haplotagging uses a segmented combinatorial barcoding system in the standard Illumina Nextera indexing positions i5 and i7 to preserve linking information. To properly convert the data, our code expects the full set of R1, I1, I2, R2 fastq files, and assigns the barcode based on the look-up table segments A, B, C and D. It then encodes the barcode as comment fields BX, QX and RX (corresponding to barcode, quality strings, and corrected barcode tags, respectively) in a standard set of paired-end fastq files with R1 and R2.

These comment fields can then be passed into a BAM file as BX, QX and RX tags using standard software like bwa with a `-C` switch.

Example bcl2fastq command:
```bash
bcl2fastq --use-bases-mask=Y150,I13,I13,Y149 --create-fastq-for-index-reads -r [INT] -w [INT] -d [INT] -p [INT] -R <run_dir, e.g. 190125_ST-J00101_0130_AHYJWTBBXX> --tiles s_[1-8] --output-dir=<output_dir> --interop-dir=<INTEROPT_DIR>  --reports-dir=<REPORT_DIR>  --stats-dir=<STATS_DIR>
```
Here the options `--use-bases-mask=Y150,I13,I13,Y149` allows the full use of all 13 positions in the index reads. Note that a single cycle is taken out of R2 to extend the I2 cycle to 13nt.
```bash
--create-fastq-for-index-reads is key here to allow our demultiplexing code to see the full, untrimmed barcodes.
```

# 2022-09-14

## Sequence data of 1000s F2 beetles

The sequence data was delivered with mask=Y151,I8,Y28,Y136, so there is no need for further conversion for the tagging scripts. The sequence data of the beetles was partially demultiplexed by the sequencing facility (NGI) into 14 pools:
```bash
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleA/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleB/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleC/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleD/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleE/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleF/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleG/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleH/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleI/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleJ/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleK/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleL/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleM/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 14:40 Sample_VF-3324-BeetleN/
drwxrws--- 2 root snic2020-2-19 4.0K Aug 19 15:15 Undetermined/
```
The `Undetermined/` directory contains files that corresponds to the approx 5-10% of reads where there is issue with the sequence, either all poly-GGGGGG-containing sequences = no signal clusters, or reads with low quality or too many mistakes in the i7/plate barcodes so it was not possible to assign them to the BeetleA-N plates:
```bash
$ ls Undetermined
total 12G
-rw-rw---- 1 root snic2020-2-19 551M Aug 18 14:04 Undetermined_S0_L004_I1_001.fastq.gz
-rw-rw---- 1 root snic2020-2-19 4.6G Aug 18 14:04 Undetermined_S0_L004_R1_001.fastq.gz
-rw-rw---- 1 root snic2020-2-19 1.3G Aug 18 14:04 Undetermined_S0_L004_R2_001.fastq.gz
-rw-rw---- 1 root snic2020-2-19 4.8G Aug 18 14:04 Undetermined_S0_L004_R3_001.fastq.gz
```
For each of the pools (single directories), there are four files:
```bash
$ ls Sample_VF-3324-BeetleA

total 36G
-rw-rw---- 1 root snic2020-2-19 1.2G Aug 18 14:04 VF-3324-BeetleA_S120_L004_I1_001.fastq.gz
-rw-rw---- 1 root snic2020-2-19  16G Aug 18 14:04 VF-3324-BeetleA_S120_L004_R1_001.fastq.gz
-rw-rw---- 1 root snic2020-2-19 4.2G Aug 18 14:04 VF-3324-BeetleA_S120_L004_R2_001.fastq.gz
-rw-rw---- 1 root snic2020-2-19  15G Aug 18 14:04 VF-3324-BeetleA_S120_L004_R3_001.fastq.gz
```

## File and script set up

Marek explained that currently there are two methods for demultiplexing, one follows a two-step approach and the other one is a two-in-one approach.

1. Two-step:
  - Demultiplex the molecular BX-barcode (adds header tag e.g. BX:Z:A01C02B03D04) in each read using `demult_fastq_VX.o`
  - Then on the R1 and R2 outputs from the first script you would run the second script `vX_R2_clipping_postTag.o` to clip-off transposon sequence from the beginning of Read2

2. Two-in-one:
 - Does both BX-tagging and R2-clipping in a single step using `demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping.o`. Additionally re-demultiplexes the plate-barcode (based on I1.fastq.gz file) to also add the plate barcode into the BX-tag, so for example as BX:Z:A01C02B03D04-N701 instead of just BX:Z:A01C02B03D04. If correctly demultiplexed and clipped, all reads will have BX-tag with rarely containing A00 or B00 or C00 or D00 and the Read2 will be shorter, clipped to around 120 bp.

Marek recommends to use the two-in-one approach as it needs only one time execution, and it also adds plate barcode in the BX tag.

## Test run

Before running these scripts in Uppmax, it is necessary to verify whether the c++ binaries match the cluster's architecture (or is it needed to re-compile them?). So I did a test run of the first script of the two step approach on a single pool data set following these steps:

1. Upload to Uppmax the three directories provided by Marek to a single directory called `scripts/`:
```bash
drwx------ 2 angela snic2020-6-128 4.0K Sep 14 12:57 demult_fastq_vX/
drwx------ 2 angela snic2020-6-128 4.0K Sep 14 11:16 demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping/
drwx------ 2 angela snic2020-6-128 4.0K Sep 14 11:16 vX_R2_clipping_postTag/
```
I changed the executable permissions of the `*.o` and `*.cpp` files with `chmod +x <filename>`
```bash
chmod +x demult_fastq_VX.cpp
chmod +x demult_fastq_VX.o
```

2. Set bash environmental variables giving paths to specific directories and files:
```bash
SCRIPT_1_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_vX’
#SCRIPT_1_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping'
PATH_2_SAMPLE='../../Sample_VF-3324-BeetleA/'
```

3. Position within the directory `demult_fastq_vX/`and run the script:
```bash
cd ${SCRIPT_1_DIR}

./demult_fastq_VX.o ${PATH_2_SAMPLE}VF-3324-BeetleA_S120_L004_ VF-3324-BeetleA_S120_L004_tag
#./demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping.o ${PATH_2_SAMPLE}VF-3324-BeetleA_S120_L004_ VF-3324-BeetleA_S120_L004_tag
```

I left that running for a few hours and I got his message printed in screen:
```bash
loaded barcodes: 96 A, 96 B, 96 C, 96 D, 3 ME
```

These files were generated:
```bash
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_vX$ ls
total 6.6G
-rw-rw-r-- 1 angela angela 1.2K Jun  9 15:12 BC_A_VX.txt
-rw-rw-r-- 1 angela angela 1.1K Jun  9 15:12 BC_B.txt
-rw-rw-r-- 1 angela angela 1.1K Jun  9 15:12 BC_C_VX.txt
-rw-rw-r-- 1 angela angela 1.1K Jun  9 15:12 BC_D.txt
-rw-rw-r-- 1 angela angela   42 Jun  9 15:12 BC_ME.txt
-rwxrwxr-x 1 angela angela  11K Jun  9 15:13 demult_fastq_VX.cpp*
-rwxrwxr-x 1 angela angela  81K Jun  9 15:13 demult_fastq_VX.o*
-rw-rw-r-- 1 angela angela    0 Sep 14 12:57 VF-3324-BeetleA_S120_L004_tag_clearBC.log
-rw-rw-r-- 1 angela angela 3.4G Sep 14 15:09 VF-3324-BeetleA_S120_L004_tag_R1_001.fastq.gz
-rw-rw-r-- 1 angela angela 3.3G Sep 14 15:09 VF-3324-BeetleA_S120_L004_tag_R2_001.fastq.gz
-rw-rw-r-- 1 angela angela    0 Sep 14 12:57 VF-3324-BeetleA_S120_L004_tag_unclearBC.log
```

Marek recommended to examine the output R1 and R2 files. The BX-tag should have ACBD with almost no 00 in it.

- For R1:
```bash
$ zcat VF-3324-BeetleA_S120_L004_tag_R1_001.fastq.gz | more

@A00181:540:HMG5KDSX3:4:1101:1090:1000 BX:Z:A01C24B40D16	RX:Z:GTACTCCTCGTTCAGGATCCCCAACTAC	QX:Z:FFFFF:FFFFFFFFFFFFFFFFFFFFFF
GNTTATATTTCTTTTGTTTATTCCTATTCTTAAAATAACTGCTTTTTTAGAAACAAACCACCTAAAATATACGTGATTGATGAATTGGCGGAACAACATGGTCACAAAGTTCTTCGTCTGCCTCCTTACCATTGTATA
TTTAACCCTATTG
+
F#FFF:F:FFFF:FFFFFFF,:FFFFFFFF:FFFFFFFF:FFFFFFF:FFFFFFFFFF,F,FFFFFFFFF::,FFF:FFFF:F,FF:F:F:FFFFFFFF:F,FFFFFF:FFFFFF,:F::F,FFF:,F,FFFFFFFFF
FFFFF,,:F::F,
@A00181:540:HMG5KDSX3:4:1101:1127:1000 BX:Z:A68C73B32D40	RX:Z:TTCTCGCACTCTGATATCGGCCGGTTCT	QX:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFFF
CNACCATACTACAACCAGACATCACAAGACATTAGCATAACAACATGAATCAGTAGCTGTAAGGAATGCATTAAAAATGAGTAAATATATAATATGAACATATACAAAGAAAAGAAACAAAAGAACAAGATATCGTTC
AGCAAAAAAAAAT
+
F#FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF
FFFFFFFFFFFFF
@A00181:540:HMG5KDSX3:4:1101:1470:1000 BX:Z:A78C90B87D10	RX:Z:ACGCTTCGGAATGAGATCGCTCGTACCT	QX:Z:FFFFFFF:FFFFF:FFFFFFFFFFF:FF
ANTATGCCCTTTGTATAAATAGAACTTGAAACCTAATCGATTCGAATGTAAGCAAAAGAATTCGGTTAGTAAAAATTTCAAATAATATGAAAAATATATTGAACACGGTCTACCCCAGGGTTCGGTACTTGGGCCGAT
TTTATCACTGATA
+
F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
FFFFFFFFFFFFF

```
- For R2:
```bash
$ zcat VF-3324-BeetleA_S120_L004_tag_R2_001.fastq.gz | more

@A00181:540:HMG5KDSX3:4:1101:1090:1000 BX:Z:A01C24B40D16	RX:Z:GTACTCCTCGTTCAGGATCCCCAACTAC	QX:Z:FFFFF:FFFFFFFFFFFFFFFFFFFFFF
TTGTCAATAAGAGACAGCTGTTATAGTAACCTTTTGCTATGCCCCAAATTAACTCAATAGGGTTAAATATACAATGGTAAGGAGGCAGACGAAGAACTTTGTGCCCATGTTGTTCCGCCAATTCATCAATCACGTA
+
F:F,,,FFF:FFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00181:540:HMG5KDSX3:4:1101:1127:1000 BX:Z:A68C73B32D40	RX:Z:TTCTCGCACTCTGATATCGGCCGGTTCT	QX:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFFF
ACAACTGTATAAGAGACAGAAGTTGCACTTTCTGGAGCACATTGAGTGATAAATTCATTTTTTTTTGCTGAACGATATCTTGTTCTTTTGTTTCTTTTCTTTGTATATGTTCATATTATATATTTACTCATTTTTA
+
FF,FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFF,FFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@A00181:540:HMG5KDSX3:4:1101:1470:1000 BX:Z:A78C90B87D10	RX:Z:ACGCTTCGGAATGAGATCGCTCGTACCT	QX:Z:FFFFFFF:FFFFF:FFFFFFFFFFF:FF
ACAAGAGTATAAGAGTCAGCCCCAACTTTGCGAACACCCTGTACCGAAATAAATTGACTATTCTACTTATATGATAAGCTTCTAAGTAAATAAACGTTTTTTTCTTAGCTTTCGTGCTGTTTCGTCAATATGCATC
+
FF,FF,F,,,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFF:FFF:FFF:

```
All looks good!

### For demult_fastq_VX.o script

Now, let's test the two-in-one script.

The original c++ script shared by Marek was tailored for 10 bp barcodes. So he shared a new version of the script, this time for 8 bp long barcodes. To upload these scripts to Uppmax, I used:
```bash
rsync -av --progress /Users/angfu103/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/code/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping.cpp angela@rackham.uppmax.uu.se:/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping/

rsync -av --progress /Users/angfu103/Dropbox/PostDoc_UU/Projects/Beetles/QTL_mapping/offspring/code/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping.o angela@rackham.uppmax.uu.se:/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping/
```

He also shared a new `Plate_BC.txt` file, the one I had was missing one barcode:
```bash
# Old version

N701	GTCGAACA
N702	CCTGATGT
N703	TCCAAGAC
N704	TTCCTCTG
N705	GGTAGGTT
N706	AACGGTTC
N707	TAGGTGCT
N708	TGTGCAAG
N709	CCAATACG
N710	TCTTGCCA
N711	GATCTAGC
N712	CTGCGAAT
N713	GCACCTTA
PhiX	GGGGGGGG

# New version

N701	GTCGAACA
N702	CCTGATGT
N703	TCCAAGAC
N704	TTCCTCTG
N705	GGTAGGTT
N706	AACGGTTC
N707	TAGGTGCT
N708	TGTGCAAG
N709	CCAATACG
N710	TCTTGCCA
N711	GATCTAGC
N712	CTGCGAAT
N713	GCACCTTA
N714	CAAGACTA
```

I changed the executable permissions of the `*.o` and `*.cpp` files within `demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping/`:
```bash
chmod +x demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping.o
chmod +x demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping.cpp
```
Submit the job `01-demultiplex_Haplo_VX.sh` to Uppmax:
```bash
#!/bin/bash
#SBATCH -A snic2022-22-595
#SBATCH -M snowy
#SBATCH -p core -n 1
#SBATCH -t 20:00:00
#SBATCH -J demult_Haplo_VX
#SBATCH -e demult_Haplo_VX_%J_%A_%a.err
#SBATCH -o demult_Haplo_VX_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-14

# Set bash environmental variables giving paths to specific directories and files:
SCRIPT_DIR='/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping'
SAMPLE_LIST=${SCRIPT_DIR}/sample.list

echo This is array job: $SLURM_ARRAY_TASK_ID
SAMPLE_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE_LIST)
#SAMPLE_target=$(sed -n "1"p $SAMPLE_LIST)
echo -e This is the SAMPLE target: ${SAMPLE_target}

SAMPLE_PREFFIX=$(echo $SAMPLE_target | tr "/" "\t" | cut -f 4)
#echo $SAMPLE_PREFFIX

# Load required programs
module load gcc/12.2.0

# Position within the directory that has the c++ script
cd ${SCRIPT_DIR}

# Run the script
./demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping.o ${SAMPLE_target}_ ${SAMPLE_PREFFIX}_tag

```
Submitted batch job 6946964 on cluster snowy
**Runtime: 00-10:59:00**


These files were generated:
```bash
/proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_vX$ ls

```
Marek recommended to examine the output R1 and R2 files. The BX-tag should have ACBD with almost no 00 in it.

- For R1:
```bash
$ zcat VF-3324-BeetleN_S133_L004_tag_R1_001.fastq.gz | more

@A00181:540:HMG5KDSX3:4:1101:1199:1000 BX:Z:A58C15B67D85-N714	RX:Z:AGTACGCAGAGCAAATGGCGCGGAAGTA+CAAGACTA	QX:Z:,,FFFFFF,FFFF:,FFFFFF,FF:FF:+FFFFFFFF
ANGCAACATTCATAAAGGTGGATCAGGCTGCGCGGCCTTATTTCACTGCCACCAAGATCGTTGGACATAAATCCAATAGTTTCTTTTTATTAGGGTTTCTTTAGAGAACTTGTATACCAAGCAGTGCTGACCACAATAACAGAATTGAGTA
+
:#FFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFF,FFFF:FFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFF::FFFFF:FFFF,FFF:FF,
@A00181:540:HMG5KDSX3:4:1101:1253:1000 BX:Z:A48C23B55D74-N714	RX:Z:GTCTACCCGTAACAACTGACCGCAGATA+CAAGACTA	QX:Z:FFFFFFFFFFFFFFFF:FFFFFFFFFFF+FFFF:FFF
ANGCAAGACCCGCCATTCCAAACAGCCCCCAAACTGAGATGGGGCGTGCTGGGAAATGTTCCCTCAAAATCTTGTTGGAAGCATTGTCCAAATAATCTAGTCGTGCAAAATATAAATATACAAAATAATTCGTGTAACGTTTACAAACTTT
+
F#FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FF:FFF:FFFFFFF:FFFFFF:FFFFFFFFFFFFFFFFFFFFF:FFF:FF:FFFFFFFFF:FFFFF:FFFFFFFFF::FFFF

```
- For R2:
```bash
$ zcat VF-3324-BeetleN_S133_L004_tag_R2_001.fastq.gz | more

@A00181:540:HMG5KDSX3:4:1101:1199:1000 BX:Z:A58C15B67D85-N714	RX:Z:AGTACGCAGAGCAAATGGCGCGGAAGTA+CAAGACTA	QX:Z:,,FFFFFF,FFFF:,FFFFFF,FF:FF:+FFFFFFFF
CTATTATACAGGGTGTCCCGAAACCGCAGGATAATTCTGTAATAATTCTAGATAGATTATGTCATGGCCTACTAGATGGTAAAAATTTGTTTTATTAAAAGTCTTATAGTTTTCGAGA
+
FF::FFFFFFFF:FFFFFF::,FFFF,FFFFFF:FFFFFF:F:FFFFF,FFFFFFFFFF::F,FF,FFF:FFFFFFFF::F::F:FFFFFFF:FFFFF:FFFFFFFFFFFFF,FFF:,
@A00181:540:HMG5KDSX3:4:1101:1253:1000 BX:Z:A48C23B55D74-N714	RX:Z:GTCTACCCGTAACAACTGACCGCAGATA+CAAGACTA	QX:Z:FFFFFFFFFFFFFFFF:FFFFFFFFFFF+FFFF:FFF
CATTGAGGCAGTGAGTTTGAATGCTGCGATCACCACGGCGTTCTTCGCGTTGCCGTTGGTTTTTTGATCGTTCTGTGAGGAGACTTCTTCATGATGATTTTCATTTCCATCCTTAAAA
+
FFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFF,FFFFF:FF:FFFFFFFF:FFFFFFFF,F:FFFF:FF,FF,FFFFFFFFFFFFF,FFFFFFFFFF:,FFF:FFF,
@A00181:540:HMG5KDSX3:4:1101:1579:1000 BX:Z:A20C12B79D37-N714	RX:Z:TCAACCCCAGGTTTGTCCTCAGACGAAC+CAAGACTA	QX:Z:FFFFFFFFFFFFFFFFFFFFFFFFFFFF+FFFFFFFF
GGTTACTCTGCTGCTCTCGTATTTGTTCTAGGATTTTGTTTTACTTCTCACGAAACTGGAAATATCCTAATCGCGAAAAATTTCGAAAATAAAATAATATGTTTTTTGGATATTATTTG
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
All looks good!

Then I transferred the processed files to a new directory:
```bash
cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads
mkdir demultiplexed_reads

cd /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/scripts/demult_fastq_VX_plus_PlateBC_New_8bpPlateBC_clipping
mv VF-3324* /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/01-raw-reads/demultiplexed_reads
```
