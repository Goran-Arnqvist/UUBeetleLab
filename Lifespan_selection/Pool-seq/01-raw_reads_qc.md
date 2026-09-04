## 2020-06-09

### Merge raw data

Pools were sequenced in 8 lanes, but there were *16 "true" samples* (8 selection lines X 2 sexes = 16), and each such sample was subjected to a single library prep and then divided up into four technical replicates which were run in four different lanes. For example, the sample of females from line Early I (EI-F) was run in lanes 5-8. 

Create symlinks to raw data:
```bash
cd <REDACTED_PATH>

ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .
ln -s <REDACTED_PATH> .

```
Verify the number of files is correct:
```bash
ls *.gz | wc -l
128
```
Concatenate reads per sample using the script `01-merge-raw-reads-per-pool.sh`:
```bash
#!/bin/bash

cd <REDACTED_PATH>

cat unmerged/EI-F_S9_L005_R1_001.fastq.gz unmerged/EI-F_S9_L006_R1_001.fastq.gz unmerged/EI-F_S9_L007_R1_001.fastq.gz unmerged/EI-F_S9_L008_R1_001.fastq.gz > EI-F_S9_L005-8_R1_001.fastq.gz
cat unmerged/EI-F_S9_L005_R2_001.fastq.gz unmerged/EI-F_S9_L006_R2_001.fastq.gz unmerged/EI-F_S9_L007_R2_001.fastq.gz unmerged/EI-F_S9_L008_R2_001.fastq.gz > EI-F_S9_L005-8_R2_001.fastq.gz
cat unmerged/EI-M_S1_L001_R1_001.fastq.gz unmerged/EI-M_S1_L002_R1_001.fastq.gz unmerged/EI-M_S1_L003_R1_001.fastq.gz unmerged/EI-M_S1_L004_R1_001.fastq.gz > EI-M_S1_L001-4_R1_001.fastq.gz
cat unmerged/EI-M_S1_L001_R2_001.fastq.gz unmerged/EI-M_S1_L002_R2_001.fastq.gz unmerged/EI-M_S1_L003_R2_001.fastq.gz unmerged/EI-M_S1_L004_R2_001.fastq.gz > EI-M_S1_L001-4_R2_001.fastq.gz
cat unmerged/EII-F_S11_L005_R1_001.fastq.gz unmerged/EII-F_S11_L006_R1_001.fastq.gz unmerged/EII-F_S11_L007_R1_001.fastq.gz unmerged/EII-F_S11_L008_R1_001.fastq.gz > EII-F_S11_L005-8_R1_001.fastq.gz
cat unmerged/EII-F_S11_L005_R2_001.fastq.gz unmerged/EII-F_S11_L006_R2_001.fastq.gz unmerged/EII-F_S11_L007_R2_001.fastq.gz unmerged/EII-F_S11_L008_R2_001.fastq.gz > EII-F_S11_L005-8_R2_001.fastq.gz
cat unmerged/EII-M_S3_L001_R1_001.fastq.gz unmerged/EII-M_S3_L002_R1_001.fastq.gz unmerged/EII-M_S3_L003_R1_001.fastq.gz unmerged/EII-M_S3_L004_R1_001.fastq.gz > EII-M_S3_L001-4_R1_001.fastq.gz
cat unmerged/EII-M_S3_L001_R2_001.fastq.gz unmerged/EII-M_S3_L002_R2_001.fastq.gz unmerged/EII-M_S3_L003_R2_001.fastq.gz unmerged/EII-M_S3_L004_R2_001.fastq.gz > EII-M_S3_L001-
4_R2_001.fastq.gz
cat unmerged/EIII-F_S13_L005_R1_001.fastq.gz unmerged/EIII-F_S13_L006_R1_001.fastq.gz unmerged/EIII-F_S13_L007_R1_001.fastq.gz unmerged/EIII-F_S13_L008_R1_001.fastq.gz > EIII-F_S13_L005-8_R1_001.fastq.gz
cat unmerged/EIII-F_S13_L005_R2_001.fastq.gz unmerged/EIII-F_S13_L006_R2_001.fastq.gz unmerged/EIII-F_S13_L007_R2_001.fastq.gz unmerged/EIII-F_S13_L008_R2_001.fastq.gz > EIII-F_S13_L005-8_R2_001.fastq.gz
cat unmerged/EIII-M_S5_L001_R1_001.fastq.gz unmerged/EIII-M_S5_L002_R1_001.fastq.gz unmerged/EIII-M_S5_L003_R1_001.fastq.gz unmerged/EIII-M_S5_L004_R1_001.fastq.gz > EIII-M_S5_L001-4_R1_001.fastq.gz
cat unmerged/EIII-M_S5_L001_R2_001.fastq.gz unmerged/EIII-M_S5_L002_R2_001.fastq.gz unmerged/EIII-M_S5_L003_R2_001.fastq.gz unmerged/EIII-M_S5_L004_R2_001.fastq.gz > EIII-M_S5_L001-4_R2_001.fastq.gz
cat unmerged/EIV-F_S15_L005_R1_001.fastq.gz unmerged/EIV-F_S15_L006_R1_001.fastq.gz unmerged/EIV-F_S15_L007_R1_001.fastq.gz unmerged/EIV-F_S15_L008_R1_001.fastq.gz > EIV-F_S15_L005-8_R1_001.fastq.gz
cat unmerged/EIV-F_S15_L005_R2_001.fastq.gz unmerged/EIV-F_S15_L006_R2_001.fastq.gz unmerged/EIV-F_S15_L007_R2_001.fastq.gz unmerged/EIV-F_S15_L008_R2_001.fastq.gz > EIV-F_S15_L005-8_R2_001.fastq.gz
cat unmerged/EIV-M_S7_L001_R1_001.fastq.gz unmerged/EIV-M_S7_L002_R1_001.fastq.gz unmerged/EIV-M_S7_L003_R1_001.fastq.gz unmerged/EIV-M_S7_L004_R1_001.fastq.gz > EIV-M_S7_L001-4_R1_001.fastq.gz
cat unmerged/EIV-M_S7_L001_R2_001.fastq.gz unmerged/EIV-M_S7_L002_R2_001.fastq.gz unmerged/EIV-M_S7_L003_R2_001.fastq.gz unmerged/EIV-M_S7_L004_R2_001.fastq.gz > EIV-M_S7_L001-4_R2_001.fastq.gz
cat unmerged/LV-F_S10_L005_R1_001.fastq.gz unmerged/LV-F_S10_L006_R1_001.fastq.gz unmerged/LV-F_S10_L007_R1_001.fastq.gz unmerged/LV-F_S10_L008_R1_001.fastq.gz > LV-F_S10_L005-8_R1_001.fastq.gz
cat unmerged/LV-F_S10_L005_R2_001.fastq.gz unmerged/LV-F_S10_L006_R2_001.fastq.gz unmerged/LV-F_S10_L007_R2_001.fastq.gz unmerged/LV-F_S10_L008_R2_001.fastq.gz > LV-F_S10_L005-8_R2_001.fastq.gz
cat unmerged/LV-M_S2_L001_R1_001.fastq.gz unmerged/LV-M_S2_L002_R1_001.fastq.gz unmerged/LV-M_S2_L003_R1_001.fastq.gz unmerged/LV-M_S2_L004_R1_001.fastq.gz > LV-M_S2_L001-4_R1_001.fastq.gz
cat unmerged/LV-M_S2_L001_R2_001.fastq.gz unmerged/LV-M_S2_L002_R2_001.fastq.gz unmerged/LV-M_S2_L003_R2_001.fastq.gz unmerged/LV-M_S2_L004_R2_001.fastq.gz > LV-M_S2_L001-4_R2_001.fastq.gz
cat unmerged/LVI-F_S12_L005_R1_001.fastq.gz unmerged/LVI-F_S12_L006_R1_001.fastq.gz unmerged/LVI-F_S12_L007_R1_001.fastq.gz unmerged/LVI-F_S12_L008_R1_001.fastq.gz > LVI-F_S12_L005-8_R1_001.fastq.gz
cat unmerged/LVI-F_S12_L005_R2_001.fastq.gz unmerged/LVI-F_S12_L006_R2_001.fastq.gz unmerged/LVI-F_S12_L007_R2_001.fastq.gz unmerged/LVI-F_S12_L008_R2_001.fastq.gz > LVI-F_S12_L005-8_R2_001.fastq.gz
cat unmerged/LVI-M_S4_L001_R1_001.fastq.gz unmerged/LVI-M_S4_L002_R1_001.fastq.gz unmerged/LVI-M_S4_L003_R1_001.fastq.gz unmerged/LVI-M_S4_L004_R1_001.fastq.gz > LVI-M_S4_L001-4_R1_001.fastq.gz
cat unmerged/LVI-M_S4_L001_R2_001.fastq.gz unmerged/LVI-M_S4_L002_R2_001.fastq.gz unmerged/LVI-M_S4_L003_R2_001.fastq.gz unmerged/LVI-M_S4_L004_R2_001.fastq.gz > LVI-M_S4_L001-4_R2_001.fastq.gz
cat unmerged/LVII-F_S14_L005_R1_001.fastq.gz unmerged/LVII-F_S14_L006_R1_001.fastq.gz unmerged/LVII-F_S14_L007_R1_001.fastq.gz unmerged/LVII-F_S14_L008_R1_001.fastq.gz > LVII-F_S14_L005-8_R1_001.fastq.gz
cat unmerged/LVII-F_S14_L005_R2_001.fastq.gz unmerged/LVII-F_S14_L006_R2_001.fastq.gz unmerged/LVII-F_S14_L007_R2_001.fastq.gz unmerged/LVII-F_S14_L008_R2_001.fastq.gz > LVII-F_S14_L005-8_R2_001.fastq.gz
cat unmerged/LVII-M_S6_L001_R1_001.fastq.gz unmerged/LVII-M_S6_L002_R1_001.fastq.gz unmerged/LVII-M_S6_L003_R1_001.fastq.gz unmerged/LVII-M_S6_L004_R1_001.fastq.gz > LVII-M_S6_L001-4_R1_001.fastq.gz
cat unmerged/LVII-M_S6_L001_R2_001.fastq.gz unmerged/LVII-M_S6_L002_R2_001.fastq.gz unmerged/LVII-M_S6_L003_R2_001.fastq.gz unmerged/LVII-M_S6_L004_R2_001.fastq.gz > LVII-M_S6_L001-4_R2_001.fastq.gz
cat unmerged/LVIII-F_S16_L005_R1_001.fastq.gz unmerged/LVIII-F_S16_L006_R1_001.fastq.gz unmerged/LVIII-F_S16_L007_R1_001.fastq.gz unmerged/LVIII-F_S16_L008_R1_001.fastq.gz > LVIII-F_S16_L005-8_R1_001.fastq.gz
cat unmerged/LVIII-F_S16_L005_R2_001.fastq.gz unmerged/LVIII-F_S16_L006_R2_001.fastq.gz unmerged/LVIII-F_S16_L007_R2_001.fastq.gz unmerged/LVIII-F_S16_L008_R2_001.fastq.gz > LVIII-F_S16_L005-8_R2_001.fastq.gz
cat unmerged/LVIII-M_S8_L001_R1_001.fastq.gz unmerged/LVIII-M_S8_L002_R1_001.fastq.gz unmerged/LVIII-M_S8_L003_R1_001.fastq.gz unmerged/LVIII-M_S8_L004_R1_001.fastq.gz > LVIII-M_S8_L001-4_R1_001.fastq.gz
cat unmerged/LVIII-M_S8_L001_R2_001.fastq.gz unmerged/LVIII-M_S8_L002_R2_001.fastq.gz unmerged/LVIII-M_S8_L003_R2_001.fastq.gz unmerged/LVIII-M_S8_L004_R2_001.fastq.gz > LVIII-M_S8_L001-4_R2_001.fastq.gz

```
Verify the merged fastq files were generated:
```bash
<REDACTED_PATH> ls

total 410G
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  12G Jun  9 12:17 EI-F_S9_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  12G Jun  9 12:17 EI-F_S9_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 15:26 EII-F_S11_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 15:28 EII-F_S11_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  14G Jun  9 12:22 EIII-F_S13_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  15G Jun  9 12:23 EIII-F_S13_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  11G Jun  9 12:24 EIII-M_S5_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  11G Jun  9 12:24 EIII-M_S5_L001-4_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE> 9.9G Jun  9 12:21 EII-M_S3_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE> 9.7G Jun  9 12:21 EII-M_S3_L001-4_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 12:18 EI-M_S1_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 12:19 EI-M_S1_L001-4_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  11G Jun  9 12:25 EIV-F_S15_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  11G Jun  9 12:25 EIV-F_S15_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  11G Jun  9 12:26 EIV-M_S7_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  12G Jun  9 12:26 EIV-M_S7_L001-4_R2_001.fastq.gz
drwxrwsr-x 2 angela <REDACTED_PROJECT_CODE> 4.0K Jun  9 15:07 FastQC/
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE> 2.6K Jun  9 15:11 fastq_files.txt
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 12:27 LV-F_S10_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 12:28 LV-F_S10_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 12:30 LVI-F_S12_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  13G Jun  9 12:31 LVI-F_S12_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  16G Jun  9 12:33 LVII-F_S14_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  16G Jun  9 12:34 LVII-F_S14_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  17G Jun  9 12:36 LVIII-F_S16_L005-8_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  17G Jun  9 12:37 LVIII-F_S16_L005-8_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  18G Jun  9 12:38 LVIII-M_S8_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  18G Jun  9 12:39 LVIII-M_S8_L001-4_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  14G Jun  9 12:34 LVII-M_S6_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  14G Jun  9 12:35 LVII-M_S6_L001-4_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  14G Jun  9 12:31 LVI-M_S4_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  14G Jun  9 12:32 LVI-M_S4_L001-4_R2_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  11G Jun  9 12:28 LV-M_S2_L001-4_R1_001.fastq.gz
-rw-rw-r-- 1 angela <REDACTED_PROJECT_CODE>  11G Jun  9 12:29 LV-M_S2_L001-4_R2_001.fastq.gz
drwxrwsr-x 2 angela <REDACTED_PROJECT_CODE>  12K Jun  9 15:25 unmerged/

```

## FastQC of raw reads

Input file `fastq_files.txt`:
```bash
ls -lh <REDACTED_PATH> | grep ".fastq.gz$" | awk '{print "<REDACTED_PATH>"$9}' > <REDACTED_PATH>

cat fastq_files.txt

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
I created a directory for the FastQC results, `mkdir FastQC_results`. The I run the script `01-raw-reads-qc.sh`:
```bash
#!/bin/bash
#SBATCH -A <REDACTED_PROJECT_CODE>
#SBATCH -p core -n 1
#SBATCH -t 5:00:00
#SBATCH -J raw-reads-qc
#SBATCH -e raw-reads-qc_%J_%A_%a.err
#SBATCH -o raw-reads-qc_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=<REDACTED_EMAIL>
#SBATCH -a 1-32

# Load modules.
module load bioinfo-tools
module load FastQC/0.11.8

# Set environment variables.
WORK_DIR='<REDACTED_PATH>'
SAMPLES_FILE='./fastq_files.txt'
OUTPUT_DIR='<REDACTED_PATH>'

cd $WORK_DIR

echo This is array job: $SLURM_ARRAY_TASK_ID
SAMPLE_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLES_FILE)
echo -e This is the sample target: ${SAMPLE_target}
echo -e This is the output file dir: ${OUTPUT_DIR}

# Run FastQC for each fastq file.
fastqc ${SAMPLE_target} -o ${OUTPUT_DIR}

```
Submitted batch job <REDACTED_JOB_ID>
**Runtime: 00-00:38:00**

Obtain a single FastQC report for all samples using MultiQC:
```bash
cd <REDACTED_PATH>

module load bioinfo-tools MultiQC/1.9
multiqc .
```
Stdout:
```bash
[INFO   ]         multiqc : This is MultiQC v1.9
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching   : <REDACTED_PATH>
Searching 64 files..  [####################################]  100%             
[INFO   ]          fastqc : Found 32 reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```

Overall, the sequence quality is good but there is a list overrepresented sequences, such:
```bash
Overrepresented sequences
Sequence	Count	Percentage	Possible Source
GTAATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGTAATAAATC	221057	0.15958714405812097	No Hit
AATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGT	220650	0.15929331953489095	No Hit
GATTTATTACGATTTATTACGATTTATTACGATTTATTACGATTTATTAC	211025	0.15234476662066784	No Hit
TAATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGTAATAAATCG	171818	0.12404015216789437	No Hit
AAATCGTAATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGTAAT	162346	0.11720205417272334	No Hit
TATTACGATTTATTACGATTTATTACGATTTATTACGATTTATTACGATT	152060	0.1097763071310923	No Hit
TTATTACGATTTATTACGATTTATTACGATTTATTACGATTTATTACGAT	150301	0.10850643652578129	No Hit
ATTTATTACGATTTATTACGATTTATTACGATTTATTACGATTTATTACG	150241	0.10846312087125108	No Hit
TACGATTTATTACGATTTATTACGATTTATTACGATTTATTACGATTTAT	145510	0.10504768151154308	No Hit
TTTATTACGATTTATTACGATTTATTACGATTTATTACGATTTATTACGA	145459	0.1050108632051924	No Hit
ATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGTA	141510	0.10215997120952829	No Hit
TTACGATTTATTACGATTTATTACGATTTATTACGATTTATTACGATTTA	140337	0.10131315016346244	No Hit
CGTAATAAATCGTAATAAATCGTAATAAATCGTAATAAATCGTAATAAAT	139198	0.10049087465496373	No Hit
```
As they do not hit with Illumina adapters and look like microsatellite repeats, perhaps they are just natural occurrences in the DNA. Let's see what happens after adapter trimming and removal of low quality bases.
