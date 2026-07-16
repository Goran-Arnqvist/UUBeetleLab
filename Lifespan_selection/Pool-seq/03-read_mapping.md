# Mapping of clean reads of the bean beetle

<!-- MDTOC maxdepth:6 firsth1:2 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Replace scaffold headers of the new HiC genome assembly](#replace-scaffold-headers-of-the-new-hic-genome-assembly)   
   - [1. Sort scaffolds by size and retrieve original scaffold names](#1-sort-scaffolds-by-size-and-retrieve-original-scaffold-names)   
   - [2. Make a key-value file listing the original and the new scaffold name pairs](#2-make-a-key-value-file-listing-the-original-and-the-new-scaffold-name-pairs)   
   - [3. Rename the fasta headers using SeqKit and the key-value file](#3-rename-the-fasta-headers-using-seqkit-and-the-key-value-file)   
   - [4. Sort-by-size the scaffolds in the new fasta file using Seqkit](#4-sort-by-size-the-scaffolds-in-the-new-fasta-file-using-seqkit)   
- [Obtain index files for the reference genome](#obtain-index-files-for-the-reference-genome)   
   - [Setup input files and scripts for read mapping](#setup-input-files-and-scripts-for-read-mapping)   
- [Generate a Multiqc report of the qualimap reports](#generate-a-multiqc-report-of-the-qualimap-reports)   
- [Get a depth of coverage histogram using BEDTools](#get-a-depth-of-coverage-histogram-using-bedtools)   
- [Plot coverage histogram in R](#plot-coverage-histogram-in-r)   

<!-- /MDTOC -->

# 2021-03-28

## Replace scaffold headers of the new HiC genome assembly

Looking at the scaffold naming, I noticed that there is not a distinction between putative chromosomes (9 autosomes + sex chromosome) and unplaced scaffolds, like in the Herring or in the horse mackerel assemblies (see below). Also, numbering does not necessarily make this distinction either, e.g. see “HiC_scaffold_7" is just 71634 bp in length.

Thinking on further analysis and plotting, I believe it'd be more convenient to sort the scaffolds by length and rename the 10 longest ones as “SUPER_” or “chr”. But before to proceed, I asked Mats and Calle for their advice.

- Bean beetle fasta.fai (head -n 20):
```bash
HiC_scaffold_1	6647717	16	80	81
HiC_scaffold_2	1187000	6730846	80	81
HiC_scaffold_3	89471000	7932700	80	81
HiC_scaffold_4	117235856	98522104	80	81
HiC_scaffold_5	96635042	217223425	80	81
HiC_scaffold_6	94445115	315066422	80	81
HiC_scaffold_7	71634	410692117	80	81
HiC_scaffold_8	74866084	410764663	80	81
HiC_scaffold_9	126082896	486566590	80	81
HiC_scaffold_10	128021113	614225540	80	81
HiC_scaffold_11	116575500	743846934	80	81
HiC_scaffold_12	108704056	861879645	80	81
HiC_scaffold_13	43142905	971942519	80	81
HiC_scaffold_14	3752528	1015624728	80	81
HiC_scaffold_15	851082	1019424180	80	81
HiC_scaffold_16	1096893	1020285918	80	81
HiC_scaffold_17	22297	1021396540	80	81
HiC_scaffold_18	27500	1021419133	80	81
HiC_scaffold_19	19500	1021446994	80	81
HiC_scaffold_20	79500	1021466755	80	81
```
- Herring fasta.fai:
```bash
chr1	33084258	6	80	81
chr2	33010319	33497824	80	81
chr3	32527562	66920778	80	81
chr4	32267647	99854941	80	81
chr5	31586861	132525940	80	81
chr6	31461554	164507643	80	81
chr7	30990621	196362473	80	81
chr8	30729556	227740483	80	81
chr9	30477381	258854165	80	81
chr10	30227731	289712521	80	81
chr11	30096327	320318106	80	81
chr12	30022480	350790645	80	81
chr13	29845739	381188413	80	81
chr14	29332771	411407231	80	81
chr15	28713521	441106669	80	81
chr16	27773822	470179117	80	81
chr17	27568510	498300119	80	81
chr18	27247294	526213243	80	81
chr19	27130643	553801136	80	81
chr20	26694162	581270920	80	81
chr21	26465981	608298767	80	81
chr22	25664052	635095580	80	81
chr23	25292897	661080440	80	81
chr24	20091098	686689506	80	81
chr25	14924191	707031750	80	81
chr26	12443209	722142501	80	81
unplaced_scaffold1	450286	734741271	80	81
unplaced_scaffold2	396991	735197206	80	81
unplaced_scaffold3	391441	735599180	80	81
unplaced_scaffold4	331978	735995535	80	81
```
- Atlantic horse mackerel fasta.fai:
```bash
SUPER_1	40754244	556216356	60	61
SUPER_2	15468355	497849993	60	61
SUPER_3	40663798	514874819	60	61
SUPER_4	36093002	271580025	60	61
SUPER_5	36888734	2086929	60	61
SUPER_6	36283071	679790459	60	61
SUPER_7	35447499	372867401	60	61
SUPER_8	33049929	41421645	60	61
SUPER_9	37400736	597649847	60	61
SUPER_10	38572790	76976386	60	61
SUPER_11	33386963	308846631	60	61
SUPER_12	29584266	342790054	60	61
SUPER_13	40695139	635673939	60	61
SUPER_14	33713437	760625655	60	61
SUPER_15	31049438	435708213	60	61
SUPER_16	35009248	721126722	60	61
SUPER_17	35709794	797516065	60	61
SUPER_18	29245829	467734632	60	61
SUPER_19	25328826	218332907	60	61
SUPER_20	25462409	244957256	60	61
SUPER_21	35833567	117860056	60	61
SUPER_22	30602059	186227198	60	61
SUPER_23	25353499	409472428	60	61
SUPER_24	29127247	155499372	60	61
scaffold_20_arrow_ctg1	46342	435248510	60	61
scaffold_61_arrow_ctg1	1035841	677047355	60	61
scaffold_66_arrow_ctg1	840293	678100485	60	61
scaffold_67_arrow_ctg1	821943	678954807	60	61
scaffold_70_arrow_ctg1	788755	716678272	60	61
scaffold_71_arrow_ctg1	759581	717480197	60	61
```

Thus, to make the header replacement I will follow these steps:
1. From the index file of the HiC assembly fasta (fai), sort scaffolds by size and retrieve original scaffolds names
2. Make a key-value file listing the original and the new scaffold name pairs
    -> Make a numerical vector from 1 to the total number of scaffolds
    -> Add the prefix string “chr_” to the 10 longest and “scaffold_” to the rest
3. Rename the fasta headers using SeqKit and the key-value file
4. Sort-by-size the scaffolds in the new fasta file using Seqkit


### 1. Sort scaffolds by size and retrieve original scaffold names

This is the total number of scaffolds in the HiC genome assembly:
```bash
cat a.obtectus_HiC_assembly.fasta | grep ">" | wc -l
3802
```
So, to create the key-value file, sort the index file of the fasta (fai) by decreasing size of scaffolds to get the scaffold names order:
```bash
cat a.obtectus_HiC_assembly.fasta.fai | sort -k2,2nr | head -n 20

HiC_scaffold_10	128021113	614225540	80	81
HiC_scaffold_9	126082896	486566590	80	81
HiC_scaffold_4	117235856	98522104	80	81
HiC_scaffold_11	116575500	743846934	80	81
HiC_scaffold_12	108704056	861879645	80	81
HiC_scaffold_5	96635042	217223425	80	81
HiC_scaffold_6	94445115	315066422	80	81
HiC_scaffold_3	89471000	7932700	80	81
HiC_scaffold_8	74866084	410764663	80	81
HiC_scaffold_13	43142905	971942519	80	81
HiC_scaffold_1	6647717	16	80	81
HiC_scaffold_14	3752528	1015624728	80	81
HiC_scaffold_2	1187000	6730846	80	81
HiC_scaffold_16	1096893	1020285918	80	81
HiC_scaffold_196	974500	1037170666	80	81
HiC_scaffold_246	960270	1045742342	80	81
HiC_scaffold_254	936000	1047557385	80	81
HiC_scaffold_195	897403	1036262027	80	81
HiC_scaffold_185	897000	1033844388	80	81
HiC_scaffold_15	851082	1019424180	80	81

cat a.obtectus_HiC_assembly.fasta.fai | sort -k2,2nr | cut -f1 > scaffold_names_sorted_by_size_a.obtectus_HiC_assembly.txt

head scaffold_names_sorted_by_size_a.obtectus_HiC_assembly.txt
HiC_scaffold_10
HiC_scaffold_9
HiC_scaffold_4
HiC_scaffold_11
HiC_scaffold_12
HiC_scaffold_5
HiC_scaffold_6
HiC_scaffold_3
HiC_scaffold_8
HiC_scaffold_13
```

### 2. Make a key-value file listing the original and the new scaffold name pairs
In Microsoft Excel, make a key-value guide file as follows:
- In col1, copy the of HiC scaffold names sorted by size
- In col2, make a series of numbers from 1 to 3802
- In col3, use the CONCATENATE() function to add the prefix "chr_" to the 10 largest scaffolds (which are the putative chromosomes in this case) and "scaffold_" to the rest
- Save this file in TXT format

To upload this file to Uppmax:
```bash
rsync -av --progress ~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/03-read_mapping/key-value_scaffold_rename_guide_a.obtectus_HiC_assembly.txt angela@rackham.uppmax.uu.se:/proj/snic2020-6-128/private/obtectus_poolseq/00-genome/

head key-value_scaffold_rename_guide_a.obtectus_HiC_assembly.txt
HiC_scaffold_10	chr_1
HiC_scaffold_9	chr_2
HiC_scaffold_4	chr_3
HiC_scaffold_11	chr_4
HiC_scaffold_12	chr_5
HiC_scaffold_5	chr_6
HiC_scaffold_6	chr_7
HiC_scaffold_3	chr_8
HiC_scaffold_8	chr_9
HiC_scaffold_13	chr_10
```

### 3. Rename the fasta headers using SeqKit and the key-value file

Replace the fasta headers using the SeqKit function [replace](https://bioinf.shenwei.me/seqkit/usage/#replace) and the key-value file just created:

replace name/sequence by regular expression.

Note that the replacement supports capture variables.
e.g. $1 represents the text of the first submatch.
ATTENTION: use SINGLE quote NOT double quotes in \*nix OS.

Examples: Adding space to all bases.
    seqkit replace -p '(.)' -r '$1 ' -s
Or use the \ escape character.
    seqkit replace -p '(.)' -r '\$1 ' -s
more on: http://bioinf.shenwei.me/seqkit/usage/#replace

Special replacement symbols (only for replacing name not sequence):
    {nr}    Record number, starting from 1
    {kv}    Corresponding value of the key (captured variable $n) by key-value file,
            n can be specified by flag -I (--key-capt-idx) (default: 1)

Usage:
  seqkit replace [flags]

Flags:
  -s, --by-seq                 replace seq
  -h, --help                   help for replace
  -i, --ignore-case            ignore case
  -K, --keep-key               keep the key as value when no value found for the key (only for sequence name)
  -I, --key-capt-idx int       capture variable index of key (1-based) (default 1)
  -m, --key-miss-repl string   replacement for key with no corresponding value
  -k, --kv-file string         tab-delimited key-value file for replacing key with value when using "{kv}" in -r (--replacement) (only for sequence name)
      --nr-width int           minimum width for {nr} in flag -r/--replacement. e.g., formating "1" to "001" by --nr-width 3 (default 1)
  -p, --pattern string         search regular expression
  -r, --replacement string     replacement. supporting capture variables.  e.g. $1 represents the text of the first submatch. ATTENTION: for \*nix OS, use SINGLE quote NOT double quotes or use the \ escape character. Record number is also supported by "{nr}".use ${1} instead of $1 when {kv} given!

Before making the header replacement on the actual fasta file, I tested the commands with toy files:
```bash
# Make toy files.
nano test.fa
>HiC_scaffold_10
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_9
CCCCAAAACCCCATGGCATCATTCA
>HiC_scaffold_4
CCCCAAAACCCCATGTTGCTACTAG
>HiC_scaffold_11
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_12
CCCCAAAACCCCATGGCATCATTCA

nano keyValue.txt
HiC_scaffold_10	chr_1
HiC_scaffold_9	chr_2
HiC_scaffold_4	chr_3
HiC_scaffold_11	chr_4
HiC_scaffold_12	chr_5

# Command to replace headers.
seqkit replace -p '(.+)$' -r '{kv}' -k keyValue.txt test.fa > new.test.fa

# Compare before and after files.
cat test.fa
>HiC_scaffold_10
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_9
CCCCAAAACCCCATGGCATCATTCA
>HiC_scaffold_4
CCCCAAAACCCCATGTTGCTACTAG
>HiC_scaffold_11
CCCCAAAACCCCATGATCATGGATC
>HiC_scaffold_12
CCCCAAAACCCCATGGCATCATTCA

cat new.test.fa
>chr_1
CCCCAAAACCCCATGATCATGGATC
>chr_2
CCCCAAAACCCCATGGCATCATTCA
>chr_3
CCCCAAAACCCCATGTTGCTACTAG
>chr_4
CCCCAAAACCCCATGATCATGGATC
>chr_5
CCCCAAAACCCCATGGCATCATTCA

```
All good. Now, make the header replacements on the actual fasta file (OTT HiC genome):
```bash
module load bioinfo-tools
module load SeqKit/0.15.0

seqkit replace -p '(.+)$' -r '{kv}' -k key-value_scaffold_rename_guide_a.obtectus_HiC_assembly.txt a.obtectus_HiC_assembly.fasta > tmp.fasta

# Stdout
[INFO] read key-value file: key-value_scaffold_rename_guide_a.obtectus_HiC_assembly.txt
[INFO] 3802 pairs of key-value loaded
```
Compare the previous and new headers:
```bash
cat a.obtectus_HiC_assembly.fasta | grep ">" | head
>HiC_scaffold_1
>HiC_scaffold_2
>HiC_scaffold_3
>HiC_scaffold_4
>HiC_scaffold_5
>HiC_scaffold_6
>HiC_scaffold_7
>HiC_scaffold_8
>HiC_scaffold_9
>HiC_scaffold_10

cat tmp.fasta | grep ">" | head
>scaffold_11
>scaffold_13
>chr_8
>chr_3
>chr_6
>chr_7
>scaffold_197
>chr_9
>chr_2
>chr_1

```
All good, the program runs extremely fast, it took just a few seconds to complete!

### 4. Sort-by-size the scaffolds in the new fasta file using Seqkit
Sort-by-size the actual scaffolds in the new fasta file with renamed headers by using the SeqKit function [sort](https://bioinf.shenwei.me/seqkit/usage/#sort):
```bash
Usage
sort sequences by id/name/sequence/length.

By default, all records will be readed into memory.
For FASTA format, use flag -2 (--two-pass) to reduce memory usage. FASTQ not
supported.
Firstly, seqkit reads the sequence head and length information.
If the file is not plain FASTA file,
seqkit will write the sequences to tempory files, and create FASTA index.
Secondly, seqkit sorts sequence by head and length information
and extracts sequences by FASTA index.

Usage:
  seqkit sort [flags]

Flags:
  -l, --by-length               by sequence length
  -n, --by-name                 by full name instead of just id
  -s, --by-seq                  by sequence
  -i, --ignore-case             ignore case
  -k, --keep-temp               keep tempory FASTA and .fai file when using 2-pass mode
  -N, --natural-order           sort in natural order, when sorting by IDs/full name
  -r, --reverse                 reverse the result
  -L, --seq-prefix-length int   length of sequence prefix on which seqkit sorts by sequences (0 for whole sequence) (default 10000)
  -2, --two-pass                two-pass mode read files twice to lower memory usage. (only for FASTA format)

Example of sort by sequence length:
$ echo -e ">seq1\nACGTNcccc\n>SEQ2\nacgtnAAAAnnn\n>seq3\nacgt" \
    | seqkit sort --quiet -l
```
To sort the scaffolds by length in descending order, use the flags `-l` and `-r`:
```bash
module load bioinfo-tools
module load SeqKit/0.15.0

seqkit sort --quiet -l -r tmp.fasta > A.obtectus_v2.0_sorted.fasta
```
Compare the previous and new headers order:
```bash
cat tmp.fasta | grep ">" | head
>scaffold_11
>scaffold_13
>chr_8
>chr_3
>chr_6
>chr_7
>scaffold_197
>chr_9
>chr_2
>chr_1

cat A.obtectus_v2.0_sorted.fasta | grep ">" | head -n 20
>chr_1
>chr_2
>chr_3
>chr_4
>chr_5
>chr_6
>chr_7
>chr_8
>chr_9
>chr_10
>scaffold_11
>scaffold_12
>scaffold_13
>scaffold_14
>scaffold_15
>scaffold_16
>scaffold_17
>scaffold_18
>scaffold_19
>scaffold_20
```
Finally, rename the final sorted fasta and eliminate temporal files:
```bash
mv A.obtectus_v2.0_sorted.fasta A.obtectus_v2.0.fasta
rm tmp.fasta
```
I created soft symbolic links to these files:
```
cd /proj/snic2020-6-128/private/beetle_genomes/acanthoscelides_obtectus/HiC_genome_assembly/A.obtectus_v2.0

cp -s /proj/snic2020-6-128/private/a_obtectus_poolseq/00-genome/*.* .
```

# 2021-03-29

## Obtain index files for the reference genome
For read mapping and variant calling using BWA and GATK, respectively, it is necessary to generate a group of indexes and dictionaries for the reference genome. For this, I used the bash script called `00-prepare_ref_genome_for_BWA_GATK.sh`:

```bash
#!/bin/bash
#SBATCH -A snic2017-7-378
#SBATCH -M snowy
#SBATCH -p core -n 8
#SBATCH -t 10:00:00
#SBATCH -J prep_genome
#SBATCH -e prep_genome_%J_%A_%a.err
#SBATCH -o prep_genome_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com

# Load required software.
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.4

# Set environment variables.
REF_FILE='/proj/snic2020-6-128/private/obtectus_poolseq/00-genome/A.obtectus_v2.0.fasta'
WORK_DIR='/proj/snic2020-6-128/private/obtectus_poolseq/00-genome'
PICARD_JAR='/sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar'

cd $WORK_DIR
# Generate BWA file index.
# where -a bwtsw specifies that we want to use the indexing algorithm that is capable of handling the whole human genome.
# Expected Result: This creates a collection of files used by BWA to perform the alignment.
bwa index -a bwtsw ${REF_FILE}

echo "################# bwa index done" ;

# Generate fasta file index.
# This creates a file called reference.fa.fai, with one record per line for each of the contigs in the FASTA reference file. Each record is composed of the contig name, size, location, basesPerLine and bytesPerLine.
samtools faidx ${REF_FILE}

echo "################# samtools faidx done" ;

# Generate sequence dictionary.
# Note that this is the new syntax for use with the latest version of Picard. Older versions used a slightly different syntax because all the tools were in separate jars, so you'd call e.g. java -jar CreateSequenceDictionary.jar directly.
# This creates a file called reference.dict formatted like a SAM header, describing the contents of your reference FASTA file.
java -Xmx48g -jar $PICARD_JAR CreateSequenceDictionary \
REFERENCE=${REF_FILE} \
OUTPUT=${REF_FILE}.dict

echo "################# picard dictionary done" ;

```
Submitted batch job 4333805 on cluster snowy
**Runtime: 00-00:25:52**

Rename genome dictionary file, because GATK will expect it to be named as `genome.dict` not `genome.fasta.dict`:
```
cd /proj/snic2020-6-128/private/obtectus_poolseq/00-genome
mv A.obtectus_v2.0.fasta.dict A.obtectus_v2.0.dict
```


### Setup input files and scripts for read mapping
I will use the read mapper [bwa](http://bio-bwa.sourceforge.net) to align clean reads against the reference genome. For this, I created various input files required by the script.

- Path to R1 files:
```bash
ls -lh /proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/*/*_R1*.gz | awk '{print $9}' > /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/R1_files_bwa.txt

cat /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/R1_files_bwa.txt

/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EI-F/EI-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EII-F/EII-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIII-F/EIII-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIII-M/EIII-M_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EII-M/EII-M_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EI-M/EI-M_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIV-F/EIV-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIV-M/EIV-M_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LV-F/LV-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVI-F/LVI-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVII-F/LVII-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVIII-F/LVIII-F_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVIII-M/LVIII-M_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVII-M/LVII-M_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVI-M/LVI-M_R1.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LV-M/LV-M_R1.fastq.gz
```
- Path to R2 files:
```bash
ls -lh /proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/*/*_R2*.gz | awk '{print $9}' > /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/R2_files_bwa.txt

cat /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/R2_files_bwa.txt

/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EI-F/EI-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EII-F/EII-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIII-F/EIII-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIII-M/EIII-M_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EII-M/EII-M_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EI-M/EI-M_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIV-F/EIV-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/EIV-M/EIV-M_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LV-F/LV-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVI-F/LVI-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVII-F/LVII-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVIII-F/LVIII-F_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVIII-M/LVIII-M_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVII-M/LVII-M_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LVI-M/LVI-M_R2.fastq.gz
/proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/LV-M/LV-M_R2.fastq.gz
```
- Sample names:
```bash
ls -lh /proj/snic2020-6-128/private/obtectus_poolseq/02-clean-reads/*/*.fastq.gz | awk '{print $9}' | sed -e 's+/+\t+g' | cut -f8 | sed 's/_R.*$//' | uniq > /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/sample_IDs_bwa.txt

cat /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/sample_IDs_bwa.txt

EI-F
EII-F
EIII-F
EIII-M
EII-M
EI-M
EIV-F
EIV-M
LV-F
LVI-F
LVII-F
LVIII-F
LVIII-M
LVII-M
LVI-M
LV-M
```

Script `03-1-read_mapping.sh` that launches 16 separate jobs (one for each sample):

```bash
#!/bin/bash
#SBATCH -A snic2017-7-378
#SBATCH -M snowy
#SBATCH -p core -n 8
#SBATCH -t 4-00:00:00
#SBATCH -J map_reads
#SBATCH -e map_reads_%J_%A_%a.err
#SBATCH -o map_reads_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH -a 1-16

# Load required software.
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.4

# Set environment variables.
# Files.
WORK_DIR='/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files'
R1_FILE='./R1_files_bwa.txt'
R2_FILE='./R2_files_bwa.txt'
ID_FILE='./sample_IDs_bwa.txt'
REF_FILE='/proj/snic2020-6-128/private/obtectus_poolseq/00-genome/A.obtectus_v2.0.fasta'
# Programs.
PICARD_JAR='/sw/apps/bioinfo/picard/2.20.4/rackham/picard.jar'

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
unset DISPLAY  # Turn display off to avoid problems with X11 in Uppmax.

if [ -d $WORK_DIR/Qualimap_results ]; then echo "Qualimap_results/ exists"; else mkdir $WORK_DIR/Qualimap_results; fi

echo -e $(date -u) ": Qualimap began..."

qualimap --java-mem-size=60G bamqc -bam $SNIC_TMP/${ID_target}.sort.MarkDup.RG.bam -ip -outdir $WORK_DIR/Qualimap_results/${ID_target} -outformat html

echo -e $(date -u) ": Qualimap ended..."

```
Submitted batch job 4334797 on cluster snowy
**Runtime: 01-03:56:00**


# 2021-03-30

## Generate a Multiqc report of the qualimap reports
```
cd /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/Qualimap_results

module load bioinfo-tools MultiQC/1.10
multiqc .

[INFO   ]         multiqc : This is MultiQC v1.10
[INFO   ]         multiqc : Template    : default
[INFO   ]         multiqc : Searching   : /crex/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/Qualimap_results
Searching   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% 800/800  
[INFO   ]        qualimap : Found 16 BamQC reports
[INFO   ]         multiqc : Compressing plot data
[INFO   ]         multiqc : Report      : multiqc_report.html
[INFO   ]         multiqc : Data        : multiqc_data
[INFO   ]         multiqc : MultiQC complete
```


## Get a depth of coverage histogram using BEDTools

Input files:
- List of BAM files:
```bash
cd /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files

ls -lh /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/*.bam | awk '{print $9}' > OTT_bam.list

cat OTT_bam.list

/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EI-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EII-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EIII-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EIII-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EII-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EI-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EIV-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/EIV-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LV-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LVI-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LVII-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LVIII-F.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LVIII-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LVII-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LVI-M.sort.MarkDup.RG.bam
/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/LV-M.sort.MarkDup.RG.bam
```
- List of sample names:
```bash
ls -lh /proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/*.bam | awk '{print $9}' | tr "/" "\t" | cut -f7 | sed 's/.sort.MarkDup.RG.bam//g' > sample_names_OTT_bams.list

cat sample_names_OTT_bams.list

EI-F
EII-F
EIII-F
EIII-M
EII-M
EI-M
EIV-F
EIV-M
LV-F
LVI-F
LVII-F
LVIII-F
LVIII-M
LVII-M
LVI-M
LV-M

```

Script `03-2-depth-of-coverage-hist.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2017-7-378
#SBATCH -M snowy
#SBATCH -p core -n 2
#SBATCH -t 10:00:00
#SBATCH -J depthCov-hist
#SBATCH -e depthCov-hist_%J_%A_%a.err
#SBATCH -o depthCov-hist_%J_%A_%a.out
#SBATCH --mail-user=angela.fuentespardo@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -a 1-16

# Load programs
module load bioinfo-tools
module load BEDTools/2.29.2

# Set environment variables.
WORK_DIR='/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files'
BAM_FILE_LIST='/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/OTT_bam.list'
ID_FILE_LIST='/proj/snic2020-6-128/private/obtectus_poolseq/03-bam-files/sample_names_OTT_bams.list'

# For a given job in the array, set the correspondent sample files.
BAM_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $BAM_FILE_LIST)
ID_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ID_FILE_LIST)

# Print current file info to stdout for future reference.
echo This is array job: $SLURM_ARRAY_TASK_ID
echo -e This is the BAM target: ${BAM_target}
echo -e This is the ID target: ${ID_target}

# Go to the working directory
cd $WORK_DIR

# Obtain depth of coverage histogram per BAM file.
bedtools genomecov -ibam $BAM_target > $ID_target.cov
grep "genome" $ID_target.cov > $ID_target.genomecov && rm $ID_target.cov

```
Submitted batch job 4351015 on cluster snowy
**Runtime: 00-00:49:00**


## Plot coverage histogram in R
I made the coverage plots using the R script `plot-depthCov.R`:
```R
# Plot coverage distribution based on stats obtained with bedtools genomecov function
# Developed by Fan Han
setwd('~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/03-read_mapping/coverage')

files <- list.files("~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/03-read_mapping/coverage/", pattern = "\\.genomecov$")

pdf('OTT_16_pools_bams_genomecov.pdf')
for(sample in files) {
  print(sample)
  path <- paste("~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/03-read_mapping/coverage/", sample, sep="")
  cov <- read.csv(path, header = FALSE, sep = "\t")
  #png(paste(path,".png", sep = ""))
  #png(paste(path,".png", sep = ""), res = 300)
  plot(cov$V2, cov$V5, xlab = "Coverage", main = sample, type = "o", cex = 0.2, xlim=c(0,150))
  #dev.off()
}
dev.off()

```
