# Data generation

## Samples and phenotypes

Naming of F2's/samples:

Samples/individuals are named by two letters and two numbers as e.g. "AB1-23". Here, "23" is just a consecutive number tagging the individual. "AB1" denotes the precise F2 family that individual "23" was from. The first letter represents the F1 male ("A") and the second the F1 female ("B") used to generate the F2 family.

The F1's, in turn, were from either of four distinct F0 families as:
A – LE1
B – EL4
C – LE2
D – EL1
Here, "LE1" means that the male was an L male and female was an E female (male first and female last; the number "1" is just a number tag for the F0 family of this type). Hence, we used 8 F0 individuals all together (4 males and 4 females).

## Library preparation

Library preparation was conducted by:

Marek Kucka, PhD
Chan-Lab Staff Scientist/Lab Manager
Friedrich Miescher Laboratory
of the Max Planck Society
Max-Planck-Ring 1
72076 Tuebingen
Germany
[website](http://fml.tuebingen.mpg.de/chan-group/group-members/)

The samples were prepared using *Haplotagging_vX* version of beads which *allow many plates sequenced in the same lane thanks to extra in-read2-barcoding* (attached pic).

Barcode structure - Ask Marek.


## Sequencing

**How to sequence:**

The sequencing needs to be done with *8-bp i7-index1*, and *13-bp of i5-index2*, so 151 + 8 + 13 + 151, and on *a single Novaseq S4 lane*.

The company could do partial demultiplexing into the 14 plates, based on the 8-bp long i7-barcode, which is here the plate barcode. The molecular and sample barcode, which has two parts in i5-index2 and another 2 parts in-line of read2 will be done by you post-data delivery using our demultiplexing script, which decodes the molecular barcode and adds it as BX:Z: tag (e.g. A01C01B01D01) for each read, where the for example here the 'D01' is the first sample barcode of each plate's 96 samples.

Partial demultiplexing can then be done by your sequencing company using mask: y151, i8, y28, y136, (where y28 is made out of 13-bp i5 index2 plus first 15-bp of read2, y136 is the leftover read2); and using a SampleSheet containing all 14 i7-barcodes which I will provide in later email.

Company will then deliver 14 sets of fastq files, delivering for each plate R1 R2 and R3 fastq.gz files.

Side note:
Sequencing companies using Novaseq S4, should be able to use up to 35 cycles for indexing or extra read-length sequencing; so with our 8 + 13 indexing setup there are at least 10 extra cycles un-used and potentially available to assign to either get longer read1 or read2 (longer than 151 bp). Because the Haplotagging_vX uses in-read2-barcoding, the final read2 is only 120 bp long after removal of the barcode information, so it could be useful to ask whether it would be possible to sequence read2 longer than 151-bp. Alternatively, because sequence quality of read2 is generally worse than that of read1, adding instead extra 10 cycles to have 161-bp read1 would be maybe even more beneficial.

In attachment:
I have also attached image explaining the difference between the traditional 'Haplotagging', you are used to hear about, and the 'vX_Haplotagging', 'vX_Haplotagging' was used to prepare your samples.


**To sequencing company:**

1 lane of Novaseq S4. 8+13 i7+i5 index sizes. Demultiplex using SampleSheet and mask: y151, i8, y28, y136; try getting extra 10-bp length for read1 or read2 if possible.


## Some relevant GitHub repos for linked reads

https://github.com/evolgenomics/haplotagging

https://github.com/wtsi-hpag/SamHaplotag

https://github.com/johandahlberg/awesome-10x-genomics

https://github.com/arshajii/ema

https://github.com/rwdavies/STITCH

https://github.com/rwdavies/QUILT

https://github.com/morispi/LEVIATHAN

https://github.com/adreau/ReMIX
