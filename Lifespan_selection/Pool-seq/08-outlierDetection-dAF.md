# Calculate pairwise delta allele frequencies (dAF) and generate Manhattan plots

<!-- MDTOC maxdepth:6 firsth1:2 numbering:0 flatten:0 bullets:1 updateOnSave:1 -->

- [Outlier SNP detection based on paired comparisons and calculation of delta allele frequencies](#outlier-snp-detection-based-on-paired-comparisons-and-calculation-of-delta-allele-frequencies)   
- [Script to calculate dAF](#script-to-calculate-daf)   
- [Obtain Manhattan plots of dAF](#obtain-manhattan-plots-of-daf)   
- [Find informative genes with high dAF differences between Early and Late reproduction](#find-informative-genes-with-high-daf-differences-between-early-and-late-reproduction)   

<!-- /MDTOC -->


# 2021-04-01
## Outlier SNP detection based on paired comparisons and calculation of delta allele frequencies

First, I defined the groupings for paired comparisons using the Excel file `~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/08-2-outlier-detection_dAF/List_paired_comparisons_OTT_2020-09-01.xlsx`. Then, I saved this file as a CVS file, which looks like this (I manually added the last new line while opening the file with Textwrangler):
```bash
Comparison;Group1;Group2
EI.II.III.IV_vs_LV.VI.VII.VIII;EI-F EI-M EII-F EII-M EIII-F EIII-M EIV-F EIV-M;LV-F LV-M LVI-F LVI-M LVII-F LVII-M LVIII-F LVIII-M
EI.II.III.IV_vs_LV;EI-F EI-M EII-F EII-M EIII-F EIII-M EIV-F EIV-M;LV-F LV-M
EI.II.III.IV_vs_LVI;EI-F EI-M EII-F EII-M EIII-F EIII-M EIV-F EIV-M;LVI-F LVI-M
EI.II.III.IV_vs_LVII;EI-F EI-M EII-F EII-M EIII-F EIII-M EIV-F EIV-M;LVII-F LVII-M
EI.II.III.IV_vs_LVIII;EI-F EI-M EII-F EII-M EIII-F EIII-M EIV-F EIV-M;LVIII-F LVIII-M
LV_vs_LVI;LV-F LV-M;LVI-F LVI-M
LV_vs_LVII;LV-F LV-M;LVII-F LVII-M
LV_vs_LVIII;LV-F LV-M;LVIII-F LVIII-M
LVI_vs_LVII;LVI-F LVI-M;LVII-F LVII-M
LVI_vs_LVIII;LVI-F LVI-M;LVIII-F LVIII-M
LVII_vs_LVIII;LVII-F LVII-M;LVIII-F LVIII-M
EI.EII.EIII.EIV.F_vs_EI.EII.EIII.EIV.M;EI-F EII-F EIII-F EIV-F;EI-M EII-M EIII-M EIV-M
EI.F_vs_EI.M;EI-F;EI-M
EII.F_vs_EII.M;EII-F;EII-M
EIII.F_vs_EIV.M;EIII-F;EIII-M
EIV.F_vs_EIV.M;EIV-F;EIV-M
LV.F_vs_LV.M;LV-F;LV-M
LVI.F_vs_LVI.M;LVI-F;LVI-M
LVII.F_vs_LVII.M;LVII-F;LVII-M
LVIII.F_vs_LVIII.M;LVIII-F;LVIII-M

```
Note that *samples within a group are separated by a single space*.

## Script to calculate dAF

The R script is called `08-2-calc-dAF_v5.R`:
```R
#####################################################################
#                                                                   
# Calculate the absolute delta allele frequency between pairs of    
# super pools                         
#                                                                   
# E-mail: apfuentesp@gmail.com, Uppsala University                  
# Date: 2020-03-25, modified: 2021-01-13                          
#                                                                   
#####################################################################

## Clean environment space.
rm(list = ls())

## ------------------------------------------------------------------
## Function definition
## ------------------------------------------------------------------

###############################################################################################
#'
#' Function that creates a CHR column of unique numbers extracted from the SNP or CHROM info
#' Date: 2020-05-05, updated on: 2020-09-01
#'
createCHR.uniqueNumbers <- function(df = NULL, chrIdentifier = NULL, scaffoldIdentifier = NULL,
                                    scaffSuffix = NULL, delimiterChrScaffNumber = NULL,
                                    delimiterChrPos = NULL) {
  ## For debugging.
  # For OTT genome v.2
  #chrIdentifier <- 'chr'
  #scaffoldIdentifier <- 'scaffold'
  #scaffSuffix <- NULL
  #delimiterChrScaffNumber <- '_'
  #delimiterChrPos <- '-'
  #df <- summary_dAF_df
  #head(df)

  # For Her genome v.2.0
  #chrIdentifier <- 'chr'
  #scaffoldIdentifier <- 'unplaced_scaffold'
  #scaffSuffix <- NULL
  #delimiterChrScaffNumber <- NULL
  #delimiterChrPos <- '-'
  #df <- summary_dAF_df

  # Load library
  library(gtools)

  if(!is.null(df$SNP) & is.null(df$CHROM) & is.null(df$POS) & !is.null(delimiterChrPos)) {  # If only the SNP column is available

    cat('\n# ------ SNP column was split into CHROM and POS columns ------ #')
    # Split the SNP name into two columns by the character '-'.
    library(dplyr)
    library(tidyr)
    df <- separate(data = df, col = SNP, into = c('CHROM', 'POS'), sep = delimiterChrPos, remove = FALSE)
    cat('\n')
    head(df)
    mixedsort(unique(df$CHROM))
  } # End, if only the SNP column is available

  if(!is.null(df$CHROM) & !is.null(df$POS)) {  # If CHROM and POS columns are available

    # Copy the CHROM and POS columns to the new CHR and BP columns (perhaps useful later for plotting).
    df$CHR <- df$CHROM
    df$BP <- df$POS

    # If scaffolds have a suffix string, remove it.
    if(!is.null(scaffSuffix)) {
      df$CHR <- gsub(scaffSuffix, '', df$CHR)
      #df$CHR <- gsub('\\|arrow\\|arrow', '', df$CHR)  # Only for OTT genome v.1.0
      #df$CHR <- gsub('\\R|arrow\\|arrow', '', df$CHR)  # Only for OTT genome v.1.0
      mixedsort(unique(df$CHR))
    }

    if(!is.null(delimiterChrScaffNumber)) {

      # Verify Chr and scaffolds numbers do not repeat (e.g. Chr_3 and scaffold_3). If so, modify scaffold numbering so it does not overlap with Chr numbering.
      tmp1 <- unique(df$CHR[grepl(chrIdentifier, df$CHR)])  # Extract Chr only
      tmp1_split <- strsplit(tmp1, delimiterChrScaffNumber, fixed = FALSE)
      tmp1_split <- matrix(unlist(tmp1_split), ncol = 2, byrow = TRUE)
      tmp1_split

      if (!is.null(scaffoldIdentifier)) {  # If there are scaffolds.
        tmp2 <- unique(df$CHR[grepl(scaffoldIdentifier, df$CHR)])  # Extract scaffolds only
        tmp2_split <- strsplit(tmp2, delimiterChrScaffNumber, fixed = FALSE)
        tmp2_split <- matrix(unlist(tmp2_split), ncol = 2, byrow = TRUE)
        tmp2_split

        if(any(tmp1_split[,2] %in% tmp2_split[,2])) {  # If any number is repeated, add '9' before all scaffolds, so 10 becomes 910, and so on.
          df$CHR <- gsub(paste0(scaffoldIdentifier, delimiterChrScaffNumber), paste0(scaffoldIdentifier, delimiterChrScaffNumber,'9'), df$CHR)
          print('\n# ------ A number 9 was added at the beginning of each scaffold because numbering was repeated between Chr and scaffolds ------ #')
          print(head(unique(df$CHR[grepl(scaffoldIdentifier, df$CHR)])))
          df$CHR <- gsub(paste0(scaffoldIdentifier, delimiterChrScaffNumber), '', df$CHR)  # Delete the scaffold name, just leave the number
        }
      }
      cat('\n# ------ Deleting chrIdentifier, if delimiterChrScaffNumber ------ #\n')
      df$CHR <- gsub(paste0(chrIdentifier, delimiterChrScaffNumber), '', df$CHR)  # Delete the chromosome name, just leave the number
      df$CHR <- gsub(paste0(scaffoldIdentifier, delimiterChrScaffNumber), '', df$CHR)  # Delete the scaffold name, just leave the number

      head(df)
      cat('\n# ------ Final numbering of Chr and scaffolds ------ #\n')
      cat(mixedsort(unique(df$CHR)))

      # Convert columns CHR and BP from character to numeric values.
      df$CHR <- as.numeric(df$CHR)
      df$BP <- as.numeric(df$BP)

      return(df)

    } else {  # If there is no separator between Chr/scaffold string to the number

      # Verify Chr and scaffolds numbers do not repeat (e.g. Chr_3 and scaffold_3). If so, modify scaffold numbering so it does not overlap with Chr numbering.
      tmp1 <- unique(df$CHR[grepl(chrIdentifier, df$CHR)])  # Extract Chr only
      tmp1_split <- gsub(chrIdentifier, '', tmp1)
      tmp1_split

      if (!is.null(scaffoldIdentifier)) {  # If there are scaffolds.
        tmp2 <- unique(df$CHR[grepl(scaffoldIdentifier, df$CHR)])  # Extract scaffolds only
        tmp2_split <- gsub(scaffoldIdentifier, '', tmp2)
        tmp2_split

        if(any(tmp1_split %in% tmp2_split)) {  # If any number is repeated, add '9' before all scaffolds, so 10 becomes 910, and so on.
          df$CHR <- gsub(scaffoldIdentifier, paste0(scaffoldIdentifier, '9'), df$CHR)
          print('\n# ------ A number 9 was added at the beginning of each scaffold because numbering was repeated between Chr and scaffolds ------ #')
          print(head(unique(df$CHR[grepl(scaffoldIdentifier, df$CHR)])))
          df$CHR <- gsub(scaffoldIdentifier, '', df$CHR)  # Delete the scaffold name, just leave the number
        }
      }  # End, if there are scaffolds.

      cat('\n# ------ Deleting chrIdentifier ------ #\n')
      df$CHR <- gsub(chrIdentifier, '', df$CHR)  # Delete the chromosome name, just leave the number

      head(df)
      cat('\n# ------ Final numbering of Chr and scaffolds ------ #\n')
      cat(mixedsort(unique(df$CHR)))

      # Convert columns CHR and BP from character to numeric values.
      df$CHR <- as.numeric(df$CHR)
      df$BP <- as.numeric(df$BP)

      return(df)
    }

  } else {  # If CHROM and POS columns are not available
    cat('\n# ------ Columns SNP, or CHROM and POS are missing in the dataframe, cannot generate unique Chr numbering without those ------ #')
  }

  # For NWHer genome v.2.0
  #df$CHR <- gsub('unplaced_scaffold', 'unplaced_scaffold30', df$CHR)  # Trick to avoid unplaced scaffolds mess up the numbering of chromosomes.
  #df$CHR <- gsub('chr', '', df$CHR)
  #df$CHR <- gsub('unplaced_scaffold', '', df$CHR)
  #head(df)
  #unique(df$CHR)

  ## OLD way to do this, assuming all scaffolds/Chr are consecutive and do not repeat numbering.
  # Replace scaffold names (characters) by numbers.
  # Obtain a list of unique scaffold names.
  #uq <- unique(snp_split$CHR)

  # Create a vector with a sequence of numbers in ascending order from 1 to the number of unique scaffolds.
  #num <- seq(1, length(uq))

  # Safe in a TXT file the scaffold = chromosome number equivalence.
  #write.table(cbind.data.frame(uq, num), file = paste0('Scaffold-chromosome_equivalence_for_plotting_PC',i,'.txt'), row.names = FALSE, sep='\t', quote = FALSE)

  # Replace the original scaffold/contig name for a number (chromosome for plotting).
  #foo <- snp_split
  #for (i in seq_along(num)) foo$CHR[foo$CHR == uq[i]] <- num[i]

  #head(foo)

  # Convert columns CHR and BP from character to numeric values.
  #foo$CHR <- as.numeric(foo$CHR)
  #foo$BP <- as.numeric(foo$BP)

}
####################################### End of function #######################################

###############################################################################################
#'
#' Function to calculate dAF
#' Returns a dataframe with SNPs and dAF for a given paired comparison
#'
calculate.dAF <- function(poolData, comparison, group1, group2, na.rm) {
  ## For debugging.
  #comparison <- 'EI.II.III.IV_vs_LV.VI.VII.VIII'
  #group1 <- c('EI-F','EI-M','EII-F','EII-M','EIII-F','EIII-M','EIV-F','EIV-M')
  #group2 <- c('LV-F','LV-M','LVI-F','LVI-M','LVII-F','LVII-M','LVIII-F','LVIII-M')
  #comparison
  #group1
  #group2

  if(na.rm == TRUE) {
    ## Subset dataframe to keep only the columns (samples) corresponding to the two superpools, and remove NAs.
    poolData <- na.omit(poolData[ ,c(group1,group2)])
    head(poolData)
  }

  ## Calculate the row mean for Group 1.
  meanAF_group1 <- rowSums(as.data.frame(poolData[ ,c(group1)]), na.rm = TRUE)/length(group1)
  head(meanAF_group1)
  length(meanAF_group1)

  ## Calculate the row mean for Group 2.
  meanAF_group2 <- rowSums(as.data.frame(poolData[ ,c(group2)]), na.rm = TRUE)/length(group2)
  head(meanAF_group2)
  length(meanAF_group2)

  ## Calculate the absolute delta allele frequency (dAF) between the two superpools.
  dAF_df <- data.frame(SNP = rownames(poolData),
                       dAF = round(abs(meanAF_group1 - meanAF_group2),4), stringsAsFactors = FALSE)
  head(dAF_df)
  str(dAF_df)

  ## Clean up some memory space.
  rm(meanAF_group1)
  rm(meanAF_group2)

  return(dAF_df)

}
####################################### End of function #######################################

###############################################################################################
#' Function to find the closest gene to a SNP of interest,
#' inspired in code developed by Mats Petterson
#'
#' The function receives a path to or a df R object listing the SNPs of interest. It requires a GTF file listing gene annotations.
#' It returns a dataframe (for dAF analysis) or of list of unique gene names (for further GO analysis) describing the
#' gene_id, gene_name and distance in bp of each SNP to the xlosest gene up to a given window size
#'
#' Created by: Angela Fuentes-Pardo, Uppsala University
#' Date: 2020-07-25, updated on. 2020-12-03
#' e-mail: apfuentesp@gmail.com
#'
#' @param SNPfile Path to the file listing SNPs of interest, required CHROM-POS or SNP columns
#' @param SNPdf R data frame listing the SNPs of interest, required CHROM-POS or SNP columns
#' @param GTFfile Path to the GTF file containing gene annotations
#' @param windowSize Size of the maximun distance searched to indentify a neighbouring gene to a SNP locus
#' @param returnUniqueGeneList Set the function to return a list of unique genes (e.g. required for GO enrichment analysis)
#' @param returnDF Set the function to return the dataframe of distance ranges (default)
#'
getClosestGene2SNPs <- function(SNPfile = NULL, SNPdf = NULL, GTFfile = NULL, windowSize = NULL, returnUniqueGeneList = NULL, returnDF = TRUE) {

  # For debugging.
  # Original GTF file from Ensembl, without gene symbols
  #gtf_file <- 'Dropbox/PostDoc_UU/Projects/Herring/genome/Clupea_harengus.Ch_v2.0.2.96.gtf'
  # GTF file with gene symbols
  #GTFfile <- '~/Dropbox/PostDoc_UU/Projects/Herring/00_genome/other-annotations/Clupea_harengus.Ch_v2.0.2.100_geneSymbol.gtf'
  #SNPfile <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/Pool-seq/analysis/06-apply-Neff-calc-AF/NWHer_15_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMode1half.excIBC-S_NTS-F.AD.Neff.AF.txt'
  #SNPfile <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/06-apply-Neff-calc-AF/Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMode1half.72_Atl.AD.Neff.AF.txt'
  #windowSize <- 40e3  # Set window size, which is the maximum distance on bp to be considered, here we are testing 40Kb as Fan Han used this value
  #returnUniqueGeneList <- TRUE
  #returnDF <- TRUE

  # Load libraries
  library(data.table)
  library(GenomicRanges)
  library(Biostrings)
  library(VariantAnnotation)
  library(rtracklayer)
  library(GenomicFeatures)

  # -----------------------------------------------------------------------------------------------------------------------
  # 1. Load SNPs file and convert it to a GRanges object
  # -----------------------------------------------------------------------------------------------------------------------

  if(!is.null(SNPfile) & is.null(SNPdf)) {
    # Load file listing SNP loci of interest (expected SNP column, CHROM-POS, but this one could be created later based on CHROM and POS columns)
    lociNames_df <- fread(SNPfile, data.table = FALSE, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    #lociNames_df <- fread(list_inputFile_contrasts[grepl(i,list_inputFile_contrasts)], data.table = FALSE, header = TRUE, sep='\t', stringsAsFactors = FALSE)
  }
  if(is.null(SNPfile) & !is.null(SNPdf)) {
    lociNames_df <- SNPdf
  }

  head(lociNames_df)

  # Sort loci by genomic position, just in case
  #library(gtools)
  #sorted_loci <- mixedsort(lociNames_df[ ,1])  # Sort loci names
  #head(sorted_loci)
  #lociNames_df <- lociNames_df[match(sorted_loci,lociNames_df[,1]), ]

  if(!"SNP" %in% colnames(lociNames_df)) {  # If SNP column does not exist
    if("CHROM-POS" %in% colnames(lociNames_df)) {  # but CHROM-POS does
      colnames(lociNames_df) <- sub("CHROM-POS", "SNP", colnames(lociNames_df))  # rename column name
    }
    if("CHROM" %in% colnames(lociNames_df) & "POS" %in% colnames(lociNames_df)) {  # but CHROM and POS columns do
      lociNames_df$SNP <- paste0(lociNames_df$CHROM, lociNames_df$POS)
    }
  }

  head(lociNames_df)
  length(lociNames_df$SNP)
  #[1] 3609921

  # Remove duplicate loci, just in case
  lociNames_df <- data.frame(SNP = unique(lociNames_df$SNP), stringsAsFactors = FALSE)
  length(lociNames_df$SNP)
  #[1] 3609921
  str(lociNames_df)

  # Split SNP into chromosome and position
  library(tidyr)
  tmp <- separate(data = lociNames_df, col = SNP, into = c("chromosome", "start"), sep = "-")
  head(tmp)

  # Create dataframe with the required columns for GRanges object creation
  lociNames_df$chromosome <- tmp$chromosome
  lociNames_df$start <- tmp$start
  lociNames_df$end <- tmp$start
  rownames(lociNames_df) <- lociNames_df$SNP
  lociNames_df$SNP <- NULL
  head(lociNames_df)
  rm(tmp)

  # Convert dataframe into a GRanges object
  SNPgr <- makeGRangesFromDataFrame(lociNames_df, keep.extra.columns = TRUE)  # Last setting makes that associated values to loci are stored as metadata columns
  SNPgr
  class(SNPgr)

  # -----------------------------------------------------------------------------------------------------------------------
  # 2. Load gene annotations from a GTF/GFF file and convert it to a GRanges object
  # -----------------------------------------------------------------------------------------------------------------------

  ## Load the GTF file (with gene annotations) as a TxDb (or GRanges) object -----------
  # In TxDB format
  #CluharV202_gtf <- makeTxDbFromGFF('~/Dropbox/PostDoc_UU/Projects/Herring/genome/Clupea_harengus.Ch_v2.0.2.96.gtf',
  #                                  format = 'gtf')

  # In GRanges format
  GTFgr <- rtracklayer::import(GTFfile)
  GTFgr
  class(GTFgr)
  #GTFgr <- rtracklayer::import('~/Dropbox/PostDoc_UU/Projects/Herring/genome/Clupea_harengus.Ch_v2.0.2.96.gtf')

  # Add the string 'chr' to match the chromosome naming in the Herring genome
  seqlevels(GTFgr) <- paste0('chr',seqlevels(GTFgr))
  GTFgr
  unique(seqlevels(GTFgr))

  # Assign Ensembl gene_id to the "names" field
  names(GTFgr) <- GTFgr$gene_id  
  GTFgr

  # Subset GTF file to only retain the gene annotations
  GTFgr <- GTFgr[GTFgr$type == 'gene', ]
  GTFgr

  # -----------------------------------------------------------------------------------------------------------------------
  # 3. Find overlapping SNPs with GTF annotations up to a certain distance
  # -----------------------------------------------------------------------------------------------------------------------

  ## Find GRanges object listing the distance to the nearest gene regardless of how far it is.
  # The function distanceToNearest() returns the distance for each range in x
  # to its nearest neighbor in the subject. NOTE: When more than one hit has equal
  # distance, the first one in genomic order has precedence.
  # (https://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/GenomicRanges/html/nearest-methods.html)
  d <- distanceToNearest(x = SNPgr, subject = GTFgr, ignore.strand = TRUE)  
  d
  dim(mcols(d))

  # Subset Granges object to retain rows (loci) with a distance <= window size
  wd <- windowSize  # window size = maximum distance on bp to be considered t search for the "closest" gene
  #wd <- 0
  print(paste0('Window size: ', windowSize))

  d_sub <- d[(elementMetadata(d)[ ,"distance"] <= wd)]
  d_sub
  dim(mcols(d_sub))

  # Extract distances (always are absolute = positive value)  # how it used to be done
  #dists <- values(d)$distance
  #dists

  # Get the index of the subject (gene annotations) and query (loci of interest) hits that passed the distance requirement
  index_subject <- subjectHits(d_sub)
  head(index_subject)
  index_query <- queryHits(d_sub)
  head(index_query)

  #index_subject <- subjectHits(d)[dists <= wd]  # how it used to be done
  #head(index_subject)
  #index_query <- queryHits(d)[dists <= wd]
  #head(index_query)

  # Assign the gene name fulfilling the distance criteria to the correspond SNP ----------
  #mcols(SNPgr)$gene_id[index_query] <- mcols(GTFgr)$gene_id[index_subject]

  # Create the new metadata columns - fill with NA (this helps to avoid problems when there are no hits)
  SNPgr$gene_id <- rep(NA, nrow(mcols(SNPgr)))
  SNPgr$gene_name <- rep(NA, nrow(mcols(SNPgr)))
  SNPgr$distance <- rep(NA, nrow(mcols(SNPgr)))
  SNPgr

  # Assign gene_ids (stored as "names") to the correspondent SNP
  SNPgr$gene_id[index_query] <- names(GTFgr[index_subject])  
  SNPgr$gene_name[index_query] <- GTFgr$gene_name[index_subject]

  SNPgr$gene_idname <- paste0(SNPgr$gene_id, " ", SNPgr$gene_name)  # make a combined columns, useful when no gene_name is available
  SNPgr$distance[index_query] <- mcols(d_sub)$distance #****
  SNPgr

  dim(mcols(SNPgr))
  #[1] 3609921       3

  if(isTRUE(returnUniqueGeneList)) {

    # Extract gene IDs and remove duplicates
    gene_list <- SNPgr$gene_id
    length(gene_list)
    #[1] 3609921
    gene_list <- gene_list[!is.na(gene_list)]
    length(gene_list)
    #[1] 2619
    #[1] 1581
    gene_list <- unique(gene_list)
    length(gene_list)
    #[1] 472
    #[1] 182
    cat("\nExecuting returnUniqueGeneList...")
    return(gene_list)
  }

  if(isTRUE(returnDF)) {

    SNP_df <- as.data.frame(SNPgr, stringsAsFactors = FALSE)
    SNP_df <- data.frame(SNP = rownames(SNP_df), SNP_df, stringsAsFactors = FALSE)
    dim(SNP_df)
    head(SNP_df)
    print("\nExecuting returnDF...")
    print("# ---------------- gene annotations added ---------------- #")
    return(SNP_df)
  }

  # Save output as a file
  #write.table(gene_list, file = paste0('gene_list-NWHer_15pools_',outputName_list[i],'.txt'), col.names = FALSE, row.names = FALSE, quote = FALSE)

}
####################################### End of function #######################################

###############################################################################################
#'                                                                   
#' Function that ranks snpEff effects and adds them to a dAF file (previously named processAnnotatedSNPs)
#' Developed by: Mats Pettersson, mats.pettersson@imbim.uu.se
#' Modifications by: Angela P. Fuentes-Pardo, apfuentesp@gmail.com
#' Uppsala University                  
#' Date: 2020-12-03, modified: 2021-01-05                          
#'
addSnpEffAnn <- function(ann_file = NULL, daf_df = NULL) {

  # For debugging
  #ann_file <- snpEff_ANN_file
  #daf_df <- summary_dAF_df

  ## This script assumes a file with the ANN field annotated with snpEff was extracted from a VCF file, using GATK-VariantsToTable function

  ## If the ANN_split file does not exist, create it
  if(!file.exists(paste0(ann_file,"_split.txt"))) {
    ## Create the awk script required
    awkfile <- file('snpEFF_awk.txt')
    writeLines(c('BEGIN{FS=","; ORS="\\t"}
  {
	match($0, /[a-z_]+[0-9]+\\t[0-9]+\\t[ATGC]\\t[ATGC]/);
	print substr($0,RSTART, RLENGTH)
  }
  {
	for(i=1; i<=NF; i++)
	{
		match($i, /[ATCG][|][35A-Za-z_]+/);
		print substr($i,RSTART, RLENGTH);
	}
	printf "%s", "\\n"
  }'), awkfile)
    close(awkfile)
    #system('cat snpEFF_awk.txt')  # verify the accesory script was created

    ## Save the filename and path
    library(tools)
    filename <- file_path_sans_ext(basename(ann_file))
    dirpath <- dirname(ann_file)
    filename
    dirpath

    ## Run the AWK script to to parse the ANN field file
    system(paste0("cat ",ann_file," | awk -f snpEFF_awk.txt > ",ann_file,"_split.txt"))
    #system('cat NWHer_15_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMode1half.excIBC-S_NTS-F.snpEff.ANN | awk -f snpEFF_awk.txt > NWHer_15_pools.DPMode1half.snpEff.ANN_split.txt')
    system('rm snpEFF_awk.txt')  # delete the awk script

  }

  ## Read in the split ANN file as character lines and create a dataframe with the format required by the addSnpEffAnn() function
  snpEff_ANN <- readLines(file(paste0(ann_file,"_split.txt")))
  snpEff_ANN <- snpEff_ANN[-1]

  ## Extract REF and ALT columns
  library(stringr)
  library(tidyr)
  #test <- snpEff_ANN[1:5]
  #str_view(test, "([ATCG].?)\t([ATCG])"))
  #str_extract(test, "([ATCG].?)\t([ATCG])"))
  df <- data.frame(COL = str_extract(snpEff_ANN, "([ATCG].?)\t([ATCG])"), stringsAsFactors = FALSE)
  df <- tidyr::separate(data = df, col = COL, into = c("REF", "ALT"), sep = "\t")
  head(df)
  snpEff_ANN_df <- data.frame(SNP = sub("([a-z0-9]+)\t([0-9]+).+", "\\1-\\2", snpEff_ANN),
                              REF = df$REF,
                              ALT = df$ALT,
                              raw = snpEff_ANN,
                              stringsAsFactors = FALSE)
  head(snpEff_ANN_df)
  #SNP REF ALT                                                                                                                      raw
  #1  chr1-757   C   T  chr1\t757\tC\tT\tT|downstream_gene_variant\tT|downstream_gene_variant\tT|downstream_gene_variant\tT|intergenic_region\t
  #2 chr1-1383   T   C chr1\t1383\tT\tC\tC|downstream_gene_variant\tC|downstream_gene_variant\tC|downstream_gene_variant\tC|intergenic_region\t
  #3 chr1-1403   T   A chr1\t1403\tT\tA\tA|downstream_gene_variant\tA|downstream_gene_variant\tA|downstream_gene_variant\tA|intergenic_region\t

  # Sort snpEff_ANN_df

  ## If SNP order is the same between dataframes
  if(all(snpEff_ANN_df$SNP == daf_df$SNP)) {

    # Bind the two dfs
    ann_df <- cbind(daf_df, snpEff_ANN_df)
    head(ann_df)

    #require(GenomicRanges)
    #ann_GR <- GRanges(seqnames = ann_df$CHROM, ranges = IRanges(start = ann_df$POS, end = ann_df$POS))
    #in_selection_reg <- findOverlaps(ann_GR, selection_reg_GR)
    # Annotate if a region is within known selection regions
    #ann_df[ ,"Selection_reg"] <- FALSE
    #ann_df[in_selection_reg@from,"Selection_reg"] <- TRUE

    # Assign snpEff effect annotations
    ann_df[ ,"NonSyn"] <- grepl("missense_variant", ann_df$raw)
    ann_df[ ,"Syn"] <- grepl("synonymous_variant", ann_df$raw)
    ann_df[ ,"Intron"] <- grepl("intron_variant", ann_df$raw)
    ann_df[ ,"FiveUTR"] <- grepl("5_prime_UTR_variant", ann_df$raw)
    ann_df[ ,"ThreeUTR"] <- grepl("3_prime_UTR_variant", ann_df$raw)
    ann_df[ ,"Downstream"] <- grepl("downstream_gene_variant", ann_df$raw)
    ann_df[ ,"Upstream"] <- grepl("upstream_gene_variant", ann_df$raw)
    ann_df[ ,"Intergenic"] <- grepl("intergenic_region", ann_df$raw)

    # ANN ranking (nr = non-redundant)
    ann_df[,"Downstream_nr"] <- ann_df[,"Downstream"] & !(ann_df[,"NonSyn"] | ann_df[,"Syn"] | ann_df[,"Intron"] | ann_df[,"FiveUTR"] | ann_df[,"ThreeUTR"])
    ann_df[,"Upstream_nr"] <- ann_df[,"Upstream"] & !(ann_df[,"NonSyn"] | ann_df[,"Syn"] | ann_df[,"Intron"] | ann_df[,"FiveUTR"] | ann_df[,"ThreeUTR"])
    ann_df[,"Intergenic_nr"] <- ann_df[,"Intergenic"] & !(ann_df[,"Downstream"] | ann_df[,"Upstream"] | ann_df[,"NonSyn"] | ann_df[,"Syn"] | ann_df[,"Intron"] | ann_df[,"FiveUTR"] | ann_df[,"ThreeUTR"])
    ann_df[,"Intron_nr"] <- ann_df[,"Intron"] & !(ann_df[,"NonSyn"] | ann_df[,"Syn"] | ann_df[,"FiveUTR"] | ann_df[,"ThreeUTR"])

    head(ann_df)
    ann_df$raw <- NULL
    print("# ---------------- snpEff annotations added ---------------- #")
    return(invisible(ann_df))

  } else {
    print("SNP mismatch between dfs, needs correction")
  }
}
############################################# END #############################################


## --------------------------------------------------------------------------------------------
## Main body of code
## --------------------------------------------------------------------------------------------

## Clean environment space.
#rm(list = ls())

library(optparse) # https://github.com/trevorld/r-optparse, https://cran.r-project.org/web/packages/optparse/optparse.pdf

## Import run parameters from the command line.
option_list <- list(
  make_option('--workDir', type='character', default = NULL, help='Set working directory.'),
  make_option('--AFfile', type='character', default = NULL, help='Path to input file of population allele frequencies (single entry for each sample, respect to major allele). Format: samples in rows and loci in columns.'),
  make_option('--contrastsFile', type='character', default = NULL, help='Path to file listing paired contrasts to perform. Items separated by spaces.'),
  make_option('--chrIdentifier', type='character', default = NULL, help='Chromosome string identifier (e.g. Chr or chr).'),
  make_option('--scaffIdentifier', type='character', default= NULL, help='Scaffold string identifier (e.g. scaffold or unplaced_scaffold).'),
  make_option('--scaffSuffix', type='character', default = NULL, help='Scaffold suffix (e.g. _grp_1).'),
  make_option('--delimiterChrScaffNumber', type='character', default= NULL, help='Character separator of Chromosome or scaffolds and their number (e.g. - or _ in Chr-1, or Chr_1).'),
  make_option('--delimiterChrPos', type='character', default = NULL, help='Character separator of Chromosome or scaffolds and the SNP position (e.g. - or _ in Chr1-298, or Chr_1-298'),
  make_option('--convertFormat', type='logical', default = FALSE, help='Choose whether the format of the file listing pair contrasts should be converted from wide to compact format. Options: TRUE or FALSE'),
  make_option('--removeNAs', type='logical', default = FALSE, help='Choose whether the loci with missing calues (NAs) should be removed from the pool-seq dataset. Options: TRUE or FALSE'),
  make_option('--addGeneAnn', type='logical', default = FALSE, help='Choose whether gene annotations should be added to the dAF file. Options: TRUE or FALSE'),
  make_option('--gtfFile', type='character', default = NULL, help='Path to the GTF file listing gene annotations.'),
  make_option('--winSize', type='numeric', default = NULL, help='Set a window size to find the closest gene to a given SNP'),
  make_option('--addsnpEffAnn', type='logical', default = FALSE, help='Choose whether snpEff annotations should be added to the dAF file. Options: TRUE or FALSE'),
  make_option('--snpEffAnnFile', type='character', default = NULL, help='Path to the snpEff ANN file after extracted from the annotated VCF file with snpEff.'),
  make_option('--addSelectionRegions', type='logical', default = FALSE, help='Choose whether snpEff annotations should be added to the dAF file. Options: TRUE or FALSE'),
  make_option('--selectionRegionsFile', type='character', default = NULL, help='Path to a GRanges file listing known selection regions.')
)

opts <- parse_args(OptionParser(option_list = option_list))
opts

## Assign parameters to local R variables.
library(stringr)
working_dir <- opts$workDir
input_file <- opts$AFfile
contrasts_file <- opts$contrastsFile
# Set string identifiers and separators of chromosomes, scaffolds, and SNP positions.
if(!is.null(opts$chrIdentifier)) {
  chr_identifier <- opts$chrIdentifier
} else { chr_identifier <- NULL }
if(!is.null(opts$scaffIdentifier)) {
  scaff_identifier <- opts$scaffIdentifier
} else { scaff_identifier <- NULL }
if(!is.null(opts$scaffSuffix)) {
  scaff_suffix <- opts$scaffSuffix
} else { scaff_suffix <- NULL }
if(!is.null(opts$delimiterChrScaffNumber)) {
  delimiter_ChrScaffNumber <- opts$delimiterChrScaffNumber
} else { delimiter_ChrScaffNumber <- NULL }
if(!is.null(opts$delimiterChrPos)) {
  delimiter_ChrPos <- opts$delimiterChrPos
} else { delimiter_ChrPos <- NULL }
if(!is.null(opts$convertFormat)) {
  convert_format <- opts$convertFormat
} else { convert_format <- FALSE }
if(!is.null(opts$removeNAs)) {
  remove_NAs <- opts$removeNAs
} else { remove_NAs <- FALSE }
if(!is.null(opts$addGeneAnn)) {
  add_gene_ann <- opts$addGeneAnn
} else { add_gene_ann <- FALSE }
gtf_file <- opts$gtfFile
if(!is.null(opts$winSize)) {
  window_size <- opts$winSize
} else { window_size <- NULL }
if(!is.null(opts$addsnpEffAnn)) {
  add_snpEff_ann <- opts$addsnpEffAnn
} else { add_snpEff_ann <- FALSE }
if(!is.null(opts$snpEffAnnFile)) {
  snpEff_ANN_file <- opts$snpEffAnnFile
} else { snpEff_ANN_file <- FALSE }
if(!is.null(opts$addSelectionRegions)) {
  add_selection_regions <- opts$addSelectionRegions
} else { add_selection_regions <- FALSE }
if(!is.null(opts$selectionRegionsFile)) {
  selectionSNPs_file <- opts$selectionRegionsFile
} else { selectionSNPs_file <- FALSE }


## For debugging.
# For Her genome v.2.0
working_dir <- '~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/08-2-outlier-detection_dAF'
input_file <- '~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/06-apply-Neff-calc-AF/OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.AD.Neff.AF.txt'
contrasts_file <- '~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/08-2-outlier-detection_dAF/List_paired_comparisons_OTT_2020-09-01.csv'
chr_identifier <- 'chr'
scaff_identifier <- 'scaffold'
scaff_suffix <- NULL
delimiter_ChrScaffNumber <- '_'
delimiter_ChrPos <- '-'
convert_format <- FALSE  # from ?? to ??
remove_NAs <- FALSE


print('# ------------- Set R objects from BASH command line: ------------- #')
working_dir
input_file
contrasts_file
chr_identifier
scaff_identifier
scaff_suffix
delimiter_ChrScaffNumber
delimiter_ChrPos
convert_format
remove_NAs

## Set working directory.
setwd(working_dir)

## ------------------------------------------------------------------
## Data loading and preprocessing
## ------------------------------------------------------------------

## Load the pool allele frequency data.

## Save filename.
library(tools)
filename <- file_path_sans_ext(basename(input_file))
filename

library(data.table)
poolData <- fread(input_file, data.table = FALSE, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
#head(poolData)  # Explore first lines of the dataframe.
row.names(poolData) <- poolData[ ,1]  # Assign first column with pool names to rownames.
#poolData <- poolData[, -1]  # Remove such first column.
head(poolData)  # Verify modification.
dim(poolData)
str(poolData)
class(poolData)

## Fix column names (optional if naming starts with a number, an X will be inserted at the beginning).
#colnames(poolData)
#tidy_colnames <- make.names(colnames(poolData), unique = TRUE)
#colnames(poolData) <- tidy_colnames
#colnames(poolData)

## Sort poolData df in ascending CHR and BP order.
#poolData <- poolData[order(poolData$CHROM.POS),]  
#head(poolData)

## Read in file listing paired comparisons (super pools).
#comparisonPairs_df <- read.csv(contrasts_file, sep = ';', stringsAsFactors = FALSE)
if(convert_format == TRUE) {
  source('~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/code/utility-code/convertFormat.pairContrasts_v2.R')
  comparisonPairs_df <- convertFormat.pairContrasts(contrasts_file)  # Function accepts TXT and CSV files

} else {

  library(tools)
  fileExt <- file_ext(contrasts_file)  # Extract the extension of the input file

  if(fileExt == 'csv') {
    comparisonPairs_df <- read.table(contrasts_file, header = TRUE, sep = ';', fill = TRUE, stringsAsFactors = FALSE)  # For a csv file
    #head(comparisonPairs_df)
  }
  if(fileExt == 'txt') {
    comparisonPairs_df <- read.table(contrasts_file, header = TRUE, sep = '\t', fill = TRUE, stringsAsFactors = FALSE)  # For a TXT file
    #head(comparisonPairs_df)
  }

  # Remove columns with NAs
  comparisonPairs_df <- comparisonPairs_df[ ,colSums(is.na(comparisonPairs_df)) != nrow(comparisonPairs_df)]
  head(comparisonPairs_df)
  str(comparisonPairs_df)
}

print('# ---------------- File listing the paired comparisons ---------------- #')
comparisonPairs_df

## ------------------------------------------------------------------
## Calculate dAF for each paired comparison
## ------------------------------------------------------------------

## Create a dataframe that will store the dAF for each pairwise comparison.
summary_dAF_df <- data.frame(SNP = rownames(poolData), stringsAsFactors = FALSE)
head(summary_dAF_df)

## Replace scaffold names (characters) by numbers.
#source('/home/afuentes/project/afuentes/code/utility-code/createCHR.uniqueNumbers.R')
summary_dAF_df <- createCHR.uniqueNumbers(df = summary_dAF_df,
                                          chrIdentifier = chr_identifier,
                                          scaffoldIdentifier = scaff_identifier,
                                          scaffSuffix = scaff_suffix,
                                          delimiterChrScaffNumber = delimiter_ChrScaffNumber,
                                          delimiterChrPos = delimiter_ChrPos)
cat('\n')
head(summary_dAF_df)
str(summary_dAF_df)

## For loop that goes row by row in the df to compute dAF for each pairwise comparison.
for(i in 1:length(rownames(comparisonPairs_df))) {
  ## Set group pairs for comparison.
  #i = 1
  comparison <- strsplit(comparisonPairs_df[i,"Comparison"], "\\s+")[[1]]  # strsplit() split character string into character vector by a separator, space(s) in this case
  group1 <- strsplit(comparisonPairs_df[i,"Group1"], "\\s+")[[1]]
  group2 <- strsplit(comparisonPairs_df[i,"Group2"], "\\s+")[[1]]

  ## Print current comparison and groups.
  print(comparison)
  print(group1)
  print(group2)

  ## Calculate dAF for the current pair group.
  dAF_tmp <- calculate.dAF(poolData = poolData,
                           comparison = comparison,
                           group1 = group1,
                           group2 = group2,
                           na.rm = remove_NAs)
  print(head(dAF_tmp))

  ## Add the dAF of the current pair comparison to the growing summary table of dAF.
  summary_dAF_df[ ,paste0('dAF_',comparison)] <- dAF_tmp$dAF[match(summary_dAF_df$SNP, dAF_tmp$SNP, nomatch = NA)] # First larger df, second smaller df.
  head(summary_dAF_df)
  summary(summary_dAF_df[ ,paste0('dAF_',comparison)])
  #summary_dAF_df[is.na(summary_dAF_df$dAF_WIreland_1b_vs_WIreland_1a.2), ]  # Explore which rows are NA = no match.
  #dAF_tmp[dAF_tmp$SNP == 'tarseq_5-3819508', ]  # abscent
  #summary_dAF_df[summary_dAF_df$SNP == 'tarseq_5-3819508', ] # present
  #Example when loci match.
  #dAF_tmp[dAF_tmp$SNP == 'tarseq_85-6299535', ]  # present
  #summary_dAF_df[summary_dAF_df$SNP == 'tarseq_85-6299535', ]  # present

  ## Generate Manhattan plot.
  #make.Manhattanplot(summary_dAF_df = summary_dAF_df,
  #                   filename = filename,
  #                   comparison = comparison)

} # End of for loop

print("")
print(head(summary_dAF_df))

## Add gene, snpEff, and selection region annotations if required.

if(isTRUE(add_gene_ann)) { # --------------- Add gene annotations to the dAF df --------------- #

  #source('/proj/snic2020-2-19/private/herring/users/angela/Pool-seq/code/utility-code/getClosestGene2SNPs_v2.R')
  #source('~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/code/utility-code/getClosestGene2SNPs_v2.R')
  SNPnames_df <- data.frame(SNP = summary_dAF_df$SNP, stringsAsFactors = FALSE)
  head(SNPnames_df)
  closestGeneAnn_df <- getClosestGene2SNPs(SNPdf = SNPnames_df, GTFfile = gtf_file, windowSize = window_size, returnDF = TRUE)
  head(closestGeneAnn_df)

  if(all.equal(SNPnames_df$SNP, closestGeneAnn_df$SNP)) {
    summary_dAF_df$gene_id <- closestGeneAnn_df$gene_id
    summary_dAF_df$gene_name <- closestGeneAnn_df$gene_name
    summary_dAF_df$gene_idname <- closestGeneAnn_df$gene_idname
    summary_dAF_df$distance <- closestGeneAnn_df$distance
  }
}

if(isTRUE(add_snpEff_ann)) { # --------------- Add snpEff effect predictions to the dAF df --------------- #

  ## the function addSnpEffAnn() extracts the more informative effects, right? <<<
  summary_dAF_df <- addSnpEffAnn(ann_file = snpEff_ANN_file,
                                 daf_df = summary_dAF_df)
} # End add_snpEff_ann
# Ignore all those unused connection warnings

if(isTRUE(add_selection_regions)) { # --------------- Annotate if SNPs are within known selection regions --------------- #

  summary_dAF_df <- annotateSelectionSNPs(daf_df = summary_dAF_df,
                                          selectionSNPs_GR = selectionSNPs_file)
} # End add_selection_regions

## Examine first lines of the final dAF file ---------
head(summary_dAF_df)

## Save the dataframe with all dAF values in a TXT file.
library(tools)
contrasts_fileName <- file_path_sans_ext(basename(contrasts_file))
filename
contrasts_fileName

write.table(summary_dAF_df, file = paste0('Summary_dAF.',contrasts_fileName,'.',filename,'.txt'),
            row.names = FALSE, sep = '\t', quote = FALSE)

```

**Note** that this script offers many possibilities (like annotating the closest gene or highlight of snpEff predictions), but for now, I only calculated the dAF since the liftover of gene annotations is pending.


## Obtain Manhattan plots of dAF

- dAF files:
```bash
cd /proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/08-outlier-loci-detection/dAF/

nano dAF_files.txt

/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/08-outlier-loci-detection/dAF/Summary_dAF.List_paired_comparisons_OTT_2020-09-01.OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.AD.Neff.AF.txt
```
- Output dir:
```bash
nano outputdir_dAF.txt

20210331
```
Script `08-2-create-ManhatPlots-dAF-allChr.R`:
```R
## -----------------------------------------------------------------------------------------
## Main body of the code 
## -----------------------------------------------------------------------------------------

## Clean environment space.
rm(list = ls())

## Receive arguments from the bash command.
args <- commandArgs(trailingOnly = TRUE)
cat(args, sep='\n')

## Assign these arguments to variables.
working_dir <- args[1]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/2020-03-25
dAF_file <- args[2]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/2020-03-25/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.AD.Neff.AF.txt
#AF_file <- args[3]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/2020-03-25/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.strictDP.exc-6a-6b-7.mono.AD.Neff.AF.txt
outputName <- args[3]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/2020-03-25/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.strictDP.exc-7.mono.AD.Neff.AF.txt
fai_sorted_file <- args[4]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/2020-03-25/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.strictDP.allSamples.mono.AD.Neff.AF.txt
rollMeanWSwide <- args[5] # /proj/snic2020-16-14/private/HorseMackerel/data/00-genome/v.2.0/fTraTra1_1.curated_primary.20200310_sorted.fa.fai
#rollMeanWSzoom <- args[6]  # /proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/Scaffolds_of_interest_HOM_dAF.txt
#sampleOrder <- args[6]
selection_regions <- args[6]
selection_regions_type <- args[7]
# Workaround since these items are not available (2021-03-31)
#selection_regions <- NULL
#selection_regions_type <- NULL

## For debugging (everything).
#working_dir <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/08-2-outlier-detection_dAF/plots'
#dAF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/08-2-outlier-detection_dAF/Summary_dAF.Herring_popGen2_contrasts_20201128_BrackishMigratory.Formatted.Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.72_Atl.AD.Neff.AF.txt'
#dAF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/Population-contrasts/comparison_chisqr_dAFplots_fig4/Summary_ChiSqr.test_Fig4.Formatted.Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMean1SD.AD.Neff.txt'
#dAF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/Population-contrasts/comparison_chisqr_dAFplots_fig4/Summary_ChiSqr.test_Fig3_4.Formatted.Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.AD.Neff.txt'
#dAF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/Population-contrasts/comparison_chisqr_dAFplots_fig4/Summary_dAF.test_Fig3_4.Formatted.Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMode1half.72_Atl.AD.Neff.AF.txt'
#dAF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/Population-contrasts/comparison_chisqr_dAFplots_fig4/Summary_dAF.test_Fig3_4.Formatted.Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.72_Atl.AD.Neff.AF.txt'
#dAF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/08-2-outlier-detection_dAF/Summary_dAF.Herring_popGen2_contrasts_20201128_Salinity.Formatted.Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.72_Atl.AD.Neff.AF.txt'
#AF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/06-apply-Neff-calc-AF/Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPMode1half.AD.Neff.AF.txt'
#AF_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/06-apply-Neff-calc-AF/Her_73_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.72_Atl.AD.Neff.AF.txt'
#outputName <- 'Her_72_pools_SNPs_DPq5-99'
#outputName <- 'Her_72_pools_SNPs_DPMode'
#outputName <- 'Her_72_pools_SNPs_DPMean'
#fai_sorted_file <- '~/Dropbox/PostDoc_UU/Projects/Herring/00_genome/Ch_v2.0.2.fasta.fai'
#regionsInterest_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/2020-03-25/Scaffolds_of_interest_HOM_dAF_2020-03-25.txt'
#regionsInterest_file <- '~/Dropbox/PhD_Dal/2019-11_Reanalysis/analysis/08-2-outlier-detection_dAF/Informative_regions_dAF_2020-06-22.txt'
#rollMeanWSwide <- 100
#rollMeanWSzoom <- 100
#sampleOrder  <- c('A_Kalix_Baltic_Spring','PB4_Hudiksvall_Baltic_Spring','PB5_Galve_Baltic_Spring','PB1_HastKar_Baltic_Spring','B_Vaxholm_Baltic_Spring','HGS1_Riga_Baltic_Spring','HGS2_Riga_Baltic_Spring','G_Gamleby_Baltic_Spring','PB11_Kalmar_Baltic_Spring','PB12_Karlskrona_Baltic_Spring','PN3_CentralBaltic_Baltic_Spring','HGS40_Gdansk_Baltic_Spring','PB6_Galve_Baltic_Summer','PB7_Galve_Baltic_Autumn','HGS3_Riga_Baltic_Autumn','HGS4_Riga_Baltic_Autumn','HGS12_BornholmBasin_Baltic_Autumn','TysklandS18_Germany_Baltic','HGS71_Rugen_Baltic_Spring','HGS72_Rugen_Baltic_Spring','HGS6_Schlei_Baltic_Spring','H_Fehmarn_Baltic_Autumn','HGS5_Schlei_Baltic_Autumn','J_Traslovslage_Baltic_Spring','PB9_Kattegat_Atlantic_Spring','HGS8_KattegatNorth_Atlantic_Spring','PB10_Skagerrak_Atlantic_Spring','O_Hamburgsund_Atlantic_Spring','HGS24_Landvik_Atlantic_Spring','HGS45_Strandfjorden_Atlantic_Spring','LandvikS17_Norway_Baltic_Spring','HGS44_Hovag_Atlantic_Spring','HGS11_RingkobingFjord_NorthSea_Spring','HGS15_NSSH_Atlantic_Spring','Q_Norway_Atlantic_Atlantic_Spring','HGS42_BergenSp_Atlantic_Spring','HGS43_BergenAut_Atlantic_Autumn','N_NorthSea_Atlantic_Autumn','HGS27_Gloppen_Atlantic_Spring','HGS26_Lusterfjorden_Atlantic_Spring','HGS25_Lindas_Atlantic_Spring','PB2_Iceland_Atlantic_Spring','HGS41_FaroeIsland_Atlantic_Autumn','HGS46_Minch_Atlantic_Spring','HGS20_CapeWrath_Atlantic_Spring','HGS23_Clyde_Atlantic_Spring','HGS36_MilfordHavenApr_Atlantic_Spring','HGS35_MilfordHavenFeb_Atlantic_Spring','HGS39_Tenby_Atlantic_Spring','HGS16_Orkney_NorthSea_Autumn','HGS32_NWCapeWrath_Atlantic_Autumn','HGS22_CapeWrath_Atlantic_Autumn','HGS21_Hebrides_Atlantic_Mixed','HGS33_Buchan_Atlantic_Autumn','HGS28_Banks_Atlantic_Autumn','HGS17_IsleOfMan_IrishSea_Autumn','HGS34_DouglasBank_Atlantic_Autumn','HGS19_TeelinBay_Atlantic_Winter','HGS30_BrucklessBay_Atlantic_Winter','HGS29_TeelinBayLateType_Atlantic_Winter','HGS31_CelticSea_Atlantic_Winter','HGS18_CelticSea_Atlantic_AutumnWinter','HGS38_Minehead_Atlantic_Winter','HGS37_Clovelly_Atlantic_Winter','HGS10_Downs_EnglishChannel_Winter','HGS9_Greenland_Atlantic_Spring','DalInB_Atlantic_Spring','DalNsS_Atlantic_Spring','DalFB_Atlantic_Spring','DalBoB_Atlantic_Autumn','DalNsF_Atlantic_Autumn','DalGeB_Atlantic_Autumn','PB8_Pacific_Pacific_Spring')
#selection_regions <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/08-2-outlier-detection_dAF/Selection-regions/SNPs_selection_regions_SalinitySpawning4INV.RData'
#selection_regions <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/08-2-outlier-detection_dAF/Selection-regions/joint_chi_dfs_selection_regions_Mats.RData'
#selection_regions <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/08-2-outlier-detection_dAF/Selection-regions/selection_regions_GRanges_Mats.RData'
#selection_regions <- '~/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/analysis/08-2-outlier-detection_dAF/Selection-regions/SNPs_selection_regions_SalinitySpawning4INV_uncorrectedPvalues.Rdata'
#selection_regions_type <- 'loci' # 'ranges'

## For debugging (heatmap).
#working_dir <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF'
#dAF_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.Neff.AF.500K.txt'
#input_file2 <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.Neff.AF.500K.txt'
#input_file3 <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.Neff.AF.500K.txt'
#input_file4 <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/Summary_dAF.List_paired_comparisons_HOM_2020-02-21.HOM_12_pools.SNPs.hf.DP10-885.GQ10.mono.miss20.maf0.05.Neff.AF.500K.txt'
#fai_sorted_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/data/00-genome/v.2.0/fTraTra1_1.curated_primary.20200310.fa.fai'
#fai_sorted_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/data/00-genome/v.1.0/fTraTra1.PB.asm1.purge1.scaff1.fa.fai'

#regionsInterest_file <- '~/Dropbox/PostDoc_UU/Projects/HorseMackerel/Pool-seq/analysis/08-2-outlier-detection_dAF/2020-03-25/Scaffolds_of_interest_HOM_dAF.txt'

print('# ------------- Set R objects from BASH command line: ------------- #')
working_dir
dAF_file
#AF_file
#regionsInterest_file
fai_sorted_file
#rollMeanWSwide
#rollMeanWSzoom
#sampleOrder

## Set working directory.
setwd(working_dir)

# Chromosome order wanted for plotting.
# For OTT genome v.2.0.
library(data.table)
fai <- fread(fai_sorted_file, data.table = FALSE, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
CHROM_order <- fai[,1]

## -----------------------------------------------------------------------------------------
## 1. Generate a Manhattan plot for ALL chromosomes
## -----------------------------------------------------------------------------------------

source('/proj/snic2020-2-19/private/herring/users/angela/Pool-seq/code/utility-code/make.ManhattanPlot_v2.R')
#source('/Users/angfu103/Dropbox/PostDoc_UU/Projects/Herring/01_PopGen2/code/utility-code/make.ManhattanPlot_v2.R')
library(data.table)

dAF_df <- fread(dAF_file, data.table = FALSE, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
paired_comparisons <- grep("dAF_", names(dAF_df), value = TRUE)
#paired_comparisons <- grep("ChiSqr_", names(dAF_df), value = TRUE)
paired_comparisons

# Subset
#paired_comparisons <- c("dAF_WIreland_1a_vs_Others_3.4.8.5a.5b.9","dAF_NorthSea_3.4_vs_Others_1a.8.5a.5b.9","dAF_Spain_8_vs_Others_1a.3.4.5a.5b.9",
#                        "dAF_NPortugal_5a_vs_Others_1a.3.4.8.5b.9","dAF_SPortugal_5b_vs_Others_1a.3.4.8.5a.9","dAF_Mediterranean_9_vs_Others_1a.3.4.8.5a.5b",
#                        "dAF_Africa_7_vs_Others_1a.3.4.8.5a.5b.9","dAF_SPortugAfrica_5b.7_vs_Others_1a.3.4.8.5a.9","dAF_WIreland_1b_vs_Others_1a.3.4.8.5a.5b.9")

# Find which SNPs are within known regions of selection associated with spawning time, salinity adaptation, and known inversions
if(selection_regions_type == 'ranges') {
  
  load(selection_regions) # “sal_p20_regions_GR”  “spwn_p20_regions_GR”
  head(sal_p20_regions_GR)  # Salinity genomic ranges
  head(spwn_p20_regions_GR)  # Spawning genomic ranges
  invlocalPops_GR <- GRanges(seqnames = c('chr6', 'chr12', 'chr17', 'chr23'), 
                             ranges = IRanges(start = c(22.2*1e+6, 17.81e+6, 25.81e+6, 16.31e+6), end = c(24.81e+6, 25.61e+6, 27.51e+6, 17.51e+6)))
  invlocalPops_GR
  #Chr6:22.2–24.8Mb
  #Chr12:17.8–25.6Mb
  #Chr17:25.8–27.5Mb
  #Chr23:16.3–17.5Mb
  
  head(dAF_df)
  
  library(GenomicRanges)
  SNPs_GR <- GRanges(seqnames = dAF_df$CHROM, ranges = IRanges(start = dAF_df$POS, end = dAF_df$POS))
  in_Salinity_region <- findOverlaps(SNPs_GR, sal_p20_regions_GR)
  in_Spawning_region <- findOverlaps(SNPs_GR, spwn_p20_regions_GR)
  in_invLocalPop <- findOverlaps(SNPs_GR, invlocalPops_GR)
  #in_Salinity_region_flank <- findOverlaps(SNPs_GR, sal_p20_regions_GR, maxgap = 1e4)  # add 1kb flank
  #in_Spawning_region_flank <- findOverlaps(SNPs_GR, spwn_p20_regions_GR, maxgap = 1e4)
  
  dAF_df[ ,"Salinity_region"] <- FALSE
  dAF_df[in_Salinity_region@from,"Salinity_region"] <- TRUE
  dAF_df[ ,"Spawning_region"] <- FALSE
  dAF_df[in_Spawning_region@from,"Spawning_region"] <- TRUE
  dAF_df[ ,"InvLocalPop_region"] <- FALSE
  dAF_df[in_invLocalPop@from,"InvLocalPop_region"] <- TRUE
  
}

if(selection_regions_type == 'loci') {
  
  library(GenomicRanges)
  
  load(selection_regions)  #  "AvS_join_pvalue10" "BvA_join_pvalue10" "sv_pvalue10"
  head(BvA_join_pvalue10)  # Salinity genomic ranges
  head(AvS_join_pvalue10)  # Spawning genomic ranges
  head(sv_pvalue10) # 4 inversions
  
  head(dAF_df)
  
  # Create Granges for each of the SNP sets
  Salinity_GR <- GRanges(seqnames = BvA_join_pvalue10$CHROM, ranges = IRanges(start = BvA_join_pvalue10$POS, end = BvA_join_pvalue10$POS))
  Spawning_GR <- GRanges(seqnames = AvS_join_pvalue10$CHROM, ranges = IRanges(start = AvS_join_pvalue10$POS, end = AvS_join_pvalue10$POS))
  Inversions_GR <- GRanges(seqnames = sv_pvalue10$CHROM, ranges = IRanges(start = sv_pvalue10$POS, end = sv_pvalue10$POS))
  
  SNPs_GR <- GRanges(seqnames = dAF_df$CHROM, ranges = IRanges(start = dAF_df$POS, end = dAF_df$POS))
  in_Salinity_region <- findOverlaps(SNPs_GR, Salinity_GR)
  in_Spawning_region <- findOverlaps(SNPs_GR, Spawning_GR)
  in_invLocalPop <- findOverlaps(SNPs_GR, Inversions_GR)
  #in_Salinity_region_flank <- findOverlaps(SNPs_GR, sal_p20_regions_GR, maxgap = 1e4)  # add 1kb flank
  #in_Spawning_region_flank <- findOverlaps(SNPs_GR, spwn_p20_regions_GR, maxgap = 1e4)
  
  dAF_df[ ,"Salinity_region"] <- FALSE
  dAF_df[in_Salinity_region@from,"Salinity_region"] <- TRUE
  dAF_df[ ,"Spawning_region"] <- FALSE
  dAF_df[in_Spawning_region@from,"Spawning_region"] <- TRUE
  dAF_df[ ,"InvLocalPop_region"] <- FALSE
  dAF_df[in_invLocalPop@from,"InvLocalPop_region"] <- TRUE
  
}

head(dAF_df)
#inputDF <- dAF_df

## Loop through all the pairwise comparisons.
for(i in paired_comparisons){
  #i = 'dAF_BalticSpring_vs_NortheastAtlanticSpring'  # 'ChiSqr_IrelandBritain_vs_NortheastAtlantic' #'dAF_brackishNonMigratory_vs_brackishMigratory'
  #for(j in chr_names) {
  #i='dAF_WIreland_1a_vs_Others_3.4.8.5a.5b.9'
  #rollMeanWSwide <- 20
  
  ## Save filename.
  #library(tools)
  #filename <- file_path_sans_ext(basename(path_to_file1))
  #filename
  
  ## Generate single plots for each dataset and all chromosomes.
  p <- make.Manhattanplot(mode = 'dAF', inputDF = dAF_df, faiFile = fai_sorted_file,    # For dAF plot >>>>>>>>>
                          #mode = 'ChiSqr', inputDF = dAF_df, faiFile = fai_sorted_file,    # For ChiSqr plot >>>>>>>>>>
                          #mode = 'ChiSqr', inputFile = dAF_file, faiFile = fai_sorted_file,    # For ChiSqr plot >>>>>>>>>>
                          #mode = 'dAF', inputFile = dAF_file, faiFile = fai_sorted_file, 
                          comparison = i, colors = c('#737373','#bdbdbd'), #c('darkgrey','skyblue'), #colors = c('#cccccc','#969696'), 'darkturquoise'
                          chromosomeOrder = CHROM_order, addBonferroni = FALSE, # *****
                          addRollmean = TRUE, windowSizeRollMean = as.numeric(rollMeanWSwide), 
                          addSelectionRegions = FALSE  # *****
                          ) +
    labs(y = 'Delta allele frequency')
  #labs(y = '-log10Pvalue')    # For ChiSqr plot >>>>>>>>>>
  #labs(y = paste0('Delta allele frequency \n (',rollMeanWSwide,' bp rollmean)'))
  #theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
  #      axis.text = element_text(size = 5), axis.title.y = element_text(size = 6))
  
  ## Save image.
  #png('test_ranges.png', height = 60, width = 169, units='mm', res = 300)
  png(paste0('ManhattanPlot_',outputName,'-',i,'.ws',rollMeanWSwide,'.allChr.png'), 
      #png(paste0('/proj/snic2020-16-14/private/HorseMackerel/analysis/08-outlier-loci-detection/dAF/',dir_file,'/Manhattan_plot_dAF.',dir_file,'.',i,'.png'), 
      height = 60, width = 169, units='mm', res = 300
      #width = 950, height = 500   # width = 1900, height = 1000   #width = 10, height = 5, units='in', res = 150
  )
  #png(paste0('ManhattanPlot_',outputName,'-',i,'.allChr_wd_15.png'), height = 5, width = 16, units='cm', res = 150)
  print(p)
  dev.off()
  
  # }
  
  ## Manhattan plot with regions of interest highlighted.
  #ph <- p + 
  
}


```

Script `08-3-dAF-plotting.sh`:
```bash
#!/bin/bash
#SBATCH -A snic2020-5-94
#SBATCH -M snowy
#SBATCH -p core -n 4
#SBATCH -t 2-00:00:00
#SBATCH -J dAF_plot
#SBATCH -e dAF_plot_%J_%A_%a.err
#SBATCH -o dAF_plot_%J_%A_%a.out
#SBATCH --mail-type=all
#SBATCH --mail-user=angela.fuentespardo@gmail.com

# Load required software.
module load bioinfo-tools
module load R/3.6.1
module load R_packages/3.6.1

echo This is array job: $SLURM_ARRAY_TASK_ID
DAF_FILE='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/08-outlier-loci-detection/dAF/dAF_files.txt'
DAF_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $DAF_FILE)
OUTPUTDIR_FILE='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/08-outlier-loci-detection/dAF/outputdir_dAF.txt'
OUTPUTDIR_target=$(sed -n "$SLURM_ARRAY_TASK_ID"p $OUTPUTDIR_FILE)
echo -e This is the DAF_target: ${DAF_FILE}
echo -e This is the OUTPUTDIR_target: ${OUTPUTDIR_FILE}

# Set environment variables.
# Files.
WORK_DIR='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/analysis/08-outlier-loci-detection/dAF/plots'
FAI_FILE='/proj/snic2020-6-128/private/obtectus_poolseq/00-genome/A.obtectus_v2.0.fasta.fai'
OUTPUT_NAME='OTT_16_pools_SNPs_DPq5-99'
ROLLMEAN_WSWIDE='100'
SELECTION_REGIONS='NA'
SELECTION_REGIONS_TYPE='NA'
#SELECTION_REGIONS='/proj/snic2020-2-19/private/herring/users/angela/Pool-seq/analysis/08-outlier-loci-detection/dAF/SNPs_selection_regions_SalinitySpawning4INV.RData'
#SELECTION_REGIONS_TYPE='loci'
# Programs.
R_CODE='/proj/snic2020-6-128/private/obtectus_poolseq/users/angela/code/utility-code/08-2-create-ManhatPlots-dAF-allChr.R'

# Create directory for output files.
if [ -d $WORK_DIR/$OUTPUTDIR_target ]; then echo "OUTPUTDIR_target/ already exists"; else mkdir $WORK_DIR/$OUTPUTDIR_target; fi

cd $WORK_DIR/$OUTPUTDIR_target

# Run the R script that calculates dAF.
Rscript $R_CODE $WORK_DIR/$OUTPUTDIR_target $DAF_target $OUTPUT_NAME $FAI_FILE $ROLLMEAN_WSWIDE $SELECTION_REGIONS $SELECTION_REGIONS_TYPE

```
Submitted batch job 4351598 on cluster snowy
**Runtime: 00-00:40:54**??



# 2021-04-06

## Find informative genes with high dAF differences between Early and Late reproduction

R script `findGenRegInterest_OTT.R`:
```R
# Generate a list of genomic regions of interest in the OTT
# Code made by: Angela P. Fuentes-Pardo, e-mail: apfuentesp@gmail.com
# Date: 2021-03-19
#

# Clean environment space
rm(list = ls())

# Set environmental variables
working_dir <- '~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/08-2-outlier-detection_dAF'
gtf_file <- '~/Dropbox/PostDoc_UU/Projects/Beetles/data/00-genome/annotation/rc3.1_obtectus.gff'
dAF_file <- '~/Dropbox/PostDoc_UU/Projects/Beetles/analysis/08-2-outlier-detection_dAF/Summary_dAF.List_paired_comparisons_OTT_2020-09-01.OTT_16_pools.UG.SNPs.hf.GQ20.mono.miss20.mac3.DPq5-99.AD.Neff.AF.20210319.txt'

# Set working directory
setwd(working_dir)

# -----------------------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------------------

# Load file of dAF contrasts --------
library(data.table)

dAF_df <- fread(dAF_file, data.table = FALSE, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
head(dAF_df)
str(dAF_df)
dim(dAF_df)
#[1] 7,852,915      10

paired_comparisons <- grep("dAF_", names(dAF_df), value = TRUE)
paired_comparisons
#'dAF_EI.II.III.IV_vs_LV.VI.VII.VIII'

# Find which loci have a dAF >= 0.9, which are the most likely loci where each of the 4 bio replicates within the
# Early and Late treatments have the same allele fixed)
dAF_df_set <- dAF_df[dAF_df$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII >= 0.9, ]
dim(dAF_df_set)
#[1] 37,250    10
head(dAF_df_set)
length(unique(dAF_df_set$CHR))
#[1] 415
length(unique(dAF_df_set$gene_id))
#[1] 1335
length(unique(dAF_df_set$gene_name))
#[1] 691

# -----------------------------------------------------------------------------------------
# Exploratory analyses
# -----------------------------------------------------------------------------------------
# I will only focus on the SNPs within genes, thus I will subset the SNPs that have a gene name assigned
dAF_df_set <- dAF_df_set[!is.na(dAF_df_set$gene_name), ]
dim(dAF_df_set)
#[1] 16,546    10

# Create a dataframe where we will store the gene metadata information
gene_name_list <- unique(dAF_df_set$gene_name)
df <- data.frame(gene_id = rep(NA, length(gene_name_list)),
                 gene_name = gene_name_list,
                 gene_idname = rep(NA, length(gene_name_list)),
                 transcript_id = rep(NA, length(gene_name_list)),
                 protein_description = rep(NA, length(gene_name_list)),
                 GOterm = rep(NA, length(gene_name_list)),
                 SNPcount = rep(NA, length(gene_name_list))
                 )
# Assign the nomber of SNPs per gene, gene name and gene id to each of the genes of interest -----------
for(gen in gene_name_list){
  #gen="slc25a3"
  df[df$gene_name == gen, "SNPcount"] <- length(dAF_df_set[dAF_df_set$gene_name == gen, "gene_name"]) # count the number of SNPs per gene
  df[df$gene_name == gen, "gene_id"] <- unique(dAF_df_set[dAF_df_set$gene_name == gen, "gene_id"]) # find unique gene_id
  df[df$gene_name == gen, "gene_idname"] <- unique(dAF_df_set[dAF_df_set$gene_name == gen, "gene_idname"]) # find unique gene_name
}
head(df)

# Make a histogram of the number of SNPs per gene with annotated gene symbol
library(ggplot2)
pdf("Histogram_SNP_count_per_gene_with_symbol_dAFgt0.9_OTT_16pools.pdf")
ggplot(df, aes(x=SNPcount)) +
  geom_histogram(color="black", fill="lightblue") +
  xlab("Number of SNPs per gene") +
  theme_bw() +
  theme(text = element_text(size=14))
#hist(df$count)
dev.off()

# Sort genes by descending number of SNPs
df <- df[order(df$SNPcount, decreasing = TRUE), ]
head(df, n = 20)

# Assign transcript id, protein description, and GOterms to each gene of interest -----------

### Option 1 ###
# Load the GTF file in GRanges format
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)

GTFgr <- rtracklayer::import(gtf_file)
GTFgr

# Select the mRNA annotations only, becasue they have the protein description and the Parent gene associated to each transcript
GTFgr_set <- GTFgr[GTFgr$type == "mRNA", ]
GTFmtDNA_df <- as.data.frame(mcols(GTFgr_set), stringsAsFactors = FALSE) # extract metadata colums as a dataframe
str(GTFmtDNA_df)
dim(GTFmtDNA_df)
#[1] 68,812    23
head(GTFmtDNA_df)

# Convert the character list columns Ontology_term and Parent to acharacter vector
GTFmtDNA_df$Parent <- unlist(GTFmtDNA_df$Parent, use.names=FALSE)

# Subset df to the target genes
geneid_list <- df$gene_id
GTFmtDNA_df <- GTFmtDNA_df[GTFmtDNA_df$Parent %in% geneid_list, ]
dim(GTFmtDNA_df)
#[1] 1989   23
head(GTFmtDNA_df)

# Assign the transcript code, gene description, and GO terms to the genes of interest
for(i in geneid_list){
  #i="ACAOBTG00000003090"
  #k="ACAOBTM00000020198"
  df[df$gene_id == i, "transcript_id"] <- paste(GTFmtDNA_df[GTFmtDNA_df$Parent == i, "ID"], collapse = ",")
  df[df$gene_id == i, "protein_description"] <- paste(unique(GTFmtDNA_df[GTFmtDNA_df$Parent == i, "product"]), collapse = ",")
  tmp <- unlist(lapply(GTFmtDNA_df[GTFmtDNA_df$Parent == i, "Ontology_term"], function(x) if(identical(x,character(0))) NA else x))
  df[df$gene_id == i, "GOterm"] <- paste(tmp[!is.na(tmp)], collapse = ",")

  #GO_df[GO_df$transcript_id == k, "Parent"] <- GTFmtDNA_df[GTFmtDNA_df$ID == k, "Parent"]
  #GO_df[GO_df$transcript_id == k, "product"] <- GTFmtDNA_df[GTFmtDNA_df$ID == k, "product"]
}

head(df)

# Save the metadata file
write.table(df, "Metadata_of_genes_with_symbol_dAFgt0.9_OTT_16pools.txt", sep="\t", quote = F, col.names = T, row.names = F)

### Option 2 ###
# Load the metadata file listing the GO terms associated to each trancript
#GO_df <- read.table('/Users/angfu103/Dropbox/PostDoc_UU/Projects/Beetles/data/00-genome/Ontology_term.txt', header=F, sep="\t")
#colnames(GO_df) <- c("transcript_id","GOterm")
#head(GO_df)
#GO_df$gene_id <- gsub("M","G",GO_df$gene_id)


# Create ranges that correspond to all the SNPs in the same gene -----------
tmp <- as.data.table(dAF_df_set)
ranges_by_geneID <- as.data.frame(tmp[, as.data.table(reduce(IRanges(min(POS), max(POS)))), by = .(gene_id, gene_name, CHROM)], stringsAsFactors = FALSE)
head(ranges_by_geneID)

# Save the genomic ranges as a RDS file
saveRDS(ranges_by_geneID, file="genomic_ranges_by_geneIDName_dAFgt0.9_E_vs_L_OTT_16pools.RDS")


############################## END #########################################

# Make ranges out of the list of SNPs of interest
library(GenomicRanges)

# Convert SNPs into GRanges
SNPs_sig_GR <- GRanges(seqnames = dAF_df_set$CHROM,
                       ranges = IRanges(start = dAF_df_set$POS, end = dAF_df_set$POS),
                       dAF = dAF_df_set$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII,
                       gene_id = dAF_df_set$gene_id,
                       gene_name = dAF_df_set$gene_name
)
length(unique(seqnames(SNPs_sig_GR)))
# 37250 ranges
#[1] 281 ranges -- for only SNPs with an annotated gene_name

# Merge close SNPs to form ranges with a minimum gap of 100 kb
# Note that reduce() returns an object of the same type as x containing reduced ranges for each distinct (seqname, strand) pairing
# min.gapwidth: Ranges separated by a gap of at least min.gapwidth positions are not merged
sig_GenomRanges_GR <- GenomicRanges::reduce(SNPs_sig_GR, min.gapwidth = 100000)
#395 ranges -- for only SNPs with an annotated gene_name
#587 ranges
#sig_GenomRanges_GR <- GenomicRanges::reduce(SNPs_sig_GR, min.gapwidth = 10000)
#1488 ranges
#sig_GenomRanges_GR <- GenomicRanges::reduce(SNPs_sig_GR, min.gapwidth = 1000)
#5501 ranges
#sig_GenomRanges_GR <- GenomicRanges::reduce(SNPs_sig_GR, min.gapwidth = 100)
#15382 ranges
#sig_GenomRanges_GR <- GenomicRanges::reduce(SNPs_sig_GR, ignore.strand=TRUE)
#36128 ranges

# Require ranges span a minimum of 100 bp
#sig_GenomRanges_GR <- sig_GenomRanges_GR[ranges(sig_GenomRanges_GR)@width >= 100]  
#sig_GenomRanges_GR
#97 ranges

## Count the number of SNPs in each range
sig_GenomRanges_GR$totalSNPcount <- countOverlaps(sig_GenomRanges_GR, SNPs_sig_GR)
sig_GenomRanges_df <- as.data.frame(sig_GenomRanges_GR, stringsAsFactors = FALSE)
sig_GenomRanges_df$seqnames <- as.character(sig_GenomRanges_df$seqnames)
#str(sig_GenomRanges_df)  # make sure the chromosomes are characters
head(sig_GenomRanges_df)
dim(sig_GenomRanges_df)
#[1] 587   6

```
