#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

#Required packages & functions
require(vcfR)
source("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/A_ob_LO_funcs.R")


 #### affected by mislablelling, correction below.
#Line origin calling of F2 beetles; version for UPPMAX###
#LO_pedigree <- read.table("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/A_ob_pedigree_v2.txt", sep = "\t", stringsAsFactors = F, header = T)
#founder_LO_df <- data.frame(ID = LO_pedigree$ID[LO_pedigree$generation == 0], line = NA)
#founder_LO_df$line[grep("[0-9]-M", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\1",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])
#founder_LO_df$line[grep("[0-9]-F", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\2",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])
####

## Adjusted version
LO_pedigree <- read.table("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/A_ob_pedigree_corrected_v1.txt", sep = "\t", stringsAsFactors = F, header = T)
founder_LO_df <- data.frame(ID = LO_pedigree$ID[LO_pedigree$generation == 0], line = NA)
founder_LO_df$line[grep("[0-9]-M", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\1",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])
founder_LO_df$line[grep("[0-9]-F", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\2",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])
##Correcting for mixed-up samples
founder_LO_df$line[founder_LO_df$ID == "Sample_VA-3193-LE1-M"] <- "E"
founder_LO_df$line[founder_LO_df$ID == "Sample_VA-3193-LE1-F"] <- "L"
##


#offspring_test_vcf <- read.vcfR(file = "./beetle_offspring_extended_test.vcf")

#founder_test_vcf <- read.vcfR(file = "./beetle_founders_extended_test.vcf")
max_variant_count <- 10933034
#max_variant_count <- 3e5
current_skip <- 0
chunk_size <- 1e5
n_chunks <- ceiling(max_variant_count/chunk_size)

for(chunk_i in 1:n_chunks){
  
  current_offspring_vcf <- read.vcfR(file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/beetle_offspring_full.vcf", nrows = chunk_size, skip = current_skip)
  offspring_test_GTs <- extract.haps(current_offspring_vcf, unphased_as_NA = F) #extract.gt(offspring_test_vcf)
  rm(current_offspring_vcf)
  
  current_founder_vcf <- read.vcfR(file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/founders/vcf/OTT_8_F0.HC.SNPs.hf.DP3.GQ20.mono.miss20.vcf.gz", nrows = chunk_size, skip = current_skip)
  founder_test_GTs <- extract.haps(current_founder_vcf, unphased_as_NA = F) #extract.gt(founder_test_vcf)
  rm(current_founder_vcf)
  
  current_skip <- current_skip + chunk_size
  
  #####
  win_seqs <- unique(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(founder_test_GTs)))
  win_size <- 100000
  for(win_seq in win_seqs){
    tmp_win_GTs <- founder_test_GTs[sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(founder_test_GTs)) == win_seq, ]
    tmp_win_starts <- seq(from = min(as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(tmp_win_GTs)))), to = max(as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(tmp_win_GTs)))), by = win_size)
    
    tmp_win_ends <- tmp_win_starts + win_size - 1
    tmp_win_ends[length(tmp_win_ends)] <- max(as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(tmp_win_GTs))))
    tmp_win_df <- data.frame(chr = win_seq, start = tmp_win_starts, end = tmp_win_ends)
    tmp_win_df$summary_string <- paste0(tmp_win_df$chr, ":", tmp_win_df$start, "-", tmp_win_df$end)
    if(win_seq == win_seqs[1]){
      LO_win_df <- tmp_win_df
    } else {
      LO_win_df <- rbind(LO_win_df, tmp_win_df)
    }
  }
  
  Lo_est_list <- list()
  for(current_win in LO_win_df$summary_string){
    #Subsetting the genotype data to current window
    win_site_filter_F0 <-  sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(founder_test_GTs)) == sub("([0-9A-Z]+)[:]([0-9e+]+)[-]([0-9e+]+)", "\\1", current_win) &
      as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(founder_test_GTs))) >= as.integer(sub("([0-9A-Z]+)[:]([0-9e+]+)[-]([0-9e+]+)", "\\2", current_win)) &
      as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(founder_test_GTs))) < as.integer(sub("([0-9A-Z]+)[:]([0-9e+]+)[-]([0-9e+]+)", "\\3", current_win))
    
    win_founder_GTs <- founder_test_GTs[win_site_filter_F0,]
    
    win_site_filter_F2 <-  sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(offspring_test_GTs)) == sub("([0-9A-Z]+)[:]([0-9e+]+)[-]([0-9e+]+)", "\\1", current_win) &
      as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(offspring_test_GTs))) >= as.integer(sub("([0-9A-Z]+)[:]([0-9e+]+)[-]([0-9e+]+)", "\\2", current_win)) &
      as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(offspring_test_GTs))) < as.integer(sub("([0-9A-Z]+)[:]([0-9e+]+)[-]([0-9e+]+)", "\\3", current_win))
    
    win_offspring_GTs <- offspring_test_GTs[win_site_filter_F2,]
    
    #Estimating LOs in current window
    LO_est_df <- data.frame(ID = LO_pedigree$ID[LO_pedigree$generation == 2], window = current_win , LL = NA, LE = NA, EE = NA)
    if("matrix" %in% class(win_offspring_GTs) & "matrix" %in% class(win_founder_GTs)){
      for(F2_ind in LO_est_df$ID){
        LO_est_df[LO_est_df$ID == F2_ind,c("LL", "LE", "EE")] <- est_LO_segment(target_offspring = F2_ind,
                                                                                founder_GT_df = win_founder_GTs, 
                                                                                offspring_GT_df = win_offspring_GTs,
                                                                                pedigree_df = LO_pedigree,
                                                                                founder_line_df = founder_LO_df)
      }
    } else{
      LO_est_df$LL <- 0
      LO_est_df$LE <- 0
      LO_est_df$EE <- 0
    }
      
    #Processing & summarizing the results
    LO_est_df$est_pos <- rowSums(LO_est_df[,c("LL", "LE", "EE")])
    LO_est_df$call <- c("LL", "LE", "EE")[max.col(LO_est_df[,c("LL", "LE", "EE")])]
    LO_est_df$call[LO_est_df$est_pos < 5] <- "None" 
    LO_est_df$proportion <- pmax(LO_est_df[,"LL"], LO_est_df[,"LE"],LO_est_df[,"EE"])/LO_est_df$est_pos
    LO_est_df$proportion[LO_est_df$est_pos < 5] <- 0
    LO_est_df$call[LO_est_df$proportion >= 0.5 & LO_est_df$proportion < 0.75] <- "LE" #Accounting for the fact that low coverage favours homzygous calls
    
    #Collecting the results
    Lo_est_list[[current_win]] <- LO_est_df
  }
  #LO_outfile <- paste0("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_results/A_ob_Lo_est_list_chunk_", chunk_i,".RData")
  LO_outfile <- paste0("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_results_corrected/A_ob_Lo_est_list_chunk_", chunk_i,".RData")
  save(Lo_est_list, file = LO_outfile)
  
}
####




