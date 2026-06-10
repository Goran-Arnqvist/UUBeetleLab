#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

#Line origin calling of F2 beetles
filtered_VCF_files <- read.table("~/Projects/Beetle/data/offspring_LO_calling/offspring_filtered_vcf_list.txt", sep = "", stringsAsFactors = F, header = F)
filtered_VCF_files$sample_ID <- sub("(.+)[.]filtered[.]vcf[.]gz", "\\1", filtered_VCF_files$V9)
filtered_VCF_files$sample_ID_file <- paste0("~/Projects/Beetle/data/offspring_LO_calling/sample_ID_files/", filtered_VCF_files$sample_ID, ".txt")
filtered_VCF_files$vcf_file <- paste0("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/offspring/filtered_vcf/",filtered_VCF_files$V9)

for (i in 1:length(filtered_VCF_files$sample_ID)){
  write(x = filtered_VCF_files$sample_ID[i], file =  filtered_VCF_files$sample_ID_file[i])
}
write(x = paste0("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/offspring/sample_ID_files/", filtered_VCF_files$sample_ID, ".txt"), file = "~/Projects/Beetle/data/offspring_LO_calling/sample_ID_file_list.txt")
write(x = filtered_VCF_files$vcf_file, file = "~/Projects/Beetle/data/offspring_LO_calling/vcf_file_list.txt")
write(x = filtered_VCF_files$sample_ID, file = "~/Projects/Beetle/data/offspring_LO_calling/sample_ID_list.txt")




#### NOTE! A sample mixup has been detected, hence the pedigree below is incorrect. New version implemented below.

#Pedigree

raw_phenotypes <- read.csv("~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/phenotypes/phenotypes_beetles.csv", stringsAsFactors = F)
raw_phenotypes$parentage_1 <- sub("([LE0-9]+)_([LE0-9]+)","\\1", raw_phenotypes$Family_cross)
raw_phenotypes$parentage_2 <- sub("([LE0-9]+)_([LE0-9]+)","\\2", raw_phenotypes$Family_cross)
raw_phenotypes$reformed_ID <- sub("-","_", raw_phenotypes$ID)

filtered_VCF_files$short_ID <- sub("([A-Z0-9]+)[.]([A-Z0-9_]+)","\\1",filtered_VCF_files$sample_ID)
raw_phenotypes$VCF_ID <- filtered_VCF_files$sample_ID[match(raw_phenotypes$reformed_ID, filtered_VCF_files$short_ID)] 

#Parental sample names:
parent_samples <- c("Sample_VA-3193-EL1-F", "Sample_VA-3193-EL1-M", 
                    "Sample_VA-3193-EL4-F", "Sample_VA-3193-EL4-M", 
                    "Sample_VA-3193-LE1-F", "Sample_VA-3193-LE1-M", 
                    "Sample_VA-3193-LE2-F", "Sample_VA-3193-LE2-M")


raw_phenotypes$p1_sire <- sapply(paste0("-", raw_phenotypes$parentage_1,"-M"), FUN = grep, x = parent_samples, value = T)
raw_phenotypes$p1_dam <- sapply(paste0("-", raw_phenotypes$parentage_1,"-F"), FUN = grep, x = parent_samples, value = T)
raw_phenotypes$p2_sire <- sapply(paste0("-", raw_phenotypes$parentage_2,"-M"), FUN = grep, x = parent_samples, value = T)
raw_phenotypes$p2_dam <- sapply(paste0("-", raw_phenotypes$parentage_2,"-F"), FUN = grep, x = parent_samples, value = T)

#Compiling the final pedigree
F2_pedigree <- raw_phenotypes[!is.na(raw_phenotypes$VCF_ID), ]
 
founder_pedigree <- data.frame(ID = unique(c(raw_phenotypes$p1_sire, raw_phenotypes$p1_dam)), Sire = 0, Dam = 0, Sex = NA, Weight = NA, Lifespan = NA, stringsAsFactors = F)
founder_pedigree$Sex[grep("-M$", founder_pedigree$ID)] <- 1
founder_pedigree$Sex[grep("-F$", founder_pedigree$ID)] <- 2

F1_IDs <- unique(c(paste(F2_pedigree$Family, F2_pedigree$parentage_1, sep = "_"),
                 paste(F2_pedigree$Family, F2_pedigree$parentage_2, sep = "_")))

F1_IDs <- paste(rep(F1_IDs, each = 2), rep(c("sire", "dam"), times = length(F1_IDs)), sep = "_")

F1_pedigree <- data.frame(ID = F1_IDs, Sire = NA, Dam = NA, Sex = NA, Weight = NA, Lifespan = NA, stringsAsFactors = F)
F1_pedigree$Sex[grep("sire", F1_pedigree$ID)] <- 1
F1_pedigree$Sex[grep("dam", F1_pedigree$ID)] <- 2

F1_pedigree$Sire <- sapply(sub("([A-Z0-9]+)_([LE0-9]+)_(sire|dam)","\\2",F1_pedigree$ID), FUN = grep, x = founder_pedigree$ID[founder_pedigree$Sex == 1], value = T)
F1_pedigree$Dam <- sapply(sub("([A-Z0-9]+)_([LE0-9]+)_(sire|dam)","\\2",F1_pedigree$ID), FUN = grep, x = founder_pedigree$ID[founder_pedigree$Sex == 2], value = T)

F2_pedigree$Sire <- sapply(paste(F2_pedigree$Family, F2_pedigree$parentage_1, sep = "_"), FUN = grep, x = F1_pedigree$ID[F1_pedigree$Sex == 1], value = T)
F2_pedigree$Dam <- sapply(paste(F2_pedigree$Family, F2_pedigree$parentage_2, sep = "_"), FUN = grep, x = F1_pedigree$ID[F1_pedigree$Sex == 2], value = T)

F2_pedigree_out <- F2_pedigree[,c("VCF_ID", "Sire", "Dam", "Sex", "Weight", "Lifespan")]
names(F2_pedigree_out)[1] <- "ID"
F2_pedigree_out$Sex[F2_pedigree_out$Sex == "M"] <- 1
F2_pedigree_out$Sex[F2_pedigree_out$Sex == "F"] <- 2
F2_pedigree_out$Sex <- as.integer(F2_pedigree_out$Sex)

output_pedigree <- rbind(founder_pedigree, F1_pedigree, F2_pedigree_out)

write(x = output_pedigree$ID, file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_samples.txt", ncolumns = 1)
write.table(x = output_pedigree, file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_pedigree.txt", sep ="\t", quote = F, row.names = F)

#Numeric version
output_pedigree$numeric_ID <- 1:dim(output_pedigree)[1]
output_pedigree$numeric_Sire <- output_pedigree$numeric_ID[match(output_pedigree$Sire, output_pedigree$ID)]
output_pedigree$numeric_Dam <- output_pedigree$numeric_ID[match(output_pedigree$Dam, output_pedigree$ID)]
output_pedigree$numeric_Sire[is.na(output_pedigree$numeric_Sire)] <- 0
output_pedigree$numeric_Dam[is.na(output_pedigree$numeric_Dam)] <- 0
output_pedigree$generation <- 2
output_pedigree$generation[grep("Sample_VA", output_pedigree$Sire)] <- 1
output_pedigree$generation[output_pedigree$Sire == 0] <- 0
write.table(x = output_pedigree[,c("numeric_ID", "numeric_Sire", "numeric_Dam", "Sex", "Weight", "Lifespan")], file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_pedigree_int.txt", sep ="\t", quote = F, row.names = F, col.names = c("ID", "Sire", "Dam", "Sex", "Weight", "Lifespan" ))

output_pedigree$mat_Grand_Dam <- output_pedigree$Dam[match(output_pedigree$Dam, output_pedigree$ID)]
output_pedigree$mat_Grand_Sire <- output_pedigree$Sire[match(output_pedigree$Dam, output_pedigree$ID)]
output_pedigree$pat_Grand_Dam <- output_pedigree$Dam[match(output_pedigree$Sire, output_pedigree$ID)]
output_pedigree$pat_Grand_Sire <- output_pedigree$Sire[match(output_pedigree$Sire, output_pedigree$ID)]
write.table(x = output_pedigree, file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_pedigree_v2.txt", sep ="\t", quote = F, row.names = F, col.names = T)



###
offspring_test_vcf <- read.vcfR(file = "~/Projects/Beetle/data/offspring_LO_calling/test_vcfs/beetle_offspring_test.vcf")
offspring_test_GTs <- extract.haps(offspring_test_vcf, unphased_as_NA = F) #extract.gt(offspring_test_vcf)
founder_test_vcf <- read.vcfR(file = "~/Projects/Beetle/data/offspring_LO_calling/test_vcfs/beetle_founders_test.vcf")
founder_test_GTs <- extract.haps(founder_test_vcf, unphased_as_NA = F) #extract.gt(founder_test_vcf)

founder_LO_df <- data.frame(ID = output_pedigree$ID[output_pedigree$generation == 0], line = NA)
founder_LO_df$line[grep("[0-9]-M", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\1",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])
founder_LO_df$line[grep("[0-9]-F", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\2",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])


###test call###
est_LO_segment(target_offspring = "DC17_31.D76_N702",   #AB14_21.D34_N705, CD11_9.D96_N701, DC17_31.D76_N702
                founder_GT_df = founder_test_GTs, 
                offspring_GT_df = offspring_test_GTs,
                pedigree_df = output_pedigree,
                founder_line_df = founder_LO_df)
####
##Main calling loop

win_seqs <- unique(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(founder_test_GTs)))
win_size <- 10000
for(win_seq in win_seqs){
  tmp_win_GTs <- founder_test_GTs[sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(founder_test_GTs)) == win_seq, ]
  tmp_win_starts <- seq(from = 1, to = max(as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(tmp_win_GTs)))), by = win_size)

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
#current_win <- "1:5000-15000"
  #Subsetting the genotype data to current window
  win_site_filter_F0 <-  sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(founder_test_GTs)) == sub("([0-9A-Z]+)[:]([0-9]+)[-]([0-9]+)", "\\1", current_win) &
                        as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(founder_test_GTs))) >= as.integer(sub("([0-9A-Z]+)[:]([0-9]+)[-]([0-9]+)", "\\2", current_win)) &
                        as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(founder_test_GTs))) < as.integer(sub("([0-9A-Z]+)[:]([0-9]+)[-]([0-9]+)", "\\3", current_win))

  win_founder_GTs <- founder_test_GTs[win_site_filter_F0,]

  win_site_filter_F2 <-  sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\1", rownames(offspring_test_GTs)) == sub("([0-9A-Z]+)[:]([0-9]+)[-]([0-9]+)", "\\1", current_win) &
                        as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(offspring_test_GTs))) >= as.integer(sub("([0-9A-Z]+)[:]([0-9]+)[-]([0-9]+)", "\\2", current_win)) &
                        as.integer(sub("[A-Za-z]+_([0-9]+)_([0-9]+)", "\\2", rownames(offspring_test_GTs))) < as.integer(sub("([0-9A-Z]+)[:]([0-9]+)[-]([0-9]+)", "\\3", current_win))

  win_offspring_GTs <- offspring_test_GTs[win_site_filter_F2,]

  #Estimating LOs in current window
  LO_est_df <- data.frame(ID = output_pedigree$ID[output_pedigree$generation == 2], window = current_win , LL = NA, LE = NA, EE = NA)
  for(F2_ind in LO_est_df$ID){
    LO_est_df[LO_est_df$ID == F2_ind,c("LL", "LE", "EE")] <- est_LO_segment(target_offspring = F2_ind,
                                                                          founder_GT_df = win_founder_GTs, 
                                                                          offspring_GT_df = win_offspring_GTs,
                                                                          pedigree_df = output_pedigree,
                                                                          founder_line_df = founder_LO_df)
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
####

est_LO_segment <- function(target_offspring, founder_GT_df, offspring_GT_df, pedigree_df, founder_line_df){
  
  
  #Proccess the pedigree to extract genotypes of the appropriate founders 
  current_founders <- as.character(pedigree_df[pedigree_df$ID == target_offspring, c("pat_Grand_Sire","pat_Grand_Dam","mat_Grand_Sire","mat_Grand_Dam")])
  
  L_founders <- founder_line_df$ID[(founder_line_df$ID %in% current_founders) & founder_line_df$line == "L"]
  L_founder_cols <- match(paste(rep(L_founders, each = 2), c(0,1), sep = "_"), colnames(founder_GT_df))
  E_founders <- founder_line_df$ID[(founder_line_df$ID %in% current_founders) & founder_line_df$line == "E"]
  E_founder_cols <- match(paste(rep(E_founders, each = 2), c(0,1), sep = "_"), colnames(founder_GT_df))
  
  #Finding positions fixed within line
  L_fixed_pos <- which(rowSums(founder_GT_df[,L_founder_cols] == founder_GT_df[,L_founder_cols[1]]) == length(L_founder_cols))
  E_fixed_pos <- which(rowSums(founder_GT_df[,E_founder_cols] == founder_GT_df[,E_founder_cols[1]]) == length(E_founder_cols))
  potentital_pos <- 1:dim(founder_GT_df)[1] %in% intersect(L_fixed_pos, E_fixed_pos)
  
  eval_pos <- rownames(founder_GT_df)[which(founder_GT_df[,E_founder_cols[1]] != founder_GT_df[,L_founder_cols[1]] & potentital_pos)]
    
  
  #Evaluating the target offspring
  eval_offspring_GTs <- offspring_GT_df[eval_pos,grep(target_offspring, colnames(offspring_GT_df))]
  eval_L_GTs <- founder_GT_df[eval_pos,L_founder_cols[1]]
  eval_E_GTs <- founder_GT_df[eval_pos,E_founder_cols[1]]
  
  LL_pos <- which(eval_offspring_GTs[,1] == eval_L_GTs & eval_offspring_GTs[,2] == eval_L_GTs)
  EE_pos <- which(eval_offspring_GTs[,1] == eval_E_GTs & eval_offspring_GTs[,2] == eval_E_GTs)
  LE_pos <- which((eval_offspring_GTs[,1] == eval_L_GTs & eval_offspring_GTs[,2] == eval_E_GTs) & (eval_offspring_GTs[,2] == eval_L_GTs & eval_offspring_GTs[,1] == eval_E_GTs))
  out_array <- c(length(LL_pos), length(LE_pos), length(EE_pos))
  names(out_array) <- c("LL", "LE", "EE")
  return(out_array)
}

###Data collection code.
LO_comp_list <- list()
#LO_data_files <- dir("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_results/", full.names = T)
LO_data_files <- dir("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_results_corrected/", full.names = T)
LO_data_files <- LO_data_files[grepl("chunk_[0-9]+[.]RData",LO_data_files)]

for(LO_data_file in LO_data_files){
  load(LO_data_file)
  LO_comp_list <- append(LO_comp_list,Lo_est_list)
  rm(Lo_est_list)
}

#save(LO_comp_list, file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_list.RData")
save(LO_comp_list, file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_list_corrected.RData")


comp_LO_df <- do.call(rbind, LO_comp_list)
#reshape_test_LO_df <- reshape(test_LO_df[, c("ID", "window", "call")], direction="wide", timevar="window", idvar = "ID")
reshape_LO_df <- reshape(comp_LO_df[, c("ID", "window", "call")], direction="wide", idvar="window", timevar = "ID")
reshape_LO_df$chr <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\1", reshape_LO_df$window))
reshape_LO_df$win_start <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\2", reshape_LO_df$window))
reshape_LO_df <- reshape_LO_df[order(reshape_LO_df$chr, reshape_LO_df$win_start), ]

#save(reshape_LO_df, file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_reshape_df.RData")
save(reshape_LO_df, file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_reshape_df_corrected.RData")

#load("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_reshape_df.RData")
load("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_reshape_df_corrected.RData")

reshape_LO_df_num <- reshape_LO_df
reshape_LO_df_num[reshape_LO_df == "None"] <- NA
reshape_LO_df_num[reshape_LO_df == "EE"] <- 0
reshape_LO_df_num[reshape_LO_df == "LE"] <- 1
reshape_LO_df_num[reshape_LO_df == "LL"] <- 2

call_cols <- grep("call", colnames(reshape_LO_df_num))
for(i in call_cols){
  reshape_LO_df_num[,i] <- as.integer(reshape_LO_df_num[,i])
}

#save(reshape_LO_df_num, file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_reshape_012_df.RData")
save(reshape_LO_df_num, file = "/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/LO_compiled_reshape_012_df_corrected.RData")


#load("~/Projects/Beetle/data/offspring_LO_calling/LO_compiled_reshape_012_df.RData") <--- Obsolete version with mixed-up parentals
#load("~/Projects/Beetle/data/offspring_LO_calling/LO_compiled_reshape_012_df_corrected.RData")

##Matching LO and Phenotypes
high_cov_windows <- which(rowSums(!is.na(reshape_LO_df_num[,grep("call[.]", colnames(reshape_LO_df_num))])) > 400)

QTL_results_df <- data.frame(win_ID = names(high_cov_windows), p_val = NA, est = NA, p_val_w = NA, est_w = NA)
QTL_results_df$chr <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\1", QTL_results_df$win_ID))
QTL_results_df$win_start <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\2", QTL_results_df$win_ID))


for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismathc detected!")

  current_model <- glm(Lifespan~LO_call + Sex + 1, data = tmp_QTL_df, na.action = na.exclude)
  currentLO_p_val <- summary(current_model)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_estimate <- summary(current_model)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "est"] <- currentLO_estimate
  QTL_results_df[i, "p_val"] <- currentLO_p_val
  
  current_model_w <- glm(Lifespan~LO_call + Sex + Weight + 1, data = tmp_QTL_df, na.action = na.exclude)
  currentLO_p_val_w  <- summary(current_model_w)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_estimate_w  <- summary(current_model_w)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "est_w"] <- currentLO_estimate_w 
  QTL_results_df[i, "p_val_w"] <- currentLO_p_val_w 
}


QTL_results_df$weight_p <- NA
QTL_results_df$weight_est <- NA

for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismathc detected!")
  
  current_model_weight <- glm(Weight~LO_call + Sex + 1, data = tmp_QTL_df, na.action = na.exclude)
  currentLO_p_weight <- summary(current_model_weight)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_weight <- summary(current_model_weight)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "weight_est"] <- currentLO_e_weight
  QTL_results_df[i, "weight_p"] <- currentLO_p_weight
}

QTL_results_df$male_p <- NA
QTL_results_df$male_est <- NA

QTL_results_df$female_p <- NA
QTL_results_df$female_est <- NA

for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismathc detected!")
  
  current_model_male <- glm(Lifespan~LO_call + Weight + 1, data = tmp_QTL_df[tmp_QTL_df$Sex == 1,], na.action = na.exclude)
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "male_est"] <- currentLO_e_male
  QTL_results_df[i, "male_p"] <- currentLO_p_male
  
  current_model_female <- glm(Lifespan~LO_call + Weight + 1, data = tmp_QTL_df[tmp_QTL_df$Sex == 2,], na.action = na.exclude)
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "female_est"] <- currentLO_e_female
  QTL_results_df[i, "female_p"] <- currentLO_p_female
}

QTL_results_df$weight_female_p <- NA
QTL_results_df$weight_female_est <- NA
QTL_results_df$weight_male_p <- NA
QTL_results_df$weight_male_est <- NA

for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismathc detected!")
  
  current_model_male_w <- glm(Weight~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex == 1,], na.action = na.exclude)
  currentLO_male_p_w <- summary(current_model_male_w)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_male_est_w <- summary(current_model_male_w)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "weight_male_est"] <- currentLO_male_est_w
  QTL_results_df[i, "weight_male_p"] <- currentLO_male_p_w
  
  current_model_female_w <- glm(Weight~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex == 2,], na.action = na.exclude)
  currentLO_female_p_w <- summary(current_model_female_w)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_female_est_w <- summary(current_model_female_w)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "weight_female_est"] <- currentLO_female_est_w
  QTL_results_df[i, "weight_female_p"] <- currentLO_female_p_w
}

QTL_results_df$male_no_w_p <- NA
QTL_results_df$male_no_w_est <- NA

QTL_results_df$female_no_w_p <- NA
QTL_results_df$female_no_w_est <- NA

for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismathc detected!")
  
  current_model_male <- glm(Lifespan~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex == 1,], na.action = na.exclude)
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "male_no_w_est"] <- currentLO_e_male
  QTL_results_df[i, "male_no_w_p"] <- currentLO_p_male
  
  current_model_female <- glm(Lifespan~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex == 2,], na.action = na.exclude)
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "female_no_w_est"] <- currentLO_e_female
  QTL_results_df[i, "female_no_w_p"] <- currentLO_p_female
}

#save(QTL_results_df, file = "~/Projects/Beetle/data/offspring_LO_calling/QTL_results_v1.RData") <-- affected by sample mix-up
save(QTL_results_df, file = "~/Projects/Beetle/data/offspring_LO_calling/QTL_results_corrected_v2.RData")

fam_list <- unique(sub("call[.]([A-Z]+)[0-9]+.+","\\1",names(reshape_LO_df_num[,call_cols])))
call_cols_full <- grep("call", colnames(reshape_LO_df_num), value = T)
AB_cols <- grep("[.]AB[0-9]+.+", call_cols_full, value = T)
BA_cols <- grep("[.]BA[0-9]+.+", call_cols_full, value = T)
CD_cols <- grep("[.]CD[0-9]+.+", call_cols_full, value = T)
DC_cols <- grep("[.]DC[0-9]+.+", call_cols_full, value = T)

QTL_df_out <- F2_pedigree_out
for(i in 1:length(high_cov_windows)){ #
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  #tmp_QTL_df <- F2_pedigree_out
  current_win <- names(high_cov_windows)[i]
  QTL_df_out[,current_win] <- NA
  QTL_df_out[,current_win] <-  as.integer(tmp_win_L0s)
}
QTL_df_out_header <- QTL_df_out[1:2,]
QTL_df_out_header[1,] <- sub("([0-9]+)[:]([0-9]+)[-].+", "\\1",names(QTL_df_out_header))
QTL_df_out_header[1,1:6] <- ""
QTL_df_out_header[2,] <- sub("([0-9]+)[:]([0-9]+)[-].+", "\\2",names(QTL_df_out_header))
QTL_df_out_header[2,1:6] <- ""

QTL_df_out <- rbind(QTL_df_out_header, QTL_df_out)
write.table(x = QTL_df_out, file = "~/Projects/Beetle/data/offspring_LO_calling/F2_LO_calls.txt", row.names = F, col.names = T, quote = F, sep = "\t")

#pdf("~/Projects/Beetle/doc/F2_QTL_profile.pdf", width = 10) <-- affected by sample mix-up
#pdf("~/Projects/Beetle/doc/F2_QTL_profile_v2.pdf", width = 10) <-- affected by sample mix-up
pdf("~/Projects/Beetle/doc/F2_QTL_profile_corrected_v1.pdf", width = 10)

hist(colSums(!is.na(reshape_LO_df_num[,call_cols]), na.rm = T), breaks = 75, ylab = "Number of F2 individuals", xlab = "Number of successfully called 100kb windows", main = "All F2s")
hist(colSums(!is.na(reshape_LO_df_num[,AB_cols]), na.rm = T), breaks = 75, ylab = "Number of F2 individuals", xlab = "Number of successfully called 100kb windows", main = "Family AB")
hist(colSums(!is.na(reshape_LO_df_num[,BA_cols]), na.rm = T), breaks = 75, ylab = "Number of F2 individuals", xlab = "Number of successfully called 100kb windows", main = "Family BA")
hist(colSums(!is.na(reshape_LO_df_num[,CD_cols]), na.rm = T), breaks = 75, ylab = "Number of F2 individuals", xlab = "Number of successfully called 100kb windows", main = "Family CD")
hist(colSums(!is.na(reshape_LO_df_num[,DC_cols]), na.rm = T), breaks = 75, ylab = "Number of F2 individuals", xlab = "Number of successfully called 100kb windows", main = "Family DC")

plot(x = -log10(QTL_results_df$p_val), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Sex")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$p_val_w), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariates: Sex & Weight")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

#current_chr <- 7
#plot(y = -log10(QTL_results_df$p_val_w[QTL_results_df$chr == current_chr]), x =  QTL_results_df$win_start[QTL_results_df$chr == current_chr], pch = 16, xlab = "Position (Chr 7)", ylab = "-log10(P)", main = "Trait: Lifespan; Covariates: Sex & Weight")
#abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$female_no_w_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Subset: Females")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$male_no_w_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Subset: Males")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$female_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Females")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$male_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Males")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$weight_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Covariate: Sex")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$weight_female_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Subset: Females")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$weight_male_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Subset: Males")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

#current_chr <- 4
#plot(y = -log10(QTL_results_df$female_p[QTL_results_df$chr == current_chr]), x =  QTL_results_df$win_start[QTL_results_df$chr == current_chr], pch = 16, xlab = paste0("Position (Chr ", current_chr, ")"), ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight", ylim = c(0,6))
#abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)
#points(y = -log10(QTL_results_df$male_p[QTL_results_df$chr == current_chr]), x =  QTL_results_df$win_start[QTL_results_df$chr == current_chr], pch = 16, col = "red")
#legend(x = "topright", legend = c("Females", "Males"), pch = 16, col = c("black", "red"))

hist(QTL_results_df$est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Covariate: Sex", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = median(QTL_results_df$est), col = "darkorchid", lwd = 2)

hist(QTL_results_df$est_w, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Covariates: Sex & Weight", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = median(QTL_results_df$est_w), col = "darkorchid", lwd = 2)

hist(QTL_results_df$female_no_w_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Subset: Females", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = median(QTL_results_df$female_no_w_est), col = "darkorchid", lwd = 2)

hist(QTL_results_df$male_no_w_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Subset: Males", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = median(QTL_results_df$male_no_w_est), col = "darkorchid", lwd = 2)


hist(QTL_results_df$female_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Covariate: Weight; Subset: Females", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = median(QTL_results_df$female_est), col = "darkorchid", lwd = 2)

hist(QTL_results_df$male_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Covariate: Weight; Subset: Males", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = median(QTL_results_df$male_est), col = "darkorchid", lwd = 2)

dev.off()

##PCA (plink2)

F2_PCA <- read.table(file = "~/Projects/Beetle/data/offspring_LO_calling/plink_PCA/plink_PCA.eigenvec", sep = "\t", header = F, stringsAsFactors = F)
names(F2_PCA) <- c("ID", paste0("PC", 1:10))
F2_PCA_eigenvals <- read.table(file = "~/Projects/Beetle/data/offspring_LO_calling/plink_PCA/plink_PCA.eigenval", sep = "\t", header = F, stringsAsFactors = F)

F2_PCA$plot_col <- NA
F2_PCA$plot_col[grep("^AB", F2_PCA$ID)] <- "firebrick1"
F2_PCA$plot_col[grep("^CD", F2_PCA$ID)] <- "slateblue1"
F2_PCA$plot_col[grep("^BA", F2_PCA$ID)] <- "firebrick4"
F2_PCA$plot_col[grep("^DC", F2_PCA$ID)] <- "slateblue4"

pdf(file = "~/Projects/Beetle/doc/F2_PCA.pdf")

plot(x = F2_PCA$PC1, F2_PCA$PC2, col = F2_PCA$plot_col, xlab = "PC1 (36%)", ylab = "PC2 (9%)", main = "All F2", pch = 20)
legend(x= "top", legend = c("AB","BA", "CD", "DC"), col = c("firebrick1", "firebrick4","slateblue1", "slateblue4"), pch = 20)
plot(x = F2_PCA$PC1[grep("firebrick", F2_PCA$plot_col)], F2_PCA$PC2[grep("firebrick", F2_PCA$plot_col)], col = F2_PCA$plot_col[grep("firebrick", F2_PCA$plot_col)], xlab = "PC1 (36%)", ylab = "PC2 (9%)", main = "Families AB & BA", pch = 20)
legend(x= "top", legend = c("AB", "BA"), col = c("firebrick1", "firebrick4"), pch = 20)
plot(x = F2_PCA$PC1[grep("slateblue", F2_PCA$plot_col)], F2_PCA$PC2[grep("slateblue", F2_PCA$plot_col)], col = F2_PCA$plot_col[grep("slateblue", F2_PCA$plot_col)], xlab = "PC1 (36%)", ylab = "PC2 (9%)", main = "Families CD & DC", pch = 20)
legend(x= "top", legend = c("CD", "DC"), col = c("slateblue1", "slateblue4"), pch = 20)

dev.off()






###########

##Revision following detection fo sample mixup/mislabelling.

# Sample Id	          Sex	EE	Pair	Correct Sample Id
#Sample_VA-3193-EL1-F	F	  L	  D	    Sample_VA-3193-EL1-F
#Sample_VA-3193-EL1-M	M	  E	  D	    Sample_VA-3193-EL1-M
#Sample_VA-3193-EL4-F	F	  L	  B	    Sample_VA-3193-EL4-F
#Sample_VA-3193-EL4-M	M	  E	  B	    Sample_VA-3193-EL4-M
#Sample_VA-3193-LE1-F	M	  L	  A	    Sample_VA-3193-LE1-M <--- These are the
#Sample_VA-3193-LE1-M	F	  E	  A	    Sample_VA-3193-LE1-F <--- mixed-up samples
#Sample_VA-3193-LE2-F	F	  E	  C	    Sample_VA-3193-LE2-F
#Sample_VA-3193-LE2-M	M	  L	  C	    Sample_VA-3193-LE2-M


#Generation the corrected Pedigree

founder_pedigree_corr <- data.frame(ID = unique(c(raw_phenotypes$p1_sire, raw_phenotypes$p1_dam)), Sire = 0, Dam = 0, Sex = NA, Weight = NA, Lifespan = NA, stringsAsFactors = F)
founder_pedigree_corr$Sex[grep("-M$", founder_pedigree_corr$ID)] <- 1
founder_pedigree_corr$Sex[grep("-F$", founder_pedigree_corr$ID)] <- 2

founder_pedigree_corr$Sex[founder_pedigree_corr$ID == "Sample_VA-3193-LE1-M"] <- 2
founder_pedigree_corr$Sex[founder_pedigree_corr$ID == "Sample_VA-3193-LE1-F"] <- 1

#F1_IDs <- unique(c(paste(F2_pedigree$Family, F2_pedigree$parentage_1, sep = "_"),
#                   paste(F2_pedigree$Family, F2_pedigree$parentage_2, sep = "_")))

#F1_IDs <- paste(rep(F1_IDs, each = 2), rep(c("sire", "dam"), times = length(F1_IDs)), sep = "_")

F1_pedigree_corr <- data.frame(ID = F1_IDs, Sire = NA, Dam = NA, Sex = NA, Weight = NA, Lifespan = NA, stringsAsFactors = F)
F1_pedigree_corr$Sex[grep("sire", F1_pedigree_corr$ID)] <- 1
F1_pedigree_corr$Sex[grep("dam", F1_pedigree_corr$ID)] <- 2

F1_pedigree_corr$Sire <- sapply(sub("([A-Z0-9]+)_([LE0-9]+)_(sire|dam)","\\2",F1_pedigree_corr$ID), FUN = grep, x = founder_pedigree_corr$ID[founder_pedigree_corr$Sex == 1], value = T)
F1_pedigree_corr$Dam <- sapply(sub("([A-Z0-9]+)_([LE0-9]+)_(sire|dam)","\\2",F1_pedigree_corr$ID), FUN = grep, x = founder_pedigree_corr$ID[founder_pedigree_corr$Sex == 2], value = T)

#F1_pedigree$Sire[F1_pedigree$Sire == "Sample_VA-3193-LE1-M"] <- "Sample_VA-3193-LE1-F"
#F1_pedigree$Sire[F1_pedigree$Dam == "Sample_VA-3193-LE1-F"] <- "Sample_VA-3193-LE1-M"

###NB! Line origin of the founders must also be corrected!
output_pedigree_corrected <- rbind(founder_pedigree_corr, F1_pedigree_corr, F2_pedigree_out)

#write(x = output_pedigree$ID, file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_samples.txt", ncolumns = 1)
#write.table(x = output_pedigree, file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_pedigree.txt", sep ="\t", quote = F, row.names = F)

#Numeric version
output_pedigree_corrected$numeric_ID <- 1:dim(output_pedigree_corrected)[1]
output_pedigree_corrected$numeric_Sire <- output_pedigree_corrected$numeric_ID[match(output_pedigree_corrected$Sire, output_pedigree_corrected$ID)]
output_pedigree_corrected$numeric_Dam <- output_pedigree_corrected$numeric_ID[match(output_pedigree_corrected$Dam, output_pedigree_corrected$ID)]
output_pedigree_corrected$numeric_Sire[is.na(output_pedigree_corrected$numeric_Sire)] <- 0
output_pedigree_corrected$numeric_Dam[is.na(output_pedigree_corrected$numeric_Dam)] <- 0
output_pedigree_corrected$generation <- 2
output_pedigree_corrected$generation[grep("Sample_VA", output_pedigree_corrected$Sire)] <- 1
output_pedigree_corrected$generation[output_pedigree_corrected$Sire == 0] <- 0
#write.table(x = output_pedigree_corrected[,c("numeric_ID", "numeric_Sire", "numeric_Dam", "Sex", "Weight", "Lifespan")], file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_pedigree_int.txt", sep ="\t", quote = F, row.names = F, col.names = c("ID", "Sire", "Dam", "Sex", "Weight", "Lifespan" ))

output_pedigree_corrected$mat_Grand_Dam <- output_pedigree_corrected$Dam[match(output_pedigree_corrected$Dam, output_pedigree_corrected$ID)]
output_pedigree_corrected$mat_Grand_Sire <- output_pedigree_corrected$Sire[match(output_pedigree_corrected$Dam, output_pedigree_corrected$ID)]
output_pedigree_corrected$pat_Grand_Dam <- output_pedigree_corrected$Dam[match(output_pedigree_corrected$Sire, output_pedigree_corrected$ID)]
output_pedigree_corrected$pat_Grand_Sire <- output_pedigree_corrected$Sire[match(output_pedigree_corrected$Sire, output_pedigree_corrected$ID)]
write.table(x = output_pedigree_corrected, file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_pedigree_corrected_v1.txt", sep ="\t", quote = F, row.names = F, col.names = T)

####
#LO_pedigree <- read.table("/proj/snic2020-6-128/private/a_obtectus_QTLmap/users/mats/custom_LO_calling/A_ob_pedigree_corrected_v1.txt", sep = "\t", stringsAsFactors = F, header = T)
LO_pedigree <- read.table("~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_pedigree_corrected_v1.txt", sep = "\t", stringsAsFactors = F, header = T)
founder_LO_df <- data.frame(ID = LO_pedigree$ID[LO_pedigree$generation == 0], line = NA)
founder_LO_df$line[grep("[0-9]-M", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\1",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])
founder_LO_df$line[grep("[0-9]-F", founder_LO_df$ID)] <- sub(".+-([EL])([EL])[0-9].+","\\2",founder_LO_df$ID[grep("[0-9]-M", founder_LO_df$ID)])
##Correcting for mixed-up samples
founder_LO_df$line[founder_LO_df$ID == "Sample_VA-3193-LE1-M"] <- "E"
founder_LO_df$line[founder_LO_df$ID == "Sample_VA-3193-LE1-F"] <- "L"
##

####

###Some alternative models, including interactions
head(high_cov_windows)
QTL_results_df[QTL_results_df$chr == 7 & QTL_results_df$female_no_w_p < 1e-5,]
target_window <- high_cov_windows["7:5106321-5206320.1"]
#QTL_results_df[QTL_results_df$chr == 7 & QTL_results_df$male_no_w_p < 1e-5,]
#target_window <- high_cov_windows["7:50479594-50579593.1"]


tmp_win_L0s <- reshape_LO_df_num[target_window,grep("call[.]", colnames(reshape_LO_df_num))]
tmp_QTL_df <- F2_pedigree_out
tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismathc detected!")

current_model_int <- glm(Lifespan~LO_call*Sex + 1, data = tmp_QTL_df, na.action = na.exclude)
summary(current_model_int)
current_model <- glm(Lifespan~LO_call + Sex + 1, data = tmp_QTL_df, na.action = na.exclude)
summary(current_model)

pdf(file = "~/Projects/Beetle/doc/F2_QTL_chr7.pdf")

plot(x = -log10(QTL_results_df$male_no_w_p)[QTL_results_df$chr %in% 1:10], y = -log10(QTL_results_df$female_no_w_p)[QTL_results_df$chr %in% 1:10], col = rainbow(10)[QTL_results_df$chr[QTL_results_df$chr %in% 1:10]], pch = 20, xlim = c(0,9), xlab = "-log10P (males)",  ylab = "-log10P (females)", main = "Chromosomes")
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

plot(x = -log10(QTL_results_df$male_no_w_p)[QTL_results_df$chr > 10], y = -log10(QTL_results_df$female_no_w_p)[QTL_results_df$chr > 10], col = (QTL_results_df$chr[QTL_results_df$chr > 10] %% 3) + 1, pch = 20, xlim = c(0,9), xlab = "-log10P (males)",  ylab = "-log10P (females)", main = "Unplaced scaffolds")

current_chr <- 7
plot(y = -log10(QTL_results_df$female_no_w_p[QTL_results_df$chr == current_chr]), x =  QTL_results_df$win_start[QTL_results_df$chr == current_chr], pch = 16, xlab = "Position (Chr 7)", ylab = "-log10(P)", main = "Trait: Lifespan", col = "firebrick") 
points(y = -log10(QTL_results_df$male_no_w_p[QTL_results_df$chr == current_chr]), x =  QTL_results_df$win_start[QTL_results_df$chr == current_chr], pch = 16, col = "slateblue") 
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)
abline(v = 5206320, col = "red", lwd = 2)

legend(x = "topright", col = c("firebrick", "slateblue") , legend = c("Females", "Males"), pch = 20)

boxplot(Lifespan~LO_call + Sex + 1, data = tmp_QTL_df)

dev.off()

#Re-estimation of observed sex
#find /proj/snic2020-6-128/private/a_obtectus_QTLmap/offspring/data/03-bam-files/split_per_sample/ -maxdepth 1 -type f -name '*.bam' > ./offspring_bam_list.txt

x_cov_file_list <- dir("~/Projects/Beetle/data/offspring_sex_estimation/cov_X/")
x_cov_file_list <- x_cov_file_list[grepl("_X_cov.txt", x_cov_file_list)]
y_cov_file_list <- dir("~/Projects/Beetle/data/offspring_sex_estimation/cov_Y/")
y_cov_file_list <- y_cov_file_list[grepl("_Y_cov.txt", y_cov_file_list)]

cov_sample_list <- sub("_X_cov.txt", "", x_cov_file_list)

offspring_cov_df <- data.frame(ID = cov_sample_list, x_cov = NA, y_cov = NA, stringsAsFactors = F)

for(cov_s in offspring_cov_df$ID){
  tmp_x_cov_file <- paste0("~/Projects/Beetle/data/offspring_sex_estimation/cov_X/", grep(cov_s, x_cov_file_list, value = T))
  tmp_x_cov <- try(read.table(tmp_x_cov_file, sep = "\t", header = F, stringsAsFactors = F))
  tmp_y_cov_file <- paste0("~/Projects/Beetle/data/offspring_sex_estimation/cov_Y/", grep(cov_s, y_cov_file_list, value = T))
  tmp_y_cov <- try(read.table(tmp_y_cov_file, sep = "\t", header = F, stringsAsFactors = F))
  if(grepl("try-error", class(tmp_x_cov))){
    offspring_cov_df$x_cov[offspring_cov_df$ID == cov_s] <- NA
  } else{
    offspring_cov_df$x_cov[offspring_cov_df$ID == cov_s] <- mean(tmp_x_cov$V7[1], na.rm  = T)
  }
  if(grepl("try-error", class(tmp_y_cov))){
    offspring_cov_df$y_cov[offspring_cov_df$ID == cov_s] <- NA
  } else{
    offspring_cov_df$y_cov[offspring_cov_df$ID == cov_s] <- mean(tmp_y_cov$V7, na.rm  = T)
  }
}

offspring_cov_df$cov_ratio <- offspring_cov_df$y_cov/offspring_cov_df$x_cov

hist(offspring_cov_df$cov_ratio, breaks = 100)

offspring_cov_df$est_sex <- NA
offspring_cov_df$est_sex[which(offspring_cov_df$cov_ratio > 0.8)] <- 1
offspring_cov_df$est_sex[which(offspring_cov_df$cov_ratio < 0.5)] <- 2
offspring_cov_df$phenotypic_sex <- NA
offspring_cov_df$phenotypic_sex <- F2_pedigree_out$Sex[match(offspring_cov_df$ID, F2_pedigree_out$ID)]
offspring_cov_df$pedigree_ID <- NA
offspring_cov_df$pedigree_ID <- F2_pedigree_out$ID[match(offspring_cov_df$ID, F2_pedigree_out$ID)]

F2_pedigree_out$cov_est_sex <- NA
F2_pedigree_out$cov_est_sex <- offspring_cov_df$est_sex[match(F2_pedigree_out$ID, offspring_cov_df$ID)]

QTL_results_df$female_COV_p <- NA
QTL_results_df$female_COV_est <- NA
QTL_results_df$male_COV_p <- NA
QTL_results_df$male_COV_est <- NA

for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismatch detected!")
  
  current_model_male <- glm(Lifespan~LO_call + Weight + 1, data = tmp_QTL_df[tmp_QTL_df$cov_est_sex == 1,], na.action = na.exclude)
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "male_COV_est"] <- currentLO_e_male
  QTL_results_df[i, "male_COV_p"] <- currentLO_p_male
  
  current_model_female <- glm(Lifespan~LO_call + Weight + 1, data = tmp_QTL_df[tmp_QTL_df$cov_est_sex == 2,], na.action = na.exclude)
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "female_COV_est"] <- currentLO_e_female
  QTL_results_df[i, "female_COV_p"] <- currentLO_p_female
}

QTL_results_df$male_COV_no_w_p <- NA
QTL_results_df$male_COV_no_w_est <- NA

QTL_results_df$female_COV_no_w_p <- NA
QTL_results_df$female_COV_no_w_est <- NA

for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismatch detected!")
  
  current_model_male <- glm(Lifespan~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$cov_est_sex == 1,], na.action = na.exclude)
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "male_COV_no_w_est"] <- currentLO_e_male
  QTL_results_df[i, "male_COV_no_w_p"] <- currentLO_p_male
  
  current_model_female <- glm(Lifespan~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$cov_est_sex == 2,], na.action = na.exclude)
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "female_COV_no_w_est"] <- currentLO_e_female
  QTL_results_df[i, "female_COV_no_w_p"] <- currentLO_p_female
}


QTL_results_df$COV_w_p <- NA
QTL_results_df$COV_w_est <- NA

QTL_results_df$COV_no_w_p <- NA
QTL_results_df$COV_no_w_est <- NA

for(i in 1:length(high_cov_windows)){
  tmp_win_L0s <- reshape_LO_df_num[high_cov_windows[i],grep("call[.]", colnames(reshape_LO_df_num))]
  tmp_QTL_df <- F2_pedigree_out
  tmp_QTL_df$LO_call <- as.integer(tmp_win_L0s)
  tmp_QTL_df$LO_ID <- colnames(tmp_win_L0s)
  if(sum(sub("call[.]", "", tmp_QTL_df$LO_ID) != tmp_QTL_df$ID) > 0) print("ID mismatch detected!")
  
  current_model <- glm(Lifespan~LO_call + cov_est_sex + 1, data = tmp_QTL_df, na.action = na.exclude)
  currentLO_p_val <- summary(current_model)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_estimate <- summary(current_model)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "COV_no_w_est"] <- currentLO_estimate
  QTL_results_df[i, "COV_no_w_p"] <- currentLO_p_val
  
  current_model_w <- glm(Lifespan~LO_call + cov_est_sex + Weight + 1, data = tmp_QTL_df, na.action = na.exclude)
  currentLO_p_val_w  <- summary(current_model_w)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_estimate_w  <- summary(current_model_w)$coefficients["LO_call","Estimate"]
  QTL_results_df[i, "COV_w_est"] <- currentLO_estimate_w 
  QTL_results_df[i, "COV_w_p"] <- currentLO_p_val_w 
}

#plot(x = -log10(QTL_results_df$female_COV_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Females (based on coverage)")
#abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

#plot(x = -log10(QTL_results_df$male_COV_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Males (based on coverage)")
#abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

pdf(file = "~/Projects/Beetle/doc/F2_sex_re_estimation.pdf")

hist(F2_pedigree_out$Lifespan, xlab = "Lifespan (days)", breaks = 50, xlim = c(0,60))
hist(F2_pedigree_out$Lifespan[F2_pedigree_out$cov_est_sex == 1], xlab = "Lifespan (days)", breaks = 50, main = "", xlim = c(0,60), col = rgb(1,0,0,0.5))
hist(F2_pedigree_out$Lifespan[F2_pedigree_out$cov_est_sex == 2], xlab = "Lifespan (days)", breaks = 50, xlim = c(0,60), col = rgb(0,0,1,0.5), add = T)
legend(legend = c("Females", "Males"), fill = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)), x = "topright")

plot(x = F2_pedigree_out$Lifespan[F2_pedigree_out$cov_est_sex == 2],y =  F2_pedigree_out$Weight[F2_pedigree_out$cov_est_sex == 2], xlab = "Lifespan (days)", xlim = c(0,60), col = rgb(0,0,1,0.5), pch = 16, ylab = "Weight (g)")
points(x = F2_pedigree_out$Lifespan[F2_pedigree_out$cov_est_sex == 1],y =  F2_pedigree_out$Weight[F2_pedigree_out$cov_est_sex == 1], xlim = c(0,60), col = rgb(1,0,0,0.5), pch = 16)
legend(legend = c("Females", "Males"), fill = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5)), x = "bottomright")


plot(x = -log10(QTL_results_df$female_COV_p), y = -log10(QTL_results_df$female_p),col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "P-values - Coverage sexing", ylab = "P-values - Phenotypic sexing", main = "Females")
plot(x = -log10(QTL_results_df$male_COV_p), y = -log10(QTL_results_df$male_p),col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "P-values - Coverage sexing", ylab = "P-values - Phenotypic sexing", main = "Males")
plot(x = -log10(QTL_results_df$female_COV_no_w_p), y = -log10(QTL_results_df$female_no_w_p),col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "P-values - Coverage sexing", ylab = "P-values - Phenotypic sexing", main = "Females (no weight correction)")
plot(x = -log10(QTL_results_df$male_COV_no_w_p), y = -log10(QTL_results_df$male_no_w_p),col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "P-values - Coverage sexing", ylab = "P-values - Phenotypic sexing", main = "Males (no weight correction)")
plot(x = -log10(QTL_results_df$COV_no_w_p), y = -log10(QTL_results_df$p_val),col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "P-values - Coverage sexing", ylab = "P-values - Phenotypic sexing", main = "All (no weight correction)")
plot(x = -log10(QTL_results_df$COV_w_p), y = -log10(QTL_results_df$p_val_w),col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "P-values - Coverage sexing", ylab = "P-values - Phenotypic sexing", main = "All (weight correction)")

plot(x = -log10(QTL_results_df$COV_no_w_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Sex; Subset: All")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$COV_w_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariates: Sex, Weight; Subset: All")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)


plot(x = -log10(QTL_results_df$female_COV_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Females (based on coverage)")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$male_COV_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Males (based on coverage)")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$female_COV_no_w_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Subset: Females (based on coverage)")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(QTL_results_df$male_COV_no_w_p), col = (QTL_results_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Subset: Males (based on coverage)")
abline(h = -log10(0.05/length(QTL_results_df$win_start)), col = "darkorchid", lwd = 2)

hist(offspring_cov_df$cov_ratio, breaks = 100, xlab = "y_coverage/x_coverage")

hist(QTL_results_df$female_COV_no_w_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Subset: Females (based on coverage)", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(QTL_results_df$female_COV_no_w_est), col = "darkorchid", lwd = 2)
abline(v =t.test(QTL_results_df$female_COV_no_w_est,alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")

hist(QTL_results_df$male_COV_no_w_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Subset: Males (based on coverage)", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(QTL_results_df$male_COV_no_w_est), col = "darkorchid", lwd = 2)
abline(v =t.test(QTL_results_df$male_COV_no_w_est,alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")

hist(QTL_results_df$female_COV_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Covariate: Weight; Subset: Females (based on coverage)", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(QTL_results_df$female_COV_est), col = "darkorchid", lwd = 2)
abline(v =t.test(QTL_results_df$female_COV_est,alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")

hist(QTL_results_df$male_COV_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Covariate: Weight; Subset: Males (based on coverage)", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(QTL_results_df$male_COV_est), col = "darkorchid", lwd = 2)
abline(v =t.test(QTL_results_df$male_COV_est,alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")

plot(x = -log10(QTL_results_df$male_COV_p)[QTL_results_df$chr %in% 1:10], y = -log10(QTL_results_df$female_COV_p)[QTL_results_df$chr %in% 1:10], col = rainbow(10)[QTL_results_df$chr[QTL_results_df$chr %in% 1:10]], pch = 20, xlim = c(0,7), ylim = c(0,7), xlab = "-log10P Males (coverage)",  ylab = "-log10P Females (coverage)", main = "Chromosomes")
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

dev.off()

write.table(x = F2_pedigree_out, file = "~/Projects/Beetle/data/offspring_LO_calling/pedigree_info/A_ob_F2_sex_re_estimate.txt", sep = "\t", row.names = F, quote = F)

EL_diff_reg <- readRDS("~/Projects/Beetle/data/EL_pool_seq/Significant-genomicRegions_OTT_16_pools.dAF_EI.II.III.IV_vs_LV.VI.VII.VIII_dAFgt0.9_20210402_sig0.9.RDS")
pdf(file = "~/Projects/Beetle/doc/Chr7_QTL_vs_divergence.pdf")
plot(y = -log10(QTL_results_df$female_COV_p[QTL_results_df$chr == 7]), x = QTL_results_df$win_start[QTL_results_df$chr == 7], pch = 16, xlab = "Position (Chr 7)", ylab = "-log10(P)",  main = "Trait: Lifespan; Covariate: Weight; Subset: Females (based on coverage)") 
abline(v = c(1.5e7, 3e7), lwd = 2, col = "red")
boxplot(x = list(a = -log10(QTL_results_df$female_COV_p[QTL_results_df$chr == 7 & QTL_results_df$win_start < 1.5e7]), 
                 b = -log10(QTL_results_df$female_COV_p[QTL_results_df$chr == 7 & QTL_results_df$win_start >= 1.5e7 & QTL_results_df$win_start < 3e7]),
                 c = -log10(QTL_results_df$female_COV_p[QTL_results_df$chr == 7 & QTL_results_df$win_start >= 3e7])), 
        ylim = c(0,7), main = "Trait: Lifespan; Covariate: Weight; Subset: Females (based on coverage)", ylab = "-log10(P)") 
text(x = c(1,2,3), y = c(6.5, 6.5, 6.5), paste("Div. frac:", round(x = c(sum(EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$width[EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$seqnames == "chr_7" & EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$start < 1.5e7])/1.5e7, 
                                                                          sum(EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$width[EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$seqnames == "chr_7" & EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$start >= 1.5e7 & EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$start < 3e7])/1.5e7,
                                                                          sum(EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$width[EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$seqnames == "chr_7" & EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$start >= 3e7])/(max(EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$start[EL_diff_reg$dAF_EI.II.III.IV_vs_LV.VI.VII.VIII$seqnames == "chr_7"]) - 3e7)), 
                                                              digits = 3)))

dev.off()

##Using imputed & corrected LO calls


imp_LO_calls <- read.csv(file = "~/Projects/Beetle/data/offspring_LO_calling/imputed_LOs/F2_LO_calls_100kb_imputed_filtered_all_chromosomes.csv", sep = ",")

adj_LO_names <- sub("X([0-9]+)[.]([0-9]+)[.]([0-9]+)[.]1", "\\1:\\2-\\3", names(imp_LO_calls))
names(imp_LO_calls) <- adj_LO_names

imp_QTL_df <- data.frame(win_ID = names(imp_LO_calls[-c(1:7)]), chr = NA, win_start = NA,  p_val = NA, est = NA, p_val_w = NA, est_w = NA)

imp_QTL_df$chr <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\1", imp_QTL_df$win_ID))
imp_QTL_df$win_start <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\2", imp_QTL_df$win_ID))


for(i in 1:dim(imp_QTL_df)[1]){
  tmp_QTL_df <- imp_LO_calls[,c(1:7,i+7)]
  names(tmp_QTL_df)[8] <- "LO_call"
 
  current_model <- glm(Lifespan~LO_call + Sex_coverage + 1, data = tmp_QTL_df, na.action = na.exclude)
  currentLO_p_val <- summary(current_model)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_estimate <- summary(current_model)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "est"] <- currentLO_estimate
  imp_QTL_df[i, "p_val"] <- currentLO_p_val
  
  current_model_w <- glm(Lifespan~LO_call + Sex_coverage + Weight + 1, data = tmp_QTL_df, na.action = na.exclude)
  currentLO_p_val_w  <- summary(current_model_w)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_estimate_w  <- summary(current_model_w)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "est_w"] <- currentLO_estimate_w 
  imp_QTL_df[i, "p_val_w"] <- currentLO_p_val_w 
}

imp_QTL_df$male_no_w_p <- NA
imp_QTL_df$male_no_w_est <- NA

imp_QTL_df$female_no_w_p <- NA
imp_QTL_df$female_no_w_est <- NA

for(i in 1:dim(imp_QTL_df)[1]){
  tmp_QTL_df <- imp_LO_calls[,c(1:7,i+7)]
  names(tmp_QTL_df)[8] <- "LO_call"
  
  current_model_male <- try(glm(Lifespan~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "M",], na.action = na.exclude))
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "male_no_w_est"] <- currentLO_e_male
  imp_QTL_df[i, "male_no_w_p"] <- currentLO_p_male
  
  current_model_female <- glm(Lifespan~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "F",], na.action = na.exclude)
  
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "female_no_w_est"] <- currentLO_e_female
  imp_QTL_df[i, "female_no_w_p"] <- currentLO_p_female
}

imp_QTL_df$male_p <- NA
imp_QTL_df$male_est <- NA

imp_QTL_df$female_p <- NA
imp_QTL_df$female_est <- NA

for(i in 1:dim(imp_QTL_df)[1]){
  tmp_QTL_df <- imp_LO_calls[,c(1:7,i+7)]
  names(tmp_QTL_df)[8] <- "LO_call"
  
  current_model_male <- try(glm(Lifespan~LO_call + Weight + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "M",], na.action = na.exclude))
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "male_est"] <- currentLO_e_male
  imp_QTL_df[i, "male_p"] <- currentLO_p_male
  
  current_model_female <- glm(Lifespan~LO_call + Weight + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "F",], na.action = na.exclude)
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "female_est"] <- currentLO_e_female
  imp_QTL_df[i, "female_p"] <- currentLO_p_female
}

imp_QTL_df$male_weight_p <- NA
imp_QTL_df$male_weight_est <- NA

imp_QTL_df$female_weight_p <- NA
imp_QTL_df$female_weight_est <- NA

for(i in 1:dim(imp_QTL_df)[1]){
  tmp_QTL_df <- imp_LO_calls[,c(1:7,i+7)]
  names(tmp_QTL_df)[8] <- "LO_call"
  
  current_model_male <- try(glm(Weight~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "M",], na.action = na.exclude))
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "male_weight_est"] <- currentLO_e_male
  imp_QTL_df[i, "male_weight_p"] <- currentLO_p_male
  
  current_model_female <- glm(Weight~LO_call + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "F",], na.action = na.exclude)
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "female_weight_est"] <- currentLO_e_female
  imp_QTL_df[i, "female_weight_p"] <- currentLO_p_female
}

imp_QTL_df$male_weight_w_life_p <- NA
imp_QTL_df$male_weight_w_life_est <- NA

imp_QTL_df$female_weight_w_life_p <- NA
imp_QTL_df$female_weight_w_life_est <- NA

for(i in 1:dim(imp_QTL_df)[1]){
  tmp_QTL_df <- imp_LO_calls[,c(1:7,i+7)]
  names(tmp_QTL_df)[8] <- "LO_call"
  
  current_model_male <- try(glm(Weight~LO_call + Lifespan + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "M",], na.action = na.exclude))
  currentLO_p_male <- summary(current_model_male)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_male <- summary(current_model_male)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "male_weight_w_life_est"] <- currentLO_e_male
  imp_QTL_df[i, "male_weight_w_life_p"] <- currentLO_p_male
  
  current_model_female <- glm(Weight~LO_call + Lifespan + 1, data = tmp_QTL_df[tmp_QTL_df$Sex_coverage == "F",], na.action = na.exclude)
  currentLO_p_female <- summary(current_model_female)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_e_female <- summary(current_model_female)$coefficients["LO_call","Estimate"]
  imp_QTL_df[i, "female_weight_w_life_est"] <- currentLO_e_female
  imp_QTL_df[i, "female_weight_w_life_p"] <- currentLO_p_female
}

#Sex-by-genotype interactions
imp_QTL_df$sex_by_geno_p <- NA
imp_QTL_df$sex_by_geno_est <- NA

for(i in 1:dim(imp_QTL_df)[1]){
  tmp_QTL_df <- imp_LO_calls[,c(1:7,i+7)]
  names(tmp_QTL_df)[8] <- "LO_call"
  
  current_model_int <- try(glm(Lifespan~LO_call*Sex_coverage + 1, data = tmp_QTL_df, na.action = na.exclude))
  currentLO_p_int <- summary(current_model_int)$coefficients["LO_call:Sex_coverageM","Pr(>|t|)"]
  currentLO_e_int <- summary(current_model_int)$coefficients["LO_call:Sex_coverageM","Estimate"]
  imp_QTL_df[i, "sex_by_geno_est"] <- currentLO_e_int
  imp_QTL_df[i, "sex_by_geno_p"] <- currentLO_p_int
}

plot(x = -log10(imp_QTL_df$sex_by_geno_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; GxS interaction")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)
plot(x = imp_QTL_df$sex_by_geno_est, col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "", main = "Trait: Lifespan; GxS interaction")




###Peak calling

bonferroni_signif <- 0.05/length(imp_QTL_df$win_start)

overall_peaks_no_weight <- beetle_QTL_peaks(qtl_df = imp_QTL_df[imp_QTL_df$chr %in% 1:10,],
                                            target_p_col = "p_val", target_est_col = "est",
                                            peak_open = bonferroni_signif*5, peak_retain = bonferroni_signif, 
                                            drop_off = 1e2)

female_peaks_no_weight <- beetle_QTL_peaks(qtl_df = imp_QTL_df[imp_QTL_df$chr %in% 1:10,],
                                            target_p_col = "female_no_w_p", target_est_col = "female_no_w_est",
                                            peak_open = bonferroni_signif*5, peak_retain = bonferroni_signif, 
                                            drop_off = 1e2)

male_peaks_no_weight <- beetle_QTL_peaks(qtl_df = imp_QTL_df[imp_QTL_df$chr %in% 1:10,],
                                            target_p_col = "male_no_w_p", target_est_col = "male_no_w_est",
                                            peak_open = bonferroni_signif*5, peak_retain = bonferroni_signif, 
                                            drop_off = 1e2)

beetle_QTL_peaks <- function(qtl_df, target_p_col, target_est_col, peak_open = 5e5, peak_retain = 1e5, drop_off = 1e3){
  peak_df_out <- data.frame(chr = NA, start = NA, prel_start = NA,  stop = NA, peak_pos = NA, peak_est = NA, min_p_val = NA, peak_index = NA)
  
  peak_idx_glob <- 0
  for(peak_chr in unique(qtl_df$chr)){
    peak_idx_chr <- 0
    tmp_peak_df <- qtl_df[qtl_df$chr == peak_chr, c("chr", "win_start", target_p_col, target_est_col)]
    in_peak <- F
    for(i in 1:dim(tmp_peak_df)[1]){
      if(!in_peak){
        if(tmp_peak_df[i,target_p_col] < peak_open){
          peak_idx_glob <- peak_idx_glob + 1
          peak_idx_chr <- peak_idx_chr + 1
          print(paste0("Preliminary peak started. Glob index: ", peak_idx_glob, "; Chr index: " , peak_idx_chr, "; Window: ", i))
          peak_df_out[peak_idx_glob,"chr"] <- peak_chr
          peak_df_out[peak_idx_glob,"prel_start"] <- tmp_peak_df[i,"win_start"]
          peak_df_out[peak_idx_glob,"peak_index"] <- peak_idx_glob
          peak_df_out[peak_idx_glob,"min_p_val"] <- tmp_peak_df[i,target_p_col]
          peak_df_out[peak_idx_glob,"peak_est"] <- tmp_peak_df[i,target_est_col]
          peak_df_out[peak_idx_glob,"peak_pos"] <- tmp_peak_df[i,"win_start"]
          in_peak <- T
        }
      } else{
        if(tmp_peak_df[i,target_p_col] < peak_df_out[peak_idx_glob,"min_p_val"]){
          peak_df_out[peak_idx_glob,"min_p_val"] <- tmp_peak_df[i,target_p_col]
          peak_df_out[peak_idx_glob,"peak_pos"] <- tmp_peak_df[i,"win_start"]
        }
        if(tmp_peak_df[i,target_p_col] > max(c(peak_open, peak_df_out[peak_idx_glob,"min_p_val"]*drop_off)) | i == dim(tmp_peak_df)[1]){
          peak_df_out[peak_idx_glob,"stop"] <- tmp_peak_df[i,"win_start"]
          
          if(any(tmp_peak_df[1:(i-1),target_p_col] > max(c(peak_open, peak_df_out[peak_idx_glob,"min_p_val"]*drop_off)))){
            tmp_peak_start_win <- max(which(tmp_peak_df[1:(i-1),target_p_col] > max(c(peak_open, peak_df_out[peak_idx_glob,"min_p_val"]*drop_off))))
            peak_df_out[peak_idx_glob,"start"] <- tmp_peak_df[tmp_peak_start_win,"win_start"]
          } else{
            peak_df_out[peak_idx_glob,"start"] <- 1
          }
          
          if(peak_idx_chr > 1){
            peak_df_out[peak_idx_glob,"start"] <- max(c(peak_df_out[peak_idx_glob,"start"],peak_df_out[peak_idx_glob - 1,"stop"]))
          }
          
          print(paste0("Preliminary peak closed. Glob index: ", peak_idx_glob, "; Chr index: ", peak_idx_chr, "; Window: ", i))
          
          if(peak_df_out[peak_idx_glob,"min_p_val"] > peak_retain){
            peak_df_out[peak_idx_glob,] <- NA
            print(paste0("Preliminary peak dropped. Glob index: ", peak_idx_glob, "; Chr index: " , peak_idx_chr, "; Window: ", i))
            peak_idx_glob <- peak_idx_glob - 1
            peak_idx_chr <- peak_idx_chr - 1
          }
          in_peak <- F
        }
      }
    }
  }  
  peak_df_out <- peak_df_out[!is.na(peak_df_out$chr),]
  return(peak_df_out)
}
### 


###Permuations

permuted_LOs <- imp_LO_calls
per_vec <- sample(dim(permuted_LOs)[1])
permuted_LOs[,1:7] <- imp_LO_calls[per_vec, 1:7]

permutation_df <- data.frame(win_ID = names(permuted_LOs[-c(1:7)]), chr = NA, win_start = NA,  p_val = NA, est = NA)

permutation_df$chr <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\1", permutation_df$win_ID))
permutation_df$win_start <- as.integer(sub("([0-9]+)[:]([0-9]+)[-]([0-9]+)", "\\2", permutation_df$win_ID))


for(i in 1:dim(permutation_df)[1]){
  tmp_perm_df <- permuted_LOs[,c(1:7,i+7)]
  names(tmp_perm_df)[8] <- "LO_call"
  
  current_model <- glm(Lifespan~LO_call + Sex_coverage + 1, data = tmp_perm_df, na.action = na.exclude)
  currentLO_p_val <- summary(current_model)$coefficients["LO_call","Pr(>|t|)"]
  currentLO_estimate <- summary(current_model)$coefficients["LO_call","Estimate"]
  permutation_df[i, "est"] <- currentLO_estimate
  permutation_df[i, "p_val"] <- currentLO_p_val
}

plot(x = -log10(permutation_df$p_val), col  = (permutation_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Sex; Subset: All")



##Plotting
pdf(file = "~/Projects/Beetle/doc/F2_LO_imputation.pdf")


plot(x = -log10(imp_QTL_df$p_val), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Sex; Subset: All", ylim = c(0,12))
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)


plot(x = -log10(imp_QTL_df$p_val_w), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariates: Sex, Weight; Subset: All")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$female_no_w_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Subset: Females")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$male_no_w_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Subset: Males")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$female_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Females")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$male_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Lifespan; Covariate: Weight; Subset: Males")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$female_weight_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Subset: Females")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$male_weight_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Subset: Males")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$female_weight_w_life_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Covariate: Lifespan; Subset: Females")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$male_weight_w_life_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Covariate: Lifespan; Subset: Males")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)



hist(imp_QTL_df$female_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan; Covariate: Weight; Subset: Females (based on coverage)", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(imp_QTL_df$female_est), col = "darkorchid", lwd = 2)
abline(v =t.test(imp_QTL_df$female_est, alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")

hist(imp_QTL_df$male_est, breaks = 50, main = "Estimates across all windows; Trait: Lifespan;  Covariate: Weight; Subset: Males (based on coverage)", xlab = "Estimate (days)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(imp_QTL_df$male_est), col = "darkorchid", lwd = 2)
abline(v =t.test(imp_QTL_df$male_est,alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")

plot(x = -log10(imp_QTL_df$male_no_w_p)[imp_QTL_df$chr %in% 1:10], y = -log10(imp_QTL_df$female_no_w_p)[imp_QTL_df$chr %in% 1:10], col = rainbow(10)[imp_QTL_df$chr[imp_QTL_df$chr %in% 1:10]], pch = 20, xlim = c(0,7), ylim = c(0,7), xlab = "-log10P Males (coverage)",  ylab = "-log10P Females (coverage)", main = "Chromosomes")
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

dev.off()

pdf(file = "~/Projects/Beetle/doc/F2_effect_by_sex.pdf")

plot(x = imp_QTL_df$male_no_w_est[imp_QTL_df$chr %in% 1:10], y = imp_QTL_df$female_no_w_est[imp_QTL_df$chr %in% 1:10], col = rainbow(10)[imp_QTL_df$chr[imp_QTL_df$chr %in% 1:10]], pch = 20, xlim = c(-2,2), ylim = c(-2,2), xlab = "Effect esitmate (Males)",  ylab = "Effect esitmate (Females)", main = "Estimates across all windows; Trait: Lifespan", axes = F)
axis(1, pos = 0, at = c(-2:2))
axis(2, pos = 0, at = c(-2:2))
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

plot(x = imp_QTL_df$male_est[imp_QTL_df$chr %in% 1:10], y = imp_QTL_df$female_est[imp_QTL_df$chr %in% 1:10], col = rainbow(10)[imp_QTL_df$chr[imp_QTL_df$chr %in% 1:10]], pch = 20, xlim = c(-2,2), ylim = c(-2,2), xlab = "Effect esitmate (Males)",  ylab = "Effect esitmate (Females)", main = "Estimates across all windows; Trait: Lifespan; Covariate: Weight", axes = F)
axis(1, pos = 0, at = c(-2:2))
axis(2, pos = 0, at = c(-2:2))
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

dev.off()


pdf(file = "~/Projects/Beetle/doc/F2_weight_models.pdf")

plot(x = -log10(imp_QTL_df$female_weight_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Subset: Females")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$male_weight_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Subset: Males")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$female_weight_w_life_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Covariate: Lifespan; Subset: Females")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)

plot(x = -log10(imp_QTL_df$male_weight_w_life_p), col = (imp_QTL_df$chr %% 2) + 1, pch = 16, xlab = "", ylab = "-log10(P)", main = "Trait: Weight; Covariate: Lifespan; Subset: Males")
abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "darkorchid", lwd = 2)


plot(x = -log10(imp_QTL_df$male_weight_p)[imp_QTL_df$chr %in% 1:10], y = -log10(imp_QTL_df$female_weight_p)[imp_QTL_df$chr %in% 1:10], xlim = c(0,6), ylim = c(0,6), col = rainbow(10)[imp_QTL_df$chr[imp_QTL_df$chr %in% 1:10]], pch = 20, xlab = "-log10P (Males)",  ylab = "-log10P (Females)", main = "P-values across all windows; Trait: Weight")
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

plot(x = -log10(imp_QTL_df$male_weight_w_life_p)[imp_QTL_df$chr %in% 1:10], y = -log10(imp_QTL_df$female_weight_w_life_p)[imp_QTL_df$chr %in% 1:10], xlim = c(0,6), ylim = c(0,6), col = rainbow(10)[imp_QTL_df$chr[imp_QTL_df$chr %in% 1:10]], pch = 20, xlab = "-log10P (Males)",  ylab = "-log10P (Females)", main = "P-values across all windows; Trait: Weight; Covariate: Lifespan")
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)


plot(x = imp_QTL_df$male_weight_est[imp_QTL_df$chr %in% 1:10], y = imp_QTL_df$female_weight_est[imp_QTL_df$chr %in% 1:10], col = rainbow(10)[imp_QTL_df$chr[imp_QTL_df$chr %in% 1:10]], pch = 20, xlab = "Effect esitmate (Males)",  ylab = "Effect esitmate (Females)", main = "Estimates across all windows; Trait: Weight", axes = F)
axis(1, pos = 0, at = c(-1:2)/10000)
axis(2, pos = 0, at = c(-1:2)/10000)
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

plot(x = imp_QTL_df$male_weight_w_life_est[imp_QTL_df$chr %in% 1:10], y = imp_QTL_df$female_weight_w_life_est[imp_QTL_df$chr %in% 1:10], col = rainbow(10)[imp_QTL_df$chr[imp_QTL_df$chr %in% 1:10]], pch = 20, xlab = "Effect esitmate (Males)",  ylab = "Effect esitmate (Females)", main = "Estimates across all windows; Trait: Weight; Covariate: Lifespan", axes = F)
axis(1, pos = 0, at = c(-1:2)/10000)
axis(2, pos = 0, at = c(-1:2)/10000)
legend(x = "topright", col = rainbow(10)[(1:10 %% 11)] , legend = 1:10, pch = 20)

hist(imp_QTL_df$female_weight_est, breaks = 50, main = "Estimates across all windows; Trait: Weight; Subset: Females", xlab = "Estimate (g)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(imp_QTL_df$female_weight_est), col = "darkorchid", lwd = 2)
abline(v =t.test(imp_QTL_df$female_weight_est, alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")

hist(imp_QTL_df$male_weight_est, breaks = 50, main = "Estimates across all windows; Trait: Weight;  Subset: Males", xlab = "Estimate (g)")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(imp_QTL_df$male_weight_est), col = "darkorchid", lwd = 2)
abline(v =t.test(imp_QTL_df$male_weight_est,alternative = "two")$conf.int, col = "darkorchid", lwd = 1.5, lty = "dashed")
dev.off()

####
pdf(file = "~/Projects/Beetle/doc/F2_QTL_peaks.pdf")

for(qtl_chr in unique(c(female_peaks_no_weight$chr, male_peaks_no_weight$chr))){

  plot(y = -log10(imp_QTL_df$male_no_w_p[imp_QTL_df$chr == qtl_chr]), x = imp_QTL_df$win_start[imp_QTL_df$chr == qtl_chr], col = "olivedrab1", pch = 16, xlab = "", ylab = "-log10(P)", main = paste0("Trait: Lifespan; Covariate: Sex, Chr: ", qtl_chr), ylim = c(0,12))
  points(y = -log10(imp_QTL_df$female_no_w_p[imp_QTL_df$chr == qtl_chr]), x = imp_QTL_df$win_start[imp_QTL_df$chr == qtl_chr], col = "darkorchid1", pch = 16)
  
  
  abline(h = -log10(0.05/length(imp_QTL_df$win_start)), col = "red", lwd = 2)
  n_peaks_m <- sum(male_peaks_no_weight$chr == qtl_chr)
  n_peaks_f <- sum(female_peaks_no_weight$chr == qtl_chr)
  if(n_peaks_m > 0){
    segments(x0 = male_peaks_no_weight$start[male_peaks_no_weight$chr == qtl_chr], 
         x1 = male_peaks_no_weight$stop[male_peaks_no_weight$chr == qtl_chr],
         y0 = 11 + c(0.2,-0.2)[1:n_peaks_m %% 2 + 1], 
         y1 = 11 + c(0.2,-0.2)[1:n_peaks_m %% 2 + 1],
         lwd = 2, col = "olivedrab4")
  }
  if(n_peaks_f > 0){
    segments(x0 = female_peaks_no_weight$start[female_peaks_no_weight$chr == qtl_chr], 
           x1 = female_peaks_no_weight$stop[female_peaks_no_weight$chr == qtl_chr],
           y0 = 10 + c(0.2,-0.2)[1:n_peaks_f %% 2 + 1], 
           y1 = 10 + c(0.2,-0.2)[1:n_peaks_f %% 2 + 1],
           lwd = 2, col = "darkorchid4")
  }
  legend(x = "topright", legend = c("Males", "Females"), col = c("olivedrab1","darkorchid1"), pch = 16)
}

dev.off()


