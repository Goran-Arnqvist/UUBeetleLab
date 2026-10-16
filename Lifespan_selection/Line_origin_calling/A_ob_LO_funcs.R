#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se
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
  eval_pos <- eval_pos[eval_pos %in% rownames(offspring_GT_df)]
  
  if(length(eval_pos) > 1){
    #Evaluating the target offspring
    eval_offspring_GTs <- offspring_GT_df[eval_pos,grep(target_offspring, colnames(offspring_GT_df))]
    eval_L_GTs <- founder_GT_df[eval_pos,L_founder_cols[1]]
    eval_E_GTs <- founder_GT_df[eval_pos,E_founder_cols[1]]
    
    LL_pos <- which(eval_offspring_GTs[,1] == eval_L_GTs & eval_offspring_GTs[,2] == eval_L_GTs)
    EE_pos <- which(eval_offspring_GTs[,1] == eval_E_GTs & eval_offspring_GTs[,2] == eval_E_GTs)
    LE_pos <- which((eval_offspring_GTs[,1] == eval_L_GTs & eval_offspring_GTs[,2] == eval_E_GTs) & (eval_offspring_GTs[,2] == eval_L_GTs & eval_offspring_GTs[,1] == eval_E_GTs))
    out_array <- c(length(LL_pos), length(LE_pos), length(EE_pos))
  } else{
    out_array <- c(0, 0, 0)
  }
  names(out_array) <- c("LL", "LE", "EE")
  return(out_array)
}
