library(dplyr)

FLiP_lib <- read.csv("./OutputData/FLiPR_NewNorm_anova_results.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
FLiP_lib <- FLiP_lib[which(FLiP_lib$q < 0.05), ]

proteins <- read.delim("./FilteredData/FLiPRFractions_tc_data_filtered.tsv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
FLiP_lib <- merge(FLiP_lib, proteins[, c(1,2)])


correct_fraction <- c()
confidence <- c()
for(i in c(1:nrow(FLiP_lib))){
  
  MW <- FLiP_lib[i, "MW"]
  cutoff <- MW/2
  cutoff <- ifelse(cutoff > 100000, 100, ifelse(cutoff > 50000, 50, ifelse(cutoff > 30000, 30, ifelse(cutoff > 10000, 10, 0))))
  
  if(cutoff == 100){
    correct_fraction <- c(correct_fraction, "FALSE")
    confidence <- c(confidence, "low")
  
  }else if(cutoff == 50 & sum(FLiP_lib[i, c("Mean30K", "Mean10K")]) > 0){
    
    correct_fraction <- c(correct_fraction, "FALSE")
    
    p_changes <- FLiP_lib[i, c(15:20)]
    p_changes_wo_small_fracs <- as.vector(p_changes[-grep("30K|10K", colnames(p_changes))] )

    if(length(which(p_changes_wo_small_fracs < 0.05)) == 0) {
      confidence <- c(confidence, "low")
      }else{
      confidence <- c(confidence, "high")
      }
    
   
    }else if(cutoff == 30 & sum(FLiP_lib[i, c("Mean10K")]) > 0){
    correct_fraction <- c(correct_fraction, "FALSE")
    
    p_changes <- FLiP_lib[i, c(15:20)]
    p_changes_wo_small_fracs <- as.vector(p_changes[-grep("10K", colnames(p_changes))] )

    if(length(which(p_changes_wo_small_fracs < 0.05)) == 0) {
      confidence <- c(confidence, "low")
      }else{
        confidence <- c(confidence, "high")
      }

    }else{
    correct_fraction <- c(correct_fraction, "TRUE")
    confidence <- c(confidence, "high")
    }
}  

FLiP_lib$confidence <- confidence
FLiP_lib$correct_fraction <- correct_fraction


FLiP_lib <- FLiP_lib[, c("Protein", "MW", "Peptide", "fully_tryptic", "Position", "tryptic", "correct_fraction", "most_sig_change", "confidence")]

protein_names <- read.delim("../../Databases/210713_Yeast_Uniprot_Names_final.tsv", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
FLiP_lib <- merge(FLiP_lib, protein_names[, c(1,2)], by.x = "Protein", by.y = "Entry")
FLiP_lib <- FLiP_lib[, c("Protein", "Gene.names", "MW", "Peptide", "fully_tryptic", "Position", "tryptic", "correct_fraction", "most_sig_change", "confidence")]

write.csv(FLiP_lib, "./final_FLiP_lib/FiLiP_lib_marker_confidence.csv")


FLiP_library_grouped  <- FLiP_lib[, -c(4,6,7,8,9)] %>% group_by(fully_tryptic) %>% summarise(confidence = ifelse("high" %in% confidence, "high", "low"))

write.csv(FLiP_library_grouped, "./final_FLiP_lib/FiLiP_lib_marker_confidence_grouped.csv", row.names = FALSE, quote = FALSE)


table(FLiP_library_grouped$confidence)



