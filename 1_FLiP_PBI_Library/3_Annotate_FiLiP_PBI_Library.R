library(dplyr)
library(readxl)
library(ggplot2)
library(protti)
library(reshape2)

source("../../Rfunctions/barplot_numbers_ggplot.R")
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# read in FLiP library
FLiP_lib <- read.csv("./OutputData/FLiP_PBI_Library_anova_results.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
# read in the detected proteins
proteins <- read.delim("./FilteredData/FLiP_PBI_library_tc_data_filtered.tsv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
proteins[proteins == 0] <- NA

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Sequence Coverage per Abundance Bin
#-----------------------------------------------------------------------------------------------------------------------------------------------------
FLiP_lib_protein_coverage <- FLiP_lib[, c("Peptide", "Protein", "q", "Sequence")]
FLiP_lib_protein_coverage <- calculate_sequence_coverage(FLiP_lib_protein_coverage, Sequence, Peptide)

proteins$mean_intensity <- log10(rowMeans(proteins[, c(3:18)], na.rm = TRUE))
proteins <- unique(merge(proteins, FLiP_lib_protein_coverage[, c("Protein", "coverage")]))

bin_size <- 0.5

proteins <- proteins %>% mutate(bin_dist = factor(mean_intensity%/%bin_size)) 
#proteins <- proteins[-which(proteins$bin_dist %in% c(7, 15)), ]

ggplot(proteins, aes(x = bin_dist, y = coverage)) + theme_classic(base_size = 20) + geom_boxplot() +
  scale_x_discrete(labels = paste0(sort(as.numeric(unique(levels(proteins$bin_dist)))*0.5), "-", sort((as.numeric(unique(levels(proteins$bin_dist))) *0.5 + 0.5) ))) +
  labs(x = "log10(Protein Intensities)", y = "Sequenc Coverage")
ggsave("./Plots/Sequence_Coverage_FLiP_complete_dataset.pdf", device = "pdf", width = 24 , height = 14, units = "cm", useDingbats=FALSE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Bar Diagram Abundance classes in filter fractions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# calculate the protein abundance in the respective fractions
proteins$mean_100K <-  log10(rowSums(proteins[, c(3:6)], na.rm = TRUE)/4)
proteins$mean_50K <-  log10(rowSums(proteins[, c(7:10)], na.rm = TRUE)/4)
proteins$mean_30K <-  log10(rowSums(proteins[, c(11:14)], na.rm = TRUE)/4)
proteins$mean_10K <-  log10(rowSums(proteins[, c(15:18)], na.rm = TRUE)/4)
proteins[sapply(proteins, is.infinite)] <- NA

# check if the protein has a FLiP marker 
proteins$is_FLiP <- rep(FALSE, nrow(proteins)) 
FLiP_proteins <- unique(FLiP_lib[which(FLiP_lib$q < 0.05), "Protein"])
proteins[which(proteins$Protein %in% FLiP_proteins), "is_FLiP"] <- TRUE

proteins <- proteins[, c(1,19,22:26)]

# split proteins into 5 equally sized bins
proteins <- proteins[order(proteins$mean_intensity, decreasing = TRUE), ]
proteins$bin <- c(rep(1, 274), rep(2, 274), rep(3,274), rep(4,274), rep(5,275))

df_detected <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df_detected) <- c(1,2,3,4,5)

for(c in c('100K', '50K', '30K', '10K')){
  
  current_fraction <- proteins[ , c(grep(c, colnames(proteins)), 8)]
  current_fraction <- current_fraction[-which(is.na(current_fraction[, 1])), ]
  bins <- as.vector(table(current_fraction$bin))
  df_detected <- rbind(df_detected, bins)
  
}

df_detected$fraction <- c('100K', '50K', '30K', '10K')
colnames(df_detected) <- c("very high","high","medium","low","very low", "fraction")

df_detected <- melt(df_detected, id.vars = 'fraction', variable.name = "Abundance")

ggplot(df_detected, aes(fill=Abundance, y=value, x=fraction)) + 
  geom_bar(position="stack", stat="identity")

# repeat the same analysis for proteins with FLiP peptides
df_FLiP <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df_FLiP) <- c(1,2,3,4,5)

proteins_FLiP <- proteins[which(proteins$is_FLiP == TRUE), ]

for(c in c('100K', '50K', '30K', '10K')){
  
  current_fraction <- proteins_FLiP[ , c(grep(c, colnames(proteins_FLiP)), 8)]
  current_fraction <- current_fraction[-which(is.na(current_fraction[, 1])), ]
  bins <- as.vector(table(current_fraction$bin))
  df_FLiP <- rbind(df_FLiP, bins)
  
}

df_FLiP$fraction <- c('100K_FLiP', '50K_FLiP', '30K_FLiP', '10K_FLiP')
colnames(df_FLiP) <- c("very high","high","medium","low","very low", "fraction")

df_FLiP <- melt(df_FLiP, id.vars = 'fraction', variable.name = "Abundance")

ggplot(df_FLiP, aes(fill=Abundance, y=value, x=fraction)) + 
  geom_bar(position="stack", stat="identity")

df_all <- rbind(df_detected, df_FLiP)
df_all$fraction <- factor(df_all$fraction, levels = c("100K", "50K", "30K","10K","100K_FLiP","50K_FLiP","30K_FLiP","10K_FLiP"))
#c("100K", "100K_FLiP", "50K", "50K_FLiP", "30K", "30K_FLiP", "10K", "10K_FLiP")
colfunc <- colorRampPalette(c("#001e64","#0057b7","#b0c7e0"))

ggplot(df_all, aes(fill=Abundance, y=value, x=fraction, label = value)) + 
  geom_bar(position="stack", stat="identity") + theme_classic(base_size = 15) + scale_fill_manual(values = colfunc(5)) + scale_color_manual(values = colfunc(5)) +
  xlab("Fraction") + ylab("Number of Proteins in Abundance Bin") + geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white")
ggsave("./Plots/BarPlots_AbundanceClass_Detection_FLiP.pdf", device = "pdf", width = 24 , height = 14, units = "cm", useDingbats=FALSE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Filter for significant peptides only
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# read in FLiP library
FLiP_lib <- read.csv("./OutputData/FLiP_PBI_Library_anova_results.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
# read in the detected proteins
proteins <- read.delim("./FilteredData/FLiP_PBI_library_tc_data_filtered.tsv", header = TRUE, stringsAsFactors = FALSE, sep = ";")


FLiP_lib <- FLiP_lib[which(FLiP_lib$q < 0.05), ]
FLiP_lib <- merge(FLiP_lib, proteins[, c(1,2)])

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# If a protein found in a fraction where the cutoff is below half of it's molecular weight, we consider it to be in the wrong fraction
# If only the protein in the wrong fraction explains the signal in the FLiP library, we consider the marker low confidence
#-----------------------------------------------------------------------------------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# add protein names
#-----------------------------------------------------------------------------------------------------------------------------------------------------
protein_names <- read.delim("../Files/210713_Yeast_Uniprot_Names_final.tsv", header = TRUE, stringsAsFactors = FALSE, sep = '\t')
FLiP_lib <- merge(FLiP_lib, protein_names[, c(1,2)], by.x = "Protein", by.y = "Entry")
FLiP_lib <- FLiP_lib[, c("Protein", "Gene.names", "MW", "Peptide", "fully_tryptic", "Position", "tryptic", "correct_fraction", "most_sig_change", "confidence")]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# add InterPro domain annotation about being protein binding, RNA binding or DNA binding
#-----------------------------------------------------------------------------------------------------------------------------------------------------
InterPro_GO <- read.csv("./OutputData/FLiP_PBI_Library_lip_data_GO_domain_repeat.csv", header = TRUE, stringsAsFactors = FALSE, sep = ';')
InterPro_GO <- InterPro_GO[which(InterPro_GO$GO_names %in% c("protein binding", "RNA binding", "DNA binding")), c("Peptide", "GO_names")]
FLiP_lib <- merge(FLiP_lib, InterPro_GO, all = TRUE)
FLiP_lib <- unique(FLiP_lib[-which(is.na(FLiP_lib$Protein)), ])

table(FLiP_lib$GO_names)
length(unique(FLiP_lib[which(FLiP_lib$GO_names == "protein binding"), "Protein"]))
length(unique(FLiP_lib[which(FLiP_lib$GO_names == "RNA binding"), "Protein"]))

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# add PDB interface information
#-----------------------------------------------------------------------------------------------------------------------------------------------------
distances_best <- read.csv("./Distance_to_interface/Average_distance_best_2_6.csv", header = TRUE, stringsAsFactors = FALSE)
distances_best$at_interface <- sapply(distances_best$distances, function(x) return(ifelse(x <= 2.6, TRUE, FALSE)))
FLiP_lib <- merge(FLiP_lib, distances_best[, c("peptide", "at_interface")], by.x = "Peptide", by.y = "peptide", all = TRUE)
FLiP_lib <- FLiP_lib[-which(is.na(FLiP_lib$Protein)), ]

colnames(FLiP_lib)[12] <- "PDB_PBI"
table(FLiP_lib$PDB_PBI)
length(unique(FLiP_lib[which(FLiP_lib$PDB_PBI == TRUE), "Protein"]))

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# any evidence for protein binding
#-----------------------------------------------------------------------------------------------------------------------------------------------------
FLiP_lib$PBI_evidence <- apply(FLiP_lib, 1, function(x) return(ifelse( (x["GO_names"] == "protein binding" | x["PDB_PBI"] == "TRUE"), "TRUE", "FALSE")))
FLiP_lib[which(is.na(FLiP_lib$PBI_evidence)), "PBI_evidence"] <- FALSE

table(FLiP_lib$PBI_evidence)
length(unique(FLiP_lib[which(FLiP_lib$PBI_evidence == TRUE), "Protein"]))

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# To how many known complexes in the ComplexPortal database do those markers map
#-----------------------------------------------------------------------------------------------------------------------------------------------------
complex_portal <- read.delim("../Files/complex_portal_559292.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
FLiP_lib$ComplexPortal <- unlist(sapply(FLiP_lib$Protein, function(x) 
  return(ifelse(TRUE %in% grepl(x, complex_portal$Identifiers..and.stoichiometry..of.molecules.in.complex), 
                paste(complex_portal[grep(x, complex_portal$Identifiers..and.stoichiometry..of.molecules.in.complex), "Recommended.name"], collapse = "|"), NA))))

FLiP_lib_complex_portal <- FLiP_lib[which(!is.na(FLiP_lib$ComplexPortal)), ]
length(unique(FLiP_lib_complex_portal$Protein))
length(unique(unlist(strsplit(FLiP_lib_complex_portal$ComplexPortal, "\\|"))))
       
       
FLiP_lib_complex_portal_PBI <- FLiP_lib_complex_portal[which(FLiP_lib_complex_portal$PBI_evidence == TRUE), ]
length(unique(FLiP_lib_complex_portal_PBI$Protein))
length(unique(unlist(strsplit(FLiP_lib_complex_portal_PBI$ComplexPortal, "\\|"))))

write.csv(FLiP_lib, "./final_FLiP_lib/FiLiP_lib_marker_confidence.csv")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Group overlapping peptides
#-----------------------------------------------------------------------------------------------------------------------------------------------------
FLiP_library_grouped  <- unique(FLiP_lib[, -c(1,4,6,7,8,9, 11, 12)] %>% group_by(fully_tryptic) %>% mutate(confidence = ifelse("high" %in% confidence, "high", "low"), 
                                                                                             PBI_evidence = ifelse("TRUE" %in% PBI_evidence, "TRUE", "FALSE"), 
                                                                                             ComplexPortal = ifelse(is.na(ComplexPortal), NA, unique(ComplexPortal))))
table(FLiP_library_grouped$PBI_evidence)
length(unique(unlist(FLiP_library_grouped[which(FLiP_library_grouped$PBI_evidence == TRUE), "Protein"])))

write.csv(FLiP_library_grouped, "./final_FLiP_lib/FiLiP_lib_marker_confidence_grouped.csv", row.names = FALSE, quote = FALSE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# plot the fraction between which the most significant changes is detected
#-----------------------------------------------------------------------------------------------------------------------------------------------------
table(FLiP_lib$most_sig_change)
barplot_numbers_ggplot(conditions = names(table(FLiP_lib$most_sig_change)), number = table(FLiP_lib$most_sig_change), 
                       file_location = "./Plots/Fractions_most_sign_change", y_lim = 2200, my_width = 30)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Copy number distribution 
#-----------------------------------------------------------------------------------------------------------------------------------------------------
copy_numbers <- read_xlsx("../Files/yeast_copy_numbers.xlsx")
colnames(copy_numbers)[9] <- "Copy_number"
copy_numbers$Copy_number <- log10(copy_numbers$Copy_number)
copy_numbers$type <- rep("all", nrow(copy_numbers))

copy_numbers_FLiP <- copy_numbers[which(copy_numbers$`Majority protein IDs` %in% FLiP_lib$Protein), ]
copy_numbers_FLiP$type <- rep("FLiP Marker", nrow(copy_numbers_FLiP))

copy_numbers_all <- rbind(copy_numbers, copy_numbers_FLiP)

ggplot(copy_numbers_all, aes(x = Copy_number, fill = type, color = type)) + geom_histogram(binwidth = 0.1) + theme_classic() + 
  scale_color_manual(values=c("grey", "#0057b7")) + scale_fill_manual(values=c("grey", "#0057b7")) + geom_vline(xintercept = quantile(copy_numbers_FLiP$Copy_number, 0.1)) +
  xlab("Log10 Copy Numbers") + ylab("Number of Proteins")
ggsave("./Plots/Copy_numbers_distribution.pdf", device = "pdf", width = 24 , height = 14, units = "cm", useDingbats=FALSE)

q_FLiP_001 <- quantile(copy_numbers_FLiP$Copy_number, 0.10, na.rm = TRUE)

length(which(copy_numbers$Copy_number > q_FLiP_001)) / nrow(copy_numbers)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# To how many known complexes do those markers map
#-----------------------------------------------------------------------------------------------------------------------------------------------------
complex_portal <- read.delim("../Files/complex_portal_559292.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
FLiP_lib$ComplexPortal <- unlist(sapply(FLiP_lib$Protein, function(x) 
  return(ifelse(TRUE %in% grepl(x, complex_portal$Identifiers..and.stoichiometry..of.molecules.in.complex), 
                paste(complex_portal[grep(x, complex_portal$Identifiers..and.stoichiometry..of.molecules.in.complex), "Recommended.name"], collapse = "|"), NA))))

FLiP_lib_complex_portal <- FLiP_lib[which(!is.na(FLiP_lib$ComplexPortal)), ]
length(unique(FLiP_lib_complex_portal$Protein))
FLiP_lib_complex_portal_PBI <- FLiP_lib_complex_portal[which(FLiP_lib_complex_portal$PBI_evidence == TRUE), ]
length(unique(FLiP_lib_complex_portal_PBI$Protein))

