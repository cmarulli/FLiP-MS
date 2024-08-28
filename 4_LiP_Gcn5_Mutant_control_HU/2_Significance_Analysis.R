library("matrixStats")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("tidyr")
library("stringr")
library("tidyverse")
library("ggsci")
library("httr")
library("reticulate")
library("EnhancedVolcano")
library("wesanderson")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
source("Rfunctions/violin_plot.R")
source("Rfunctions/calc_mean_sd_pep_prot.R")
source("Rfunctions/anova_pep_prot.R")
source('Rfunctions/tukey_range_test.R')
source("Rfunctions/semi_to_fully.R")
source("Rfunctions/barplot_numbers_ggplot.R")
source("../1_FLiP_PBI_Library/Rfunctions/get_GO_annotation.R")
source("../1_FLiP_PBI_Library/Rfunctions/fisher_GO_test.R")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Color Palette
#-----------------------------------------------------------------------------------------------------------------------------------------------------
pal <- wes_palette("Darjeeling1", 40, type = "continuous")
pal_dis <- wes_palette("Darjeeling1", 5, type = "discrete")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------------------------------------------------------------------------------------
conditions <- c("Mutant", "MutantHU")
n <- length(conditions)
filename = "LiP_Gcn5_Mutant_control_HU"
rep = 4

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data
#-----------------------------------------------------------------------------------------------------------------------------------------------------
lip_data <- read.table(paste("./FilteredData/", filename, "_lip_data_filtered.tsv", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tc_data <- read.delim(paste("./FilteredData/", filename, "_tc_data_filtered.tsv", sep = ""), sep = ";", header = TRUE, stringsAsFactors = FALSE)
FLiP_marker <-  read.csv("../1_FLiP_PBI_Library/final_FLiP_lib/FiLiP_lib_marker_confidence.csv",  header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# Remove non-proteotypic peptides if present
if(length(which(lip_data$Proteotypic == "False")) > 0){
  lip_data <- lip_data[-which(lip_data$Proteotypic == "False"), ]
}

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the mean peptide/protein ratio and the corresponding standard error
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# subset the datasets to only contain relevant information
lip_sub <- lip_data[, c(3:(3+rep*n-1), 1, 2)]
tc_sub <- tc_data[, c(2:(3+rep*n-2), 1)]

# replace the proteins that have not been measured in a fraction by NA
tc_sub[tc_sub == 0] <- NA

list_pep_prot_mean_sd <- calc_mean_sd_pep_prot(lip_sub, tc_sub, n, column_names = conditions, filename = filename, rep = rep)
mean_peptides_proteins <- list_pep_prot_mean_sd[[1]]
sd_peptides_proteins <- list_pep_prot_mean_sd[[2]]
nb_replicates <- list_pep_prot_mean_sd[[3]]

#keep only peptides which have a mean peptide/protein ratio measured in at least 2 fractions
index <- apply(mean_peptides_proteins, 1, function(x) length(which(x != 0)))
mean_peptides_proteins <- mean_peptides_proteins[which(index >= 2), ]
sd_peptides_proteins <- sd_peptides_proteins[which(index >= 2), ]
nb_replicates <- nb_replicates[which(index >= 2), ]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# ANOVA test and Tukey's Range Test 
#-----------------------------------------------------------------------------------------------------------------------------------------------------
anova_summary <- anova_pep_prot(mean_peptides_proteins, sd_peptides_proteins, nb_replicates, lip_data, c("Mutant", "MutantHU"))
anova_summary <- merge(anova_summary, lip_data[,c('Peptide', 'Position', 'Protein')])

# add the information about the fully tryptic parent of a semi-tryptic peptide
anova_summary <- semi_to_fully(anova_summary)

# add Position stop information
anova_summary$Stop_Position <- apply(anova_summary, 1, function(x) return(as.numeric(x["Position"]) + nchar(as.character(x["Peptide"]))))
FLiP_marker$Stop_Position <- apply(FLiP_marker, 1, function(x) return(as.numeric(x["Position"]) + nchar(x["Peptide"])))

# add information about whether a peptide is a FLiP marker
has_FLiP <- c()
marker <- c()
confidence <- c()

for(i in c(1:nrow(anova_summary))){
  
  # if a peptide is present more than once in a protein sequence
  if(grepl(',', anova_summary[i, "Position"])){
    start <- strsplit(anova_summary[i,"Position"], ",")
    stop <- sapply(start, function(x) return(as.numeric(x)+nchar(anova_summary[i, "Peptide"])))
    df <- data.frame(start, stop)
    pep_pos <- as.vector((apply(df, 1, function(x) return(c(x[1]:x[2])))))
  }else{
    pep_pos <- c(anova_summary[i, "Position"]:anova_summary[i, "Stop_Position"])
  }
  
  markers <- FLiP_marker[which(FLiP_marker$Protein == anova_summary[i, "Protein"]), ]
  if(nrow(markers) > 1){
    marker_positions <- apply(markers, 1, function(x) return(c(x["Position"]:x["Stop_Position"])))
    # check whether at least half of the peptide is overlapping with a marker
    intersect <- unlist(lapply(marker_positions, function(x) return(length(intersect(pep_pos, x)) > length(pep_pos)/2)))
  }else if(nrow(markers == 1)){
    marker_positions <- c(markers[1,"Position"]:markers[1,"Stop_Position"])
    intersect <- length(intersect(pep_pos, marker_positions)) > length(pep_pos)/2
  }else{
    intersect <- FALSE
  }
  
  if(TRUE %in% intersect){
    
    has_FLiP <- c(has_FLiP, TRUE)
    
    # if there is a high confidence marker, report this one
    if("high" %in% markers[which(intersect == TRUE), "confidence"]){
      marker <- c(marker, markers[which(intersect == TRUE & markers$confidence == "high"), "Peptide"][1])
      confidence <- c(confidence, markers[which(intersect == TRUE & markers$confidence == "high"), "confidence"][1])
    }else{
      marker <- c(marker, markers[which(intersect == TRUE), "Peptide"][1])
      confidence <- c(confidence, markers[which(intersect == TRUE), "confidence"][1])
    }
    
  }else{
    has_FLiP <- c(has_FLiP, FALSE)
    marker <- c(marker, NA)
    confidence <- c(confidence, NA)
  }
}

anova_summary$has_FLiP <- has_FLiP
anova_summary$marker <- marker
anova_summary$confidence <- confidence

# add gene name information 
gene_names <- read.delim("../Files/210713_Yeast_Uniprot_Names_final.tsv", header = TRUE, stringsAsFactors = FALSE)
anova_summary$Gene <- sapply(anova_summary$Protein, function(x) return(gene_names[which(gene_names$Entry == x), "Gene.names"]))

# store information
write.table(anova_summary, paste("./OutputData/", filename, "_anova_results.csv", sep = ""),  row.names = FALSE, sep = ";")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# P-Body Analysis
#-----------------------------------------------------------------------------------------------------------------------------------------------------
p_body_name<- c("DCP2", "EDC3", "PAT1", "XRN1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", "LSM6", "LSM7", "UPF1", "DHH1")

anova_pbodies <- anova_summary[which(anova_summary$Gene %in% p_body_name), ]
nb_p_body_detected <- length(unique(anova_pbodies$Protein))
anova_pbodies_changed_FLiP <- anova_pbodies[which(anova_pbodies$q < 0.05 & anova_pbodies$has_FLiP == TRUE), ]
nb_p_body_changing <- length(unique(anova_pbodies_changed_FLiP$Gene))

anova_summary$type <- rep("all", nrow(anova_summary))
anova_summary[which(anova_summary$has_FLiP == TRUE), "type"] <- "FLiP"
anova_summary[which(anova_summary$Gene %in% p_body_name), "type"] <- "P-Body"
anova_summary[which(anova_summary$Gene %in% p_body_name & anova_summary$has_FLiP == TRUE), "type"] <- "P-Body FLiP"
anova_summary<- anova_summary[order(anova_summary$type), ] 


# only for plotting get rid of the p body proteins that are no FLiP markers 
anova_summary_plot <- anova_summary[-which(anova_summary$type == "P-Body"), ]

keyvals <- ifelse(anova_summary_plot$type == "P-Body FLiP", '#0057b7',
                  'grey')

names(keyvals)[keyvals == 'grey'] <- 'all'
names(keyvals)[keyvals == '#0057b7'] <- 'P-body'

labels <- anova_summary_plot$Gene
labels[which(keyvals == "grey")] <- ""

EnhancedVolcano(anova_summary_plot, x = "L2FC", y = "q",  pCutoff = 0.05, FCcutoff = 0,   colCustom = keyvals,
                lab = anova_summary_plot$Gene, 
                selectLab = anova_summary_plot[which(anova_summary_plot$type == "P-Body FLiP" & anova_summary_plot$q < 0.05), "Gene"],
                ylim = c(0,2.7), xlim = c(-5, 5), 
                title = "", subtitle = "", 
                pointSize = c(ifelse(anova_summary_plot$type %in% c("P-Body FLiP"), 2, 1)),
                colAlpha = c(ifelse(anova_summary_plot$type %in% c("P-Body FLiP"), 1, 0.3)),
                labSize = 2,
                labCol = 'black',
                gridlines.major = FALSE, gridlines.minor = FALSE, ylab = bquote(~-log[10]~italic(q)))
ggsave(paste("./Plots/", filename, "_volcano_anova_pbody_proteins.pdf", sep =""),  device = "pdf", width = 14 , height = 16, units = "cm")


# Volcano Plot without P-Body proteins highlighted
keyvlas <- ifelse(
  anova_summary$tryptic == "Specific", '#0057b7',
  ifelse(anova_summary$tryptic != "Specific", '#b0c7e0',
         'black'))
keyvlas[is.na(keyvlas)] <- 'black'

names(keyvlas)[keyvlas == '#0057b7'] <- 'fully'
names(keyvlas)[keyvlas == '#b0c7e0'] <- 'semi'

EnhancedVolcano(anova_summary, x = "L2FC", y = "q", lab = "", ylab = bquote(~-log[10]~italic(q)), pCutoff = 0.05, FCcutoff = 0.5,
                colCustom = keyvlas,
                ylim = c(0,2.6), xlim = c(-6, 6), 
                subtitle = "", pointSize = 0.5, gridlines.major = FALSE, gridlines.minor = FALSE, title = "", legendPosition = "")
ggsave(paste("./Plots/", filename, "_volcano_anova_semi_fully_tryptic.pdf", sep =""),  device = "pdf", width = 14 , height = 16, units = "cm")

anova_significiant_sign <- anova_summary[which(anova_summary$q < 0.05), ]
nb_sign_peptides <- length(unique(anova_significiant_sign$Peptide))
nb_sign_proteins <- length(unique(anova_significiant_sign$Protein))
nb_sign_peptides_grouped <- length(unique(anova_significiant_sign$fully_tryptic))
barplot_numbers_ggplot(c("Peptides", "Positions", "Proteins"), c(nb_sign_peptides, nb_sign_peptides_grouped, nb_sign_proteins), file_location = paste("./Plots/", filename, "NBSignificant", sep =""), y_lim = nb_sign_peptides + 500)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Protein Abundance Changes
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# replace 0 with NA
tc_data[tc_data == 0] <- NA

# keep only proteins that have been measured in both conditions
index <- apply(tc_data[, c(2:9)], 1, function(x) return(length(which(is.na(x)))))
tc_data <- tc_data[which(index < 4), ]

p_values <- c()
L2FCs <- c()

for(i in c(1:nrow(tc_data))){
  control_means <- mean(as.numeric(tc_data[i, c(2:5)]), na.rm = TRUE)
  HU_means <- mean(as.numeric(tc_data[i, c(6:9)]), na.rm = TRUE)     
  
  L2FC <- log2(HU_means/control_means)
  p <- t.test(as.numeric(tc_data[i, c(6:9)]), as.numeric(tc_data[i, c(2:5)]), na.rm = TRUE)$p.value
  
  p_values <- c(p_values, p)
  L2FCs <- c(L2FCs, L2FC)
}

q <- p.adjust(p_values, method = "BH")
summary <- data.frame(Protein = tc_data$Protein, L2FC = L2FCs, p = p_values, q = q)
write.csv(summary, paste("./OutputData/", filename, "_Protein_Abundance_Changes.csv", sep = ""), row.names = FALSE, quote = FALSE)

EnhancedVolcano(summary, x = "L2FC", y = "q", lab = "", pCutoff = 0.1, FCcutoff = 0.5, ylim = c(0,2.5), xlim = c(-3, 3), 
                col = c("gray40", "gray40","gray40","#0057b7"), ylab = bquote(~-log[10]~italic(q)),
                title = "", subtitle = "", pointSize = 1, gridlines.major = FALSE, gridlines.minor = FALSE, legendPosition = "") 
ggsave(paste("./Plots/", filename, "_Protein_Abundance_Change.pdf", sep = ''), device = "pdf", width = 14 , height = 16, units = "cm")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# GO Annotation from InterPro
#-----------------------------------------------------------------------------------------------------------------------------------------------------
GO_annotation <- get_GO_annotation(unique(anova_summary$Protein))
colnames(GO_annotation) <- c("Protein", "IPRs", "GO_names", "GO_identifiers", "Start_GO", "End_GO", "GO_type")
anova_summary <- merge(anova_summary, GO_annotation, by = 'Protein', all = TRUE)

# Reduce the dataset to GO-annotated peptides
anova_summary <- anova_summary[-which(is.na(anova_summary$GO_type)), ]

# get rid of peptides with a non uniquely assigned position
anova_summary <- anova_summary[-grep(",", anova_summary$Position), ]

# Keep only the GO annotations for which the location in the protein overlaps with the detected peptide 
anova_summary$index <- c(1:nrow(anova_summary))
peptides <- unique(anova_summary$Peptide)
index <- c()
for(i in c(1:length(peptides))){
  
  current_peptide <- anova_summary[which(anova_summary$Peptide == peptides[i]), c("Position", "Stop_Position", "Start_GO", "End_GO", "index")]
  intersect <- as.vector(apply(current_peptide, 1, function(x) length(intersect(c(x[1]:x[2]), c(x[3]:x[4])))))
  
  if(sum(intersect) > 1){
    index <- c(index, current_peptide[which(intersect!=0), "index"])
  }
}

anova_summary <- anova_summary[which(anova_summary$index %in% index), ]
write.table(anova_summary, "./OutputData/LiP_Mutant_control_HU_anova_results_GO.csv", quote = FALSE, row.names = FALSE, sep = ";")


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# GO Enrichment Analysis
# note that these plots can slightly differ from the ones published in the manuscript because the GO annotation is downloaded on a different data
# GO enrichemnt of significant peptides
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# add a binary column for significance
anova_summary_GO <- anova_summary
anova_summary_GO$sign <- rep(FALSE, nrow(anova_summary_GO))
anova_summary_GO[which(anova_summary_GO$q < 0.05), "sign"] <- TRUE

p_values_GO_terms <- fisher_GO_test(anova_summary_GO)
p_values_GO_terms <- p_values_GO_terms[which(p_values_GO_terms$p_values < 0.05), ]
p_values_GO_terms <- p_values_GO_terms[order(p_values_GO_terms$p_values, decreasing = TRUE), ]

p_values_GO_terms$count <- sapply(p_values_GO_terms$GO_Terms, function(x) length(which(anova_summary_GO$GO_names == x & anova_summary_GO$sign == TRUE)))
p_values_GO_terms <- p_values_GO_terms[order(p_values_GO_terms$count, decreasing = FALSE), ]
p_values_GO_terms$GO_Terms <- factor(p_values_GO_terms$GO_Terms, levels = p_values_GO_terms$GO_Terms)

ggplot(p_values_GO_terms, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave("./Plots/Mutant_HU_control_LiP_Interpro_GO_enrichment_all_sign_peptides.jpg",  device = "jpg", width = 25 , height = 18, units = "cm", )

# domain based only
anova_summary_GO_domain <- anova_summary_GO[which(anova_summary_GO$GO_type %in% c("domain", "repeat")), ]
p_values_GO_terms_domain <- fisher_GO_test(anova_summary_GO_domain)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(anova_summary_GO_domain$GO_names == x & anova_summary_GO_domain$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave("./Plots/Mutant_HU_control_LiP_Interpro_GO_enrichment_domain_sign_peptides.jpg",  device = "jpg", width = 20 , height = 16, units = "cm")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# GO Enrichment Analysis of shortlist FLiP Hits
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# add a binary column for significance
anova_summary_GO <- anova_summary
anova_summary_GO$sign <- rep(FALSE, nrow(anova_summary_GO))
anova_summary_GO[which(anova_summary_GO$q < 0.05 & anova_summary_GO$has_FLiP == TRUE), "sign"] <- TRUE

p_values_GO_terms <- fisher_GO_test(anova_summary_GO)
p_values_GO_terms <- p_values_GO_terms[which(p_values_GO_terms$p_values < 0.05), ]
p_values_GO_terms <- p_values_GO_terms[order(p_values_GO_terms$p_values, decreasing = TRUE), ]

p_values_GO_terms$count <- sapply(p_values_GO_terms$GO_Terms, function(x) length(which(anova_summary_GO$GO_names == x & anova_summary_GO$sign == TRUE)))
p_values_GO_terms <- p_values_GO_terms[order(p_values_GO_terms$count, decreasing = FALSE), ]
p_values_GO_terms$GO_Terms <- factor(p_values_GO_terms$GO_Terms, levels = p_values_GO_terms$GO_Terms)

ggplot(p_values_GO_terms, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave("./Plots/Mutant_HU_control_LiP_Interpro_GO_enrichment_all_sign_FLiP_peptides.jpg",  device = "jpg", width = 25 , height = 18, units = "cm", )

# domain based only
anova_summary_GO_domain <- anova_summary_GO[which(anova_summary_GO$GO_type %in% c("domain", "repeat")), ]
p_values_GO_terms_domain <- fisher_GO_test(anova_summary_GO_domain)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(anova_summary_GO_domain$GO_names == x & anova_summary_GO_domain$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)
 
ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave("./Plots/Mutant_HU_control_LiP_Interpro_GO_enrichment_domain_sign_FLiP_peptides.jpg",  device = "jpg", width = 20 , height = 16, units = "cm")

nrow(anova_summary_GO_domain[which(anova_summary_GO_domain$q < 0.05 & anova_summary_GO_domain$GO_names == "protein binding"), ])
nrow(anova_summary_GO_domain[which(anova_summary_GO_domain$q < 0.05 & anova_summary_GO_domain$has_FLiP == TRUE & anova_summary_GO_domain$GO_names == "protein binding"), ])

