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
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Color Palette
#-----------------------------------------------------------------------------------------------------------------------------------------------------
pal <- wes_palette("Darjeeling1", 40, type = "continuous")
pal_dis <- wes_palette("Darjeeling1", 5, type = "discrete")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------------------------------------------------------------------------------------
conditions <- c("Control", "HU")
n <- length(conditions)
filename = "LiP_WT_control_HU"
rep = 4

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data
#-----------------------------------------------------------------------------------------------------------------------------------------------------
lip_data <- read.table(paste("./FilteredData/", filename, "_lip_data_filtered.tsv", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tc_data <- read.delim(paste("./FilteredData/", filename, "_tc_data_filtered.tsv", sep = ""), sep = ";", header = TRUE, stringsAsFactors = FALSE)
FLiP_marker <-  read.csv("../1_FLiP_PBI_Library/final_FLiP_lib/FiLiP_lib_marker_confidence.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

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
# ANOVA test 
#-----------------------------------------------------------------------------------------------------------------------------------------------------
anova_summary <- anova_pep_prot(mean_peptides_proteins, sd_peptides_proteins, nb_replicates, lip_data, c("Control", "HU"))
anova_summary <- merge(anova_summary, lip_data[,c('Peptide', 'Position', 'Protein')])

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Annotate Anova Results
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# add the information about the fully tryptic parent of a semi-tryptic peptide
anova_summary <- semi_to_fully(anova_summary)

# add information about stop position
anova_summary$Stop_Position <- apply(anova_summary, 1, function(x) return(as.numeric(x[12])+nchar(x[2])))
FLiP_marker$Stop_Position <- apply(FLiP_marker, 1, function(x) return(as.numeric(x[6]) + nchar(x[4])))

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
    marker_positions <- apply(markers, 1, function(x) return(c(x[6]:x[11])))
    # check whether at least half of the peptide is overlapping with a marker
    intersect <- unlist(lapply(marker_positions, function(x) return(length(intersect(pep_pos, x)) > length(pep_pos)/2)))
  }else if(nrow(markers == 1)){
    marker_positions <- c(markers[1,6]:markers[1,11])
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