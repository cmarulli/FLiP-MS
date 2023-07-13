library("ggplot2")
library("ggpubr")
library("reshape2")
library("tidyr")
library("stringr")
library("tidyverse")
library("ggsci")
library("httr")
library("reticulate")
library("dplyr")
library("UpSetR")
library("wesanderson")
#-----------------------------------------------------------------------------------------------------------------------------------------------------
#paramters
conditions <- c("Mutant", "MutantHU")
n <- length(conditions)
filename = "LiP_Gcn5_Mutant_control_HU"

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Color Palette
#-----------------------------------------------------------------------------------------------------------------------------------------------------
pal <- wes_palette("Darjeeling1", 40, type = "continuous")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Functions loaded from files in Rfunctions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
source("Rfunctions/semi_fully_tryptic.R")
source("Rfunctions/violin_plot.R")
source("Rfunctions/triplicates_lip.R")
source("Rfunctions/triplicates_tc.R")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in Data and rename columns
#-----------------------------------------------------------------------------------------------------------------------------------------------------
lip_data <- read.delim("./Spectronaut/20230303_135638_CD014_SAGA_Mutant_LiP_Report_mutant_only.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
tc_data <- read.delim("./Spectronaut/20230303_135638_CD014_SAGA_Mutant_TC_Report_mutant_only.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# original annotation
colnames(lip_data) <- c("Protein", "Peptide", "Proteotypic", "Position", "tryptic",
                        "Mutant_2", "Mutant_4", "MutantHU_1", "MutantHU_4",
                        "Mutant_1", "Mutant_3", "MutantHU_2", "MutantHU_3")

colnames(tc_data) <- c("Protein",
                       "Mutant_2", "Mutant_4", "MutantHU_1", "MutantHU_4",
                       "Mutant_1", "Mutant_3", "MutantHU_2", "MutantHU_3")
                       
lip_data <- lip_data[, c("Peptide", "Protein",
                         "Mutant_1", "Mutant_2", "Mutant_3", "Mutant_4", 
                         "MutantHU_1",  "MutantHU_2",  "MutantHU_3",  "MutantHU_4", 
                         "tryptic", "Position", "Proteotypic")]

tc_data <- tc_data[, c( "Protein",
                        "Mutant_1", "Mutant_2", "Mutant_3", "Mutant_4", 
                        "MutantHU_1",  "MutantHU_2",  "MutantHU_3",  "MutantHU_4" 
                        )]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Keep only LiP peptides measured in at least triplicates
#-----------------------------------------------------------------------------------------------------------------------------------------------------
lip_result_list <- triplicates_lip(lip_data, n, conditions, filename)
lip_data <- lip_result_list[[1]]
nb_identified_peptides <- lip_result_list[[2]]
identified_peptides <- lip_result_list[[3]]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# LiP-Data Summary: Plot of peptides and their nature in distinct fractions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# plot the number of peptides in each fraction
# need a folder called "Plots" where all the plots are being stored
df_nb_peptides <- data.frame(fraction = conditions, number = nb_identified_peptides)
df_nb_peptides$fraction <- factor(df_nb_peptides$fraction, levels = df_nb_peptides$fraction)

ggplot(df_nb_peptides, aes(x = fraction, y = number, fill = c(1:n))) + geom_bar( stat="identity", show.legend = FALSE) +
  scale_fill_gradientn(colours = pal) + theme_classic(base_size = 26) + xlab("") + ylab ("Number of Peptides") 

ggsave(paste("./Plots/", filename, "_Number_Peptides.jpg", sep =""), device = "jpg", width = 16 , height = 12, units = "cm")

# plot the percentage of summed intensity of semi and fully tryptic peptides with a barplot
lip_list <- list()
tryptic_column <- which(colnames(lip_data) == "tryptic")
for(i in c(1:n)){
  lip_list[[i]] <- lip_data[, c(1,2, ((i-1)*4+3):((i-1)*4+6), tryptic_column)]
}

intesisty_semi_trytic_peptides <- unlist(lapply(lip_list, function(x) sum(colSums(x[which(x$tryptic != "Specific"), c(3:6)]))))
intensity_fully_trytic_peptides <- unlist(lapply(lip_list, function(x) sum(colSums(x[which(x$tryptic == "Specific"), c(3:6)]))))

percentage_semi_tryptic_peptide_intensities <- as.data.frame(intesisty_semi_trytic_peptides /(intesisty_semi_trytic_peptides+intensity_fully_trytic_peptides))
colnames(percentage_semi_tryptic_peptide_intensities) <- "percentage"
percentage_semi_tryptic_peptide_intensities$fraction <- conditions
percentage_semi_tryptic_peptide_intensities$fraction <- factor(percentage_semi_tryptic_peptide_intensities$fraction, levels = unique(percentage_semi_tryptic_peptide_intensities$fraction))

ggplot(percentage_semi_tryptic_peptide_intensities, aes(x = fraction, y = percentage, fill = c(1:n))) + geom_bar( stat="identity", show.legend = FALSE) +
  scale_fill_gradientn(colours = pal) + theme_classic(base_size = 18) + xlab("") + ylab ("Percentage") + ggtitle("Intensities of Semi-Tryptic Peptides") 

ggsave(paste("./Plots/", filename, "_Percentage_Semi_typtic_LiP.jpg", sep =""), device = "jpg", width = 16 , height = 12, units = "cm")

# heatmap of proteins in each fraction (! even though tc filtered data is exactly the same as for Master Thesis Analysis the Heatmap looks different)
jpeg(paste("./Plots/", filename, "_Peptide_Heatmap.jpg", sep =""),  width = 24, height = 14,units = "cm", res = 1000)
heatmap(as.matrix(lip_data[, c(3:10) ]), col =  heat.colors(25))
dev.off()

# UpSet Plot gives intersection of measured proteins between all conditions
jpeg(paste("./Plots/UpSetR_",filename,  "_Peptides", ".jpg", sep = ""), width = 16, height = 12,units = "cm", res = 1000)
  upset(fromList(identified_peptides), order.by = "freq")
dev.off()


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# TC-Data Filtering: keep only proteins that have been measured in at least triplicates
#-----------------------------------------------------------------------------------------------------------------------------------------------------
tc_result_list <- triplicates_tc(tc_data, n, conditions, filename)
tc_data <- tc_result_list[[1]]
nb_identified_proteins <- tc_result_list[[2]]
identified_proteins <- tc_result_list[[3]]
if(length(which(duplicated(tc_data$Protein))) > 0){
  tc_data <- tc_data[-which(duplicated(tc_data$Protein)), ]
}
write.table(tc_data, paste("./FilteredData/", filename, "_tc_data_filtered.tsv", sep = ""),  row.names = FALSE, quote = FALSE, sep = ";")


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# TC Data Summary Plots
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# number of proteins per fraction
df_nb_proteins <- data.frame(fraction = conditions , number = nb_identified_proteins)
df_nb_proteins$fraction <- factor(df_nb_proteins$fraction, levels = df_nb_proteins$fraction)

# number of proteins that has been identified in each fraction
ggplot(df_nb_proteins, aes(x = fraction, y = number, fill = c(1:n))) + geom_bar( stat="identity", show.legend = FALSE) +
  scale_fill_gradientn(colours = pal) + theme_classic(base_size = 26) + xlab("") + ylab ("Number of Proteins") 

ggsave(paste("./Plots/", filename, "_Number_Proteins.jpg", sep = ""), device = "jpg", width = 16 , height = 12, units = "cm")

# heatmap of proteins in each fraction (! even though tc filtered data is exactly the same as for Master Thesis Analysis the Heatmap looks different)
jpeg(paste("./Plots/", filename, "_Protein_Heatmap.jpg", sep =""),  width = 24, height = 14,units = "cm", res = 1000)
  heatmap(as.matrix(tc_data[, c(2: (ncol(tc_data))) ]), col =  heat.colors(25))
dev.off()

# UpSet Plot gives intersection of measured proteins between all conditions
jpeg(paste("./Plots/UpSetR_",filename, "_Proteins.jpg", sep = ""), width = 16, height = 12,units = "cm", res = 1000)
  upset(fromList(identified_proteins), order.by = "freq")
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Subset LiP Data to only contain the Proteins present in the TC data and safe
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# merge lip data with the molecular weight of the tryptic control data - like this we can get rid of peptides from proteins which have not been quantified in the tryptic control
index <- which(lip_data$Protein %in% tc_data$Protein)
lip_data <- lip_data[index, ]
lip_data <- unique(lip_data)

# export the normalized and filtered lip_data
write.table(lip_data, paste("./FilteredData/", filename, "_lip_data_filtered.tsv", sep = ""),  row.names = FALSE, quote = FALSE, sep = "\t")