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
#parameters
#-----------------------------------------------------------------------------------------------------------------------------------------------------
conditions <- c("WT_Control", "WT_HU")
n <- length(conditions)
filename <- "Acetylation_WT_control_HU"

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Color Palette
#-----------------------------------------------------------------------------------------------------------------------------------------------------
pal <- wes_palette("Darjeeling1", 40, type = "continuous")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Functions loaded from files in Rfunctions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
source("Rfunctions/semi_fully_tryptic.R")
source("Rfunctions/violin_plot.R")
source("Rfunctions/triplicates_lip_CD_CD.R")
source("Rfunctions/triplicates_tc_CD.R")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in Data and rename columns
#-----------------------------------------------------------------------------------------------------------------------------------------------------

# WT Control vs HU
acetyl_data <- read.delim("./Spectronaut/230504_CD010_actylation_WT_control_HU_Report_merged_library_peptides_imputed.xls", header = TRUE, stringsAsFactors = FALSE)
tc_data <- read.delim("./Spectronaut/230504_CD010_actylation_WT_control_HU_Report_merged_library_proteins_imputed.xls", header = TRUE, stringsAsFactors = FALSE)

# WT vs Mutant Control
#acetyl_data <- read.delim("./Spectronaut/20230516_104849_230505_CD010_CD014_acetylation_WT_mutant_control_CD010_CD011_CD014_merged_library_Report_peptides.xls", header = TRUE, stringsAsFactors = FALSE)
#tc_data <- read.delim("./Spectronaut/20230516_104941_230505_CD010_CD014_acetylation_WT_mutant_control_CD010_CD011_CD014_merged_library_Report_proteins.xls", header = TRUE, stringsAsFactors = FALSE)

# WT vs Mutant HU
#acetyl_data <- read.delim("./Spectronaut/20230516_104548_230505_CD010_CD014_acetylation_WT_mutant_HU_CD010_CD011_CD014_library_Report_peptides.xls", header = TRUE, stringsAsFactors = FALSE)
#tc_data <- read.delim("./Spectronaut/20230516_104734_230505_CD010_CD014_acetylation_WT_mutant_HU_CD010_CD011_CD014_library_Report_proteins.xls", header = TRUE, stringsAsFactors = FALSE)

# Mutant Control vs HU
#acetyl_data <- read.delim("./Spectronaut/20230516_114912_230516_CD014_aceylation_mutant_HU_control_CD010_CD011_CD014_merged_lib_Report_peptides.xls", header = TRUE, stringsAsFactors = FALSE)
#tc_data <- read.delim("./Spectronaut/20230516_114942_230516_CD014_aceylation_mutant_HU_control_CD010_CD011_CD014_merged_lib_Report_proteins.xls", header = TRUE, stringsAsFactors = FALSE)

# remove duplicated peptide rows
acetyl_data <- acetyl_data[-which(duplicated(acetyl_data)), ]

# doerichr 10
colnames(acetyl_data) <- c("Protein", "Peptide", "Proteotypic", "Position", "mod_Peptide", 
                        "WT_Control_4", "WT_Control_3",
                        "WT_HU_4",
                        "WT_Control_1", "WT_Control_2",
                        "WT_HU_1", "WT_HU_3", "WT_HU_2")
                        #"Mutant_HU_2",
                        #"Mutant_Control_2", "Mutant_Control_4", 
                        #"Mutant_HU_3", "Mutant_HU_1", "Mutant_HU_4", 
                        #"Mutant_Control_1", "Mutant_Control_3")

colnames(tc_data) <- c("Protein",
                       "WT_Control_4", "WT_Control_3",
                       "WT_HU_4",
                       "WT_Control_1", "WT_Control_2",
                       "WT_HU_1", "WT_HU_3", "WT_HU_2")
                        #"Mutant_HU_2",
                        #"Mutant_Control_2", "Mutant_Control_4", 
                        #"Mutant_HU_3", "Mutant_HU_1", "Mutant_HU_4", 
                        #"Mutant_Control_1", "Mutant_Control_3")

acetyl_data <- acetyl_data[, c("mod_Peptide", "Protein", 
                               "WT_Control_1", "WT_Control_2", "WT_Control_3", "WT_Control_4",
                               "WT_HU_1", "WT_HU_2", "WT_HU_3", "WT_HU_4",
                               #"Mutant_Control_1", "Mutant_Control_2", "Mutant_Control_3", "Mutant_Control_4",
                               #"Mutant_HU_1", "Mutant_HU_2", "Mutant_HU_3", "Mutant_HU_4", 
                               "Peptide", "Position")]

tc_data <- tc_data[, c( "Protein", 
                        "WT_Control_1", "WT_Control_2", "WT_Control_3", "WT_Control_4", 
                        "WT_HU_1", "WT_HU_2", "WT_HU_3", "WT_HU_4")]
                        #"Mutant_Control_1", "Mutant_Control_2", "Mutant_Control_3", "Mutant_Control_4", 
                        #"Mutant_HU_1", "Mutant_HU_2", "Mutant_HU_3", "Mutant_HU_4")]
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Keep only acetylated peptides
#-----------------------------------------------------------------------------------------------------------------------------------------------------
acetyl_data <-  acetyl_data[grepl( "Acetyl (K)", acetyl_data$mod_Peptide, fixed = TRUE),]
acetyl_data <- unique(acetyl_data)

print(paste("Number of acetylated peptides", nrow(acetyl_data)))
print(paste("Coming from", length(unique(acetyl_data$Protein)), "Proteins"))

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Keep only LiP peptides measured in at least triplicates
#-----------------------------------------------------------------------------------------------------------------------------------------------------
lip_result_list <- triplicates_lip(acetyl_data, n, conditions, filename)
acetyl_data_temp <- lip_result_list[[1]]
acetyl_data <- merge(acetyl_data_temp, acetyl_data[, c("mod_Peptide", "Position")])
nb_identified_peptides <- lip_result_list[[2]]
identified_peptides <- lip_result_list[[3]]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# LiP-Data Summary: Plot of Peptides and their nature in distinct fractions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# plot the number of pepetides in each fraction
# need a folder called "Plots" where all the plots are being stored
df_nb_peptides <- data.frame(fraction = conditions, number = nb_identified_peptides)
df_nb_peptides$fraction <- factor(df_nb_peptides$fraction, levels = df_nb_peptides$fraction)

ggplot(df_nb_peptides, aes(x = fraction, y = number, fill = c(1:n))) + geom_bar( stat="identity", show.legend = FALSE) +
  scale_fill_gradientn(colours = pal) + theme_classic(base_size = 26) + xlab("") + ylab ("Number of Acetylated Peptides") 

ggsave(paste("./Plots/", filename, "_Number_Peptides.jpg", sep =""), device = "jpg", width = 16 , height = 12, units = "cm")

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
jpeg(paste("./Plots/", filename, "Protein_Heatmap.jpg", sep =""),  width = 24, height = 14,units = "cm", res = 1000)
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
index <- which(acetyl_data$Protein %in% tc_data$Protein)
acetyl_data <- acetyl_data[index, ]
acetyl_data <- unique(acetyl_data)

# export the normalized and filtered acetyl_data
write.table(acetyl_data, paste("./FilteredData/", filename, "_acetyl_data_filtered.tsv", sep = ""),  row.names = FALSE, quote = FALSE, sep = "\t")