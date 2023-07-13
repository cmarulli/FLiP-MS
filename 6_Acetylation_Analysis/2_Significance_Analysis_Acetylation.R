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
source("Rfunctions/anova_pep_prot_CD.R")
source('Rfunctions/tukey_range_test.R')
source('Rfunctions/get_combinations_CD.R')
source("Rfunctions/semi_to_fully_CD.R")
source("Rfunctions/barplot_numbers_ggplot.R")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Color Palette
#-----------------------------------------------------------------------------------------------------------------------------------------------------
pal <- wes_palette("Darjeeling1", 40, type = "continuous")
pal_dis <- wes_palette("Darjeeling1", 5, type = "discrete")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Parameters
#-----------------------------------------------------------------------------------------------------------------------------------------------------
conditions <-  c("WT_Control", "WT_HU")
n <- length(conditions)
filename <- "Acetylation_WT_control_HU"
import_filename <- "Acetylation_WT_control_HU"

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data
#-----------------------------------------------------------------------------------------------------------------------------------------------------
acetyl_data <- read.table(paste("./FilteredData/",  import_filename, "_acetyl_data_filtered.tsv", sep = ""), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tc_data <- read.delim(paste("./FilteredData/", import_filename, "_tc_data_filtered.tsv", sep = ""), sep = ";", header = TRUE, stringsAsFactors = FALSE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the mean peptide/protein ratio and the corresponding standard error
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# subset the datasets to only contain relevant information
lip_sub <- acetyl_data[, c(as.vector(sapply(conditions, function(x) return(paste(x, c(1:4), sep = "_")))),"mod_Peptide", "Protein") ]
tc_sub <- tc_data[, c(as.vector(sapply(conditions, function(x) return(paste(x, c(1:4), sep = "_")))), "Protein") ]

# replace the proteins that have not been measured in a fraction by NA
tc_sub[tc_sub == 0] <- NA

# store the mean peptides/protein intensities and the corresponding propagated error in mean_sd for each fraction
mean_peptides_proteins <- data.frame(row.names = lip_sub$mod_Peptide)
sd_peptides_proteins <- data.frame(row.names = lip_sub$mod_Peptide)
nb_replicates <- data.frame(row.names = lip_sub$mod_Peptide)

# loop over each fraction/condition
for (i in c(0:(n-1))){
  
  # get the peptides and proteins measured in one fraction
  peptides <- lip_sub[, c((i*4 +1) : (i*4 + 4), ncol(lip_sub) -1 , ncol(lip_sub))]
  proteins <- tc_sub[, c((i*4 +1) : (i*4 + 4), ncol(tc_sub))]
  
  # use the protein name as row names
  rownames(proteins) <- proteins[,ncol(proteins)]
  proteins <- proteins[, -ncol(proteins)]
  
  replicates <- c()
  mean_pep_prot <- c()
  propagated_sd <- c() 
  
  # for each peptide check whether it was measured in 3 or 4 replicates and calculate the mean and sd based on this
  for(j in c(1:nrow(peptides))){
    # get the number of replicates for which a peptide has been measured
    current_peptide <- peptides[j, ]
    index_replicates <- which(current_peptide[,c(1:4)]!=0)
    # store in how many repliates a peptide has been measured
    replicates <- c(replicates,  length(index_replicates))
    
    # if the peptide was measured at least in triplicates compute the mean and sd
    if(length(index_replicates) >= 3){
      
      mean_peptides <- mean(as.numeric(current_peptide[, index_replicates]))
      mean_proteins <- mean(as.numeric(proteins[current_peptide$Protein, index_replicates ]))
      
      mean_pep_prot_current <- mean_peptides / mean_proteins
      mean_pep_prot <- c(mean_pep_prot, mean_pep_prot_current)
      
      sd_peptides <- sd(current_peptide[, index_replicates])
      sd_proteins <- sd(proteins[current_peptide$Protein, index_replicates ])
      
      propagated_sd_current <- mean_pep_prot_current * sqrt((sd_peptides/mean_peptides)^2 + (sd_proteins/mean_proteins)^2)
      propagated_sd <- c(propagated_sd, propagated_sd_current)
      
    } else{
      mean_pep_prot <- c(mean_pep_prot, 0)
      propagated_sd <- c(propagated_sd, 0)
    }
  }
  
  mean_peptides_proteins <- cbind.data.frame(mean_peptides_proteins, mean_pep_prot)
  sd_peptides_proteins <- cbind.data.frame(sd_peptides_proteins, propagated_sd)
  nb_replicates <- cbind.data.frame(nb_replicates, replicates)
}

# adjust column names
colnames(mean_peptides_proteins) <- conditions
colnames(sd_peptides_proteins) <- conditions
colnames(nb_replicates) <- conditions

# replace 0 with NA for violin plots
violin_df_mean_peptide_protein <- mean_peptides_proteins
violin_df_mean_peptide_protein[violin_df_mean_peptide_protein == 0] <- NA 

violin_df_sd_peptides_proteins <- sd_peptides_proteins
violin_df_sd_peptides_proteins[violin_df_sd_peptides_proteins == 0] <- NA 

violin_plot(violin_df_mean_peptide_protein, my_title = "Mean Peptide/Protein", my_ylab = "log10(mean(pep/prot))", my_filename = paste("mean_peptide_proteins_", filename, sep =""))
violin_plot(violin_df_sd_peptides_proteins, my_title = "Standard Deviation Peptide/Protein", my_ylab = "log10(sd(pep/prot))", my_filename = paste("sd_peptide_proteins_", filename, sep =""))

list_pep_prot_mean_sd <- list(mean_peptides_proteins, sd_peptides_proteins, nb_replicates)
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
# iterate over each peptide and calculate the F statistic
F_values <- c()
p_values <- c()
MS_errors <- c()
MS_groups <- c()
L2FCs <- c()

# ANOVA: calculate F value and corresponding probability 
for (i in c(1:nrow(mean_peptides_proteins))){
  
  # extract the means and the standard deviations of the peptide/protein ratio in all fractions for the current peptide
  means <- as.vector(as.numeric(mean_peptides_proteins[i, ]))
  sds <- as.vector(as.numeric(sd_peptides_proteins[i, ]))
  rep <- as.vector(as.numeric(nb_replicates[i, ]))
  
  # if index != 0 then peptitde was not measured in this fraction
  index <- sapply(means, function(x) length(which(x == 0)))
  index <- which(index == 0)
  current_means <- means[index]
  current_sds <- sds[index]
  current_rep <- rep[index]
  
  # define the number of groups that we compare (only take the groups that have non-zero peptide intensities)
  nb_groups <- length(current_means)
  total_nb_observations <- sum(current_rep)
  
  # mean of means
  grand_mean <- mean(current_means)
  # within standard error
  MS_group <- sum(sapply(current_means, function(x) (x - grand_mean)^2) * current_rep)/ (nb_groups-1)
  # between standard error
  MS_error <-  sum(sapply(current_sds, function(x) x^2) * (current_rep-1)) / (total_nb_observations - nb_groups) 
  
  # calculate the F value and the corresponding probability to observe such an event or a more extreme event
  F_value <- MS_group/MS_error
  p_value <- pf(F_value, nb_groups -1 , total_nb_observations - nb_groups, lower.tail = FALSE)

  # store the variables
  F_values <- c(F_values, F_value)
  p_values <- c(p_values, p_value)
  MS_errors <- c(MS_errors, MS_error)
  MS_groups <- c(MS_groups, MS_group)
  L2FCs <- c(L2FCs, log2(current_means[2]/current_means[1]))
  
}

# correct for multiple testing
q_values <- p.adjust(p_values, method = "BH")

# store the summary of the results in anova_summary
anova_summary <- cbind(rownames(mean_peptides_proteins), F_values, p_values, q_values, L2FCs,  MS_errors, MS_groups, mean_peptides_proteins)
colnames(anova_summary) <- c("mod_Peptide", "F", "p", "q", "L2FC", "MSError", "MSGroup", paste("Mean", conditions))
anova_summary <- merge(anova_summary, acetyl_data[, c("Protein","mod_Peptide", "Position")])
write.table(anova_summary, paste("./OutputData/", filename, "_anova_results.csv", sep = ""),  row.names = FALSE)

# annotate the dataset with gene names
gene_names <- read.csv("../Files/210713_Yeast_Uniprot_Names.csv", header = TRUE, stringsAsFactors = FALSE)
gene_names$Gene.names <- sapply(gene_names$Gene.names, function(x) return(strsplit(x, split = ' ')[[1]][1]))

anova_summary_named <- merge(anova_summary, gene_names, by.x = "Protein", by.y = "Entry")
anova_summary_named <- merge(anova_summary_named, acetyl_data[, c("mod_Peptide", "Peptide")])
anova_summary_sign <- anova_summary[which(anova_summary$q < 0.01), ]

# add information aboud whether the changing acetylation protein is a known gcn5 target
gcn5_targets <- read.delim("../Files/GCN5_acetylation_targets.csv", header = TRUE, stringsAsFactors = FALSE, sep = ',')
gcn5_targets <- unique(c(gcn5_targets$GENE, "HHF1", "NUP60", "UME6"))

anova_summary_named$saga_regulated <- rep(FALSE, nrow(anova_summary_named))
anova_summary_named[which(anova_summary_named$Gene.names %in% gcn5_targets), "saga_regulated"] <- TRUE

anova_summary_named[which(anova_summary_named$Gene.names == "HHF1"), "Gene.names"] <- "H4"
anova_summary_named[which(anova_summary_named$Gene.names == "HHT1"), "Gene.names"] <- "H3"

protein_complexes <- read.delim("../3_Network_Analysis_WT_control_HU/final/Complex_Summary_WT_control_HU.tsv", header = TRUE, stringsAsFactors = FALSE)

anova_summary_named$complex_changes <- sapply(anova_summary_named$Protein, function(x){
  if(x %in% protein_complexes$Protein){
    return(protein_complexes[grep(x, protein_complexes$Protein), "Complex"])
  }else{
    return(0)
  }
})
write.csv(anova_summary_named, paste("./OutputData/", filename, "_anova_results_named.csv", sep = ""),  row.names = FALSE, quote = FALSE)

keyvals <- ifelse(
  anova_summary_named$saga_regulated == TRUE, '#a7171a',
  ifelse(anova_summary_named$q < 0.1 & abs(anova_summary_named$L2FC) > 0.5, '#0057B7',
         'grey'))

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == '#a7171a'] <- 'Gcn5 Target'
names(keyvals)[keyvals == '#0057B7'] <- 'Non Gcn5 Target'
 
keyvals.shape <- ifelse(anova_summary_named$complex_changes != 0, 8,16)
names(keyvals.shape)[keyvals.shape == 8] <- 'PPI Change'
names(keyvals.shape)[keyvals.shape == 16] <- 'No PPI Change'

EnhancedVolcano(anova_summary_named, x = "L2FC", y = "q",  ylab = bquote(~-log[10]~italic(q)), pCutoff = 0.1, FCcutoff = 1,
                ylim = c(0, 2.1), xlim = c(-21, 21), 
                colCustom = keyvals,
                shapeCustom = keyvals.shape,
                subtitle = "", pointSize = 2, gridlines.major = FALSE, gridlines.minor = FALSE, title = "", shape = 16, 
                #col = c('grey75', 'grey75','grey75', "#0057b7"),
                lab = anova_summary_named$Gene.names, labSize =  2, legendPosition = 'none')
ggsave(paste("./Plots/", filename, "_volcano_anova.pdf", sep =""),  device = "pdf", width = 14 , height = 16, units = "cm")