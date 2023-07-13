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

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Information about the data to be defined each time the script is run
#-----------------------------------------------------------------------------------------------------------------------------------------------------
filename <- paste("FLiPR_NewNorm")

conditions <- c('100K', '50K', '30K', '10K')
fractions <- c("Mean100K", "Mean50K", "Mean30K", "Mean10K")
n = length(conditions)
 
lip_data <-read.table("./FilteredData/FLiPRFractions_lip_data_filtered.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tc_data <- read.delim("./FilteredData/FLiPRFractions_tc_data_filtered.tsv", sep = ";", header = TRUE, stringsAsFactors = FALSE)

lip_data$MW <- as.numeric(lip_data$MW)
lip_data <- unique(lip_data)
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Color Palette
#-----------------------------------------------------------------------------------------------------------------------------------------------------
library("wesanderson")
pal <- wes_palette("Darjeeling1", 40, type = "continuous")
pal_dis <- wes_palette("Darjeeling1", 5, type = "discrete")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/violin_plot.R")
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/calc_mean_sd_pep_prot.R")
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/anova_pep_prot.R")
source('/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/tukey_range_test.R')
source('/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/get_combinations.R')
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/semi_to_fully.R")
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/barplot_numbers_ggplot.R")
source('/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/get_GO_annotation.R')
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/fisher_GO_test.R")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Remove non-proteotypic peptides if present
#-----------------------------------------------------------------------------------------------------------------------------------------------------
if(length(which(lip_data$Proteotypic == "False")) > 0){
  lip_data <- lip_data[-which(lip_data$Proteotypic == "False"), ]
}

lip_data$Position <- as.numeric(lip_data$Position)
lip_data <- lip_data[-which(is.na(lip_data$Position)), ]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the mean peptide/protein ratio and the corresponding standard error
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# subset the datasets to only contain relevant information
lip_sub <- lip_data[, c(3:(3+4*n-1), 1, 2)]
tc_sub <- tc_data[, c(3:(3+4*n-1), 1)]

# replace the proteins that have not been measured in a fraction by NA
tc_sub[tc_sub == 0] <- NA


# calculate the peptide per protein intensity in every fraction
peptide_protein <- data.frame()

for(i in c(1:nrow(lip_sub))){
  peptide_protein <- rbind(peptide_protein, lip_sub[i, c(1:16)]/tc_sub[which(tc_sub$Protein == lip_sub[i, "Protein"]), c(1:16)])
}

# replace all 0 intensities with NA, otherwise in the median normalization, this will lead to wrong medians
peptide_protein[peptide_protein == 0] <- NA
violin_plot(peptide_protein, my_title = "Peptide/Protein", my_ylab = "log10((pep/prot))", my_filename = paste("peptide_proteins_", filename, sep =""))

# median normalize the peptide per protein values
median_all <- median(as.matrix(peptide_protein), na.rm = TRUE)
peptide_protein_norm <- as.data.frame(apply(peptide_protein, 2, function(x) return(x/median(x, na.rm = TRUE) * median_all)))

peptide_protein_norm <- cbind(peptide_protein_norm, lip_sub$Peptide)

violin_plot(peptide_protein_norm[, c(1:16)], my_title = "Norm Peptide/Protein", my_ylab = "log10((pep/prot))", 
            my_filename = paste("norm_peptide_proteins_", filename, sep =""))



list_pep_prot_mean_sd <- calc_mean_sd_pep_prot(peptide_protein_norm, lip_sub, tc_sub, n, column_names = conditions, filename = filename)
mean_peptides_proteins <- list_pep_prot_mean_sd[[1]]
sd_peptides_proteins <- list_pep_prot_mean_sd[[2]]
nb_replicates <- list_pep_prot_mean_sd[[3]]

#keep only peptides which have a mean peptide/protein ratio measured in at least 2 fractions
index <- apply(mean_peptides_proteins, 1, function(x) length(which(x != 0)))
mean_peptides_proteins <- mean_peptides_proteins[which(index >= 2), ]
sd_peptides_proteins <- sd_peptides_proteins[which(index >= 2), ]
nb_replicates <- nb_replicates[which(index >= 2), ]

covered_peptides <- rownames(mean_peptides_proteins)
covered_peptides <- lip_data[which(lip_data$Peptide %in% covered_peptides), ]
covered_peptides <- covered_peptides[, c("Peptide", "Protein", "Position")]

write.csv(covered_peptides, "./OutputData/MS_detected_peptides_2fractions.csv", row.names = FALSE, quote = FALSE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# ANOVA test and Tukey's Range Test 
#-----------------------------------------------------------------------------------------------------------------------------------------------------
anova_summary <- anova_pep_prot(mean_peptides_proteins, sd_peptides_proteins, nb_replicates, lip_data, fractions)

tukey_list <- tukey_range_test(anova_summary, n, nb_replicates, conditions)

anova_summary <- tukey_list[[1]]
tukey_q <- tukey_list[[2]]
turkey_L2FC <- tukey_list[[3]]
anova_summary <- merge(anova_summary, lip_data[,c('Peptide', 'Position', 'Protein')])

# add the information about the fully tryptic parent of a semi-tryptic peptide
anova_summary <- semi_to_fully(anova_summary)
write.table(anova_summary, paste("./OutputData/", filename, "_anova_results.csv", sep = ""),  row.names = FALSE, sep = ";")

anova_summary <- read.delim(paste("./OutputData/", filename, "_anova_results.csv", sep = ""), sep = ";",  header = TRUE, stringsAsFactors = FALSE)

keyvals <- ifelse(
  anova_summary$tryptic == "Specific", '#B0C7E0',
  ifelse(anova_summary$tryptic != "Specific", '#0057B7',
         'black'))
keyvals[is.na(keyvals)] <- 'black'

names(keyvals)[keyvals == '#B0C7E0'] <- 'fully'
names(keyvals)[keyvals == '#0057B7'] <- 'semi'

EnhancedVolcano(anova_summary, x = "L2FC", y = "q", lab = "", ylab = bquote(~-log[10]~italic(q)), pCutoff = 0.05, FCcutoff = 0,
                colCustom = keyvals, ylim = c(0,8), xlim = c(-11, 11), subtitle = "", pointSize = 0.5, gridlines.major = FALSE, gridlines.minor = FALSE, title = "")
ggsave(paste("./Plots/", filename, "_volcano_anova.pdf", sep =""),  device = "pdf", width = 12 , height = 16, units = "cm")

keyvals <- rep('#0057b7', nrow(anova_summary))
names(keyvals)[keyvals == '#0057b7'] <- 'all'
EnhancedVolcano(anova_summary, x = "L2FC", y = "q", lab = "", ylab = bquote(~-log[10]~italic(q)), pCutoff = 0.05, FCcutoff = 0,
                colCustom = keyvals, ylim = c(0,8), xlim = c(-11, 11), subtitle = "", pointSize = 0.5, gridlines.major = FALSE, gridlines.minor = FALSE, title = "")
ggsave(paste("./Plots/", filename, "_volcano_anova_one_color.pdf", sep =""),  device = "pdf", width = 12 , height = 16, units = "cm")


anova_significiant <- anova_summary[which(anova_summary$q < 0.05), ]
nb_sign_peptides <- length(unique(anova_significiant$Peptide))
nb_sign_proteins <- length(unique(anova_significiant$Protein))
nb_sign_peptides_grouped <- length(unique(anova_significiant$fully_tryptic))

df <- data.frame(conditions = c("Peptides", "Positions", "Proteins"), number = c(nb_sign_peptides, nb_sign_peptides_grouped, nb_sign_proteins))
df$conditions <- factor(df$conditions, levels = df$conditions)

ggplot(df, aes(x = conditions, y = number, fill = c(1:length(conditions)))) + geom_bar(stat="identity", show.legend = FALSE) +
    scale_fill_gradientn(colours = c("#b0c7e0","#0057b7","#001e64")) + theme_classic(base_size = 28) + xlab("") + ylab("") +
  geom_text(aes(label = round(number, 1)), vjust = -0.5, size = 8) + ylim(0, 9500)

ggsave("./Plots/FLiP_Library_Summary.jpg", device = "jpg", width = 16 , height = 14, units = "cm")




#-----------------------------------------------------------------------------------------------------------------------------------------------------
# ANOVA group by peptides and overlapping peptides
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# if there are multiple overlapping fully tryptic peptides (missed cleavages), assign them to the longest fully tryptic parent
ft_peptides <- unique(anova_summary$fully_tryptic)
for(ft in ft_peptides){
  if(length(grep(ft, ft_peptides))  > 1){

    overlapping_ft_pep <- ft_peptides[grep(ft, ft_peptides)]
    # if there is indeed another ft peptide from the same protein that is overlapping replace it with the longest one
    if(length(unique(overlapping_ft_pep)) > 1 && length(unique(anova_summary[which(anova_summary$fully_tryptic == ft), ]$Protein)) == 1){
      anova_summary[which(anova_summary$fully_tryptic == ft), ]$fully_tryptic <-  overlapping_ft_pep[which.max(nchar(overlapping_ft_pep))]
    }
  }
}

# group overlapping peptides by median or min q-value



library("dplyr")
anova_median_p <- anova_summary %>% group_by(fully_tryptic) %>% summarize(median_p = median(as.numeric(p), na.rm = TRUE))
anova_min_p <- anova_summary %>% group_by(fully_tryptic) %>% summarize(min_p = min(as.numeric(p), na.rm = TRUE))


anova_summary_grouped <- merge(anova_summary, anova_median_p)
anova_summary_grouped <- merge(anova_summary_grouped, anova_min_p)

anova_summary_grouped <- unique(anova_summary_grouped[, c("fully_tryptic", "Protein", "median_p", "min_p")])

anova_summary_grouped$median_q <- p.adjust(anova_summary_grouped$median_p, method = "BH")
anova_summary_grouped$min_q <- p.adjust(anova_summary_grouped$min_p, method = "BH")


anova_grouped_min_significiant <- anova_summary_grouped[which(anova_summary_grouped$min_q < 0.05), ]
nb_sign_peptides_grouped <- length(unique(anova_grouped_min_significiant$fully_tryptic))
nb_sign_proteins <- length(unique(anova_grouped_min_significiant$Protein))

df <- data.frame(conditions = c("Positions", "Proteins"), number = c(nb_sign_peptides_grouped, nb_sign_proteins))
df$conditions <- factor(df$conditions, levels = df$conditions)

ggplot(df, aes(x = conditions, y = number, fill = c(1:length(conditions)))) + geom_bar(stat="identity", show.legend = FALSE) +
  scale_fill_gradientn(colours = c( "slategray3", "steelblue2")) + theme_classic(base_size = 28) + xlab("") + ylab("") +
  geom_text(aes(label = round(number, 1)), vjust = -0.5, size = 8) + ylim(0, 9500)

ggsave("./Plots/FLiP_Library_Summary_grouped_min_q.jpg", device = "jpg", width = 16 , height = 14, units = "cm")



anova_grouped_median_significiant <- anova_summary_grouped[which(anova_summary_grouped$median_q < 0.05), ]
nb_sign_peptides_grouped <- length(unique(anova_grouped_median_significiant$fully_tryptic))
nb_sign_proteins <- length(unique(anova_grouped_median_significiant$Protein))

df <- data.frame(conditions = c("Positions", "Proteins"), number = c(nb_sign_peptides_grouped, nb_sign_proteins))
df$conditions <- factor(df$conditions, levels = df$conditions)

ggplot(df, aes(x = conditions, y = number, fill = c(1:length(conditions)))) + geom_bar(stat="identity", show.legend = FALSE) +
  scale_fill_gradientn(colours = c( "slategray3", "steelblue2")) + theme_classic(base_size = 28) + xlab("") + ylab("") +
  geom_text(aes(label = round(number, 1)), vjust = -0.5, size = 8) + ylim(0, 9500)

ggsave("./Plots/FLiP_Library_Summary_grouped_median_q.jpg", device = "jpg", width = 16 , height = 14, units = "cm")

write.table(anova_summary_grouped, "./OutputData/FLiPRFraction_anova_results_grouped.csv", row.names = FALSE, quote = FALSE, sep =";")



#-----------------------------------------------------------------------------------------------------------------------------------------------------
# GO Annotation from InterPro
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# GO_annotation <- get_GO_annotation(unique(anova_summary$Protein))
GO_annotation <- get_GO_annotation(unique(lip_data$Protein))
colnames(GO_annotation) <- c("Protein", "IPRs", "GO_names", "GO_identifiers", "Start_GO", "End_GO", "GO_type")
lip_data <- merge(lip_data, GO_annotation, by = 'Protein', all = TRUE)

#anova_summary <- merge(anova_summary, GO_annotation, by = "Protein", all = TRUE)

# Reduce the dataset to GO-annotated peptides
lip_data <- lip_data[-which(is.na(lip_data$GO_type)), ]


lip_data$Position <- as.numeric(lip_data$Position)
#lip_data <- lip_data[-which(is.na(lip_data$Position)), ]
lip_data$index <- c(1:nrow(lip_data))
lip_data$PositionStop <- lip_data$Position + nchar(lip_data$Peptide)

# Keep only the GO annotations for which the location in the protein overlaps with the detected peptide 
peptides <- unique(lip_data$Peptide)
index <- c()
for(i in c(1:length(peptides))){
  
  current_peptide <- lip_data[which(lip_data$Peptide == peptides[i]), c("Position", "PositionStop", "Start_GO", "End_GO", "index")]
  intersect <- as.vector(apply(current_peptide, 1, function(x) length(intersect(c(x[1]:x[2]), c(x[3]:x[4])))))
  
  if(sum(intersect) > 1){
    index <- c(index, current_peptide[which(intersect!=0), "index"])
  }
}

lip_data <- lip_data[which(lip_data$index %in%  index), ]
write.table(lip_data, paste("./OutputData/", filename, "_lip_data_GO.csv", sep =""), quote = FALSE, row.names = FALSE, sep = ";")



lip_data_GO <- read.delim("./OutputData/FLiPRFraction_lip_data_GO.csv", sep = ";",  header = TRUE, stringsAsFactors = FALSE)
lip_data_GO <- lip_data

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# GO Enrichment Analysis
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# add a binary column for significance
lip_data_GO$sign <- rep(FALSE, nrow(lip_data_GO))
anova_sign <- anova_summary[which(anova_summary$q < 0.05), ]
lip_data_GO[which(lip_data_GO$Peptide %in% anova_sign$Peptide),]$sign <- TRUE

p_values_GO_terms <- fisher_GO_test(lip_data_GO)
p_values_GO_terms <- p_values_GO_terms[which(p_values_GO_terms$p_values < 0.05), ]
p_values_GO_terms <- p_values_GO_terms[order(p_values_GO_terms$p_values, decreasing = TRUE), ]

p_values_GO_terms$count <- sapply(p_values_GO_terms$GO_Terms, function(x) length(which(lip_data_GO$GO_names == x & lip_data_GO$sign == TRUE)))
p_values_GO_terms$GO_Terms <- factor(p_values_GO_terms$GO_Terms, levels = p_values_GO_terms$GO_Terms)

ggplot(p_values_GO_terms, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 10) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/", filename, "_GO_enrichment_all.jpg", sep =""),  device = "jpg", width = 24 , height = 12, units = "cm", )



# domain based only
lip_data_GO_domain <- lip_data_GO[which(lip_data_GO$GO_type %in% c("domain", "repeat")), ]
p_values_GO_terms_domain <- fisher_GO_test(lip_data_GO_domain)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.01), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(lip_data_GO_domain$GO_names == x & lip_data_GO_domain$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 12) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/", filename, "_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 24 , height = 12, units = "cm", )


ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/", filename, "_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 21 , height = 24, units = "cm", )


ggplot(p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 1.03e-06), ], aes(x= count, y = GO_Terms, fill = count)) +
  geom_bar(stat = "identity", fill = c(rep("#b0c7e0",12))) +  theme_classic(base_size = 16) + theme(axis.text.y =element_text(size=18)) +
  xlab("Significant Peptide Count") + ylab("") 
ggsave(filename = paste("./Plots/", filename, "_GO_enrichment_domain_short.pdf", sep =""),  device = "pdf", width = 35 , height = 12, units = "cm", )

# store lip_GO_data_domain for Lorenzo

lip_data_GO_domain_reduced <- lip_data_GO_domain[, c("Protein", "Peptide", "tryptic", "Position", "PositionStop", "IPRs","GO_names","GO_identifiers","Start_GO","End_GO","GO_type","sign" )]
write.table(lip_data_GO_domain_reduced, "./OutputData/FLiPR_lip_data_domain_annotation.csv", sep = ";", row.names=FALSE, quote = FALSE)

