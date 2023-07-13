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

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Color Palette
#-----------------------------------------------------------------------------------------------------------------------------------------------------
library("wesanderson")
pal <- wes_palette("Darjeeling1", 40, type = "continuous")
library(ggsci)
library(RColorBrewer)

colfunc<-colorRampPalette(c("#b0c7e0","#0057b7","#001e64"))
colfunc_inverse <- colorRampPalette(c("#001e64","#0057b7","#b0c7e0"))
colfunc_grey_blue <- colorRampPalette(c( "#001e64","#b0c7e0", "grey95"))

plot(rep(1,50),col=(colfunc(50)), pch=19,cex=2)

plot(rep(1,50),col=(colfunc_grey_blue(50)), pch=19,cex=2)


# highlights: yellow: #ffc900 red:#a7171a
# greys: grey75 
# colour gradient: barplot(1:12, col = wes_palette("FantasticFox1", 12, type = "continuous"))

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Functions loaded from files in Rfunctions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/semi_fully_tryptic.R")
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/violin_plot.R")
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/triplicates_lip.R")
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/triplicates_tc.R")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Input Data
#-----------------------------------------------------------------------------------------------------------------------------------------------------

lip_100K <- read.delim("./Input/20220203_093918_CD010_100K_LiP_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
lip_50K <- read.delim("./Input/20220203_100025_CD010_50K_LiP_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
lip_30K <- read.delim("./Input/20220203_100458_CD010_30K_LiP_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
lip_10K <- read.delim("./Input/20220203_102911_CD010_10K_LiP_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

tc_100K <- read.delim("./Input/20220203_094101_CD010_100K_TC_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
tc_50K <- read.delim("./Input/20220203_100054_CD010_50K_TC_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
tc_30K <- read.delim("./Input/20220203_102849_CD010_30K_TC_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
tc_10K <- read.delim("./Input/20220203_102935_CD010_10K_TC_Report.xls", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

colnames(lip_100K) <- c("Protein", "Peptide", "Proteotypic", "Position", "tryptic", "100K2", "100K3", "100K4", "100K1")
colnames(lip_50K) <- c("Protein", "Peptide", "Proteotypic", "Position", "tryptic", "50K1", "50K3", "50K2", "50K4")
colnames(lip_30K) <- c("Protein", "Peptide", "Proteotypic", "Position", "tryptic", "30K4", "30K3", "30K2", "30K1")
colnames(lip_10K) <- c("Protein", "Peptide", "Proteotypic", "Position", "tryptic", "10K1", "10K2", "10K4", "10K3")


colnames(tc_100K) <- c("MW", "Protein", "100K2", "100K3", "100K4", "100K1")
colnames(tc_50K) <- c("MW", "Protein", "50K1", "50K3", "50K2", "50K4")
colnames(tc_30K) <- c("MW", "Protein", "30K4", "30K3", "30K2", "30K1")
colnames(tc_10K) <- c("MW", "Protein", "10K1", "10K2", "10K4", "10K3")


lip_data <- merge(lip_100K, lip_50K, all = TRUE)
lip_data <- merge(lip_data, lip_30K, all = TRUE)
lip_data <- merge(lip_data, lip_10K, all = TRUE)


tc_data <- merge(tc_100K, tc_50K, all = TRUE)
tc_data <- merge(tc_data, tc_30K, all = TRUE)
tc_data <- merge(tc_data, tc_10K, all = TRUE)

conditions <-  c("100K", "50K", "30K", "10K")
n <- length(conditions)
filename <- "FLiPRFractions"




#only uniquely identified proteins
tc_data <- tc_data[-grep(";", tc_data$Protein), ]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# LiP-Data Filtering: keep only peptides that have been measured in at least triplicates
#-----------------------------------------------------------------------------------------------------------------------------------------------------


lip_data <- lip_data[, c("Peptide", "Protein",
                                    "100K1", "100K2", "100K3", "100K4", "50K1", "50K2", "50K3", "50K4", 
                                    "30K1", "30K2", "30K3", "30K4", "10K1", "10K2", "10K3", "10K4",
                                    "tryptic", "Position", "Proteotypic")]


lip_result_list <- triplicates_lip(lip_data, n, conditions, filename)
lip_data <- lip_result_list[[1]]
nb_identified_peptides <- lip_result_list[[2]]
identified_peptides <- lip_result_list[[3]]
nb_identified_proteins_lip <- lip_result_list[[4]]
identified_proteins_lip <- lip_result_list[[5]]

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# LiP-Data Summary: Plot of Peptides and their nature in distinct fractions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# plot the number of peptides in each fraction
# need a folder called "Plots" where all the plots are being stored
df_nb_peptides <- data.frame(fraction = c(conditions, "all"), number = c(nb_identified_peptides, length(unique(unlist(identified_peptides)))))
df_nb_peptides$fraction <- factor(df_nb_peptides$fraction, levels = df_nb_peptides$fraction)

library("extrafont")

ggplot(df_nb_peptides, aes(x = fraction, y = number, fill = as.factor(c(1:(n+1))))) + geom_bar( stat="identity", show.legend = FALSE) + ylim(c(0, 40000)) +theme_classic(base_size = 22) + xlab("") +
  ylab ("Number of Peptides") + scale_fill_manual(values = c(colfunc(5))) + geom_text(aes(label = round(number, 1)), vjust = -0.5, size = 8)

ggsave(paste("./Plots/", filename, "_Number_Peptides_LiP.pdf", sep =""), device = "pdf", width = 16 , height = 14, units = "cm", useDingbats=FALSE)


# plot the number of proteins for which LiP peptides have been identified in each fraction

df_nb_proteins <- data.frame(fraction = c(conditions, "all"), number = c(nb_identified_proteins_lip, length(unique(unlist(identified_proteins_lip)))))
df_nb_proteins$fraction <- factor(df_nb_proteins$fraction, levels = df_nb_proteins$fraction)


ggplot(df_nb_proteins, aes(x = fraction, y = number, fill = as.factor(c(1:(n+1))))) + geom_bar( stat="identity", show.legend = FALSE) + ylim(c(0, 2200)) +
  scale_fill_manual(values = colfunc(5)) + theme_classic(base_size = 22) + xlab("") + ylab ("Number of Proteins LiP") + geom_text(aes(label = round(number, 1)), vjust = -0.5, size = 8) 


ggsave(paste("./Plots/", filename, "_Number_Proteins_LiP.pdf", sep =""), device = "pdf", width = 16 , height = 14, units = "cm", useDingbats=FALSE)



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

ggplot(percentage_semi_tryptic_peptide_intensities, aes(x = fraction, y = percentage, fill = c(1:n))) + geom_bar( stat="identity", show.legend = FALSE) + ylim(c(0, 0.75))+
  scale_fill_gradientn(colours = colfunc(5)) + theme_classic(base_size = 24) + xlab("") + ylab ("Percentage")+ geom_text(aes(label = round(percentage, 2)), vjust = -0.5, size = 8) 

ggsave(paste("./Plots/", filename, "_Percentage_Semi_typtic_LiP.pdf", sep =""), device = "pdf", width = 16 , height = 14, units = "cm")


# UpSet Plot gives intersection of measured proteins between all conditions
pdf(paste("./Plots/UpSetR_",filename,  "_Peptides_LiP", ".pdf", sep = ""), width = 16, height = 14, compress = FALSE )
upset(fromList(identified_peptides), order.by = "freq", nintersects = 10, text.scale = 5, keep.order = TRUE, main.bar.color = colfunc_inverse(15), sets.bar.color = colfunc(4), 
      set_size.show = FALSE, point.size = 10, line.size = 2)
dev.off()

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# TC-Data Filtering: keep only proteins that have been measured in at least triplicates
#-----------------------------------------------------------------------------------------------------------------------------------------------------

tc_data <- tc_data[, c( "Protein", "MW", 
                        "100K1", "100K2", "100K3", "100K4",
                        "50K1", "50K2", "50K3", "50K4", "30K1", "30K2", "30K3", "30K4",
                        "10K1", "10K2", "10K3", "10K4")]

tc_result_list <- triplicates_tc(tc_data, n, conditions, filename)

tc_data <- tc_result_list[[1]]
nb_identified_proteins <- tc_result_list[[2]]
identified_proteins <- tc_result_list[[3]]

if(length(which(duplicated(tc_data$Protein))) > 0){
  tc_data <- tc_data[-which(duplicated(tc_data$Protein)), ]
}


write.table(tc_data, paste("./FilteredData/", filename, "_tc_data_filtered.tsv", sep = ""),  row.names = FALSE, quote = FALSE, sep = ";")
tc_data <- read.delim("./FilteredData/FLiPRFractions_tc_data_filtered.tsv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# TC Data Summary Plots
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# number of proteins per fraction
df_nb_proteins <- data.frame(fraction = c(conditions, "all") , number = c(nb_identified_proteins, length(unique(unlist(identified_proteins)))))
df_nb_proteins$fraction <- factor(df_nb_proteins$fraction, levels = df_nb_proteins$fraction)

# number of proteins that has been identified in each fraction
ggplot(df_nb_proteins, aes(x = fraction, y = number, fill = c(1:(n+1)))) + geom_bar( stat="identity", show.legend = FALSE) + ylim(c(0, 2700))+
  scale_fill_gradientn(colours = colfunc(5)) + theme_classic(base_size = 24) + xlab("") + ylab ("Number of Proteins") + geom_text(aes(label = round(number, 1)), vjust = -0.5, size = 8) 

ggsave(paste("./Plots/", filename, "_Number_Proteins.pdf", sep = ""), device = "pdf", width = 16 , height = 14, units = "cm")

# heatmap of proteins in each fraction (! even though tc filtered data is exactly the same as for Master Thesis Analysis the Heatmap looks different)
jpeg(paste("./Plots/", filename, "Protein_Heatmap.pdf", sep =""),  width = 24, height = 14,units = "cm", res = 1000)
heatmap(as.matrix(tc_data[, c(3: (ncol(tc_data)-1)) ]), col =  colfunc_grey_blue(6))
dev.off()

# UpSet Plot gives intersection of measured proteins between all conditions
pdf(paste("./Plots/UpSetR_",filename, "_Proteins_TC.pdf", sep = ""), width = 16, height = 14, compress = FALSE)
upset(fromList(identified_proteins), order.by = "freq", nintersects = 10, text.scale = 5, keep.order = TRUE, main.bar.color = colfunc_inverse(15), sets.bar.color = colfunc(4), 
      set_size.show = FALSE, point.size = 10, line.size = 2)
dev.off()


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# TC Data Molecular Weight Distribution Across Fractions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# summed proteins intensities over each bin for each condition
ggplots <- list()

# classify proteins based on bin sizes
level_intervals <- levels(tc_data$bin)
level_intervals <- unlist(strsplit(level_intervals, "\\["))
level_intervals <- unlist(strsplit(level_intervals, "\\("))
level_intervals <- unlist(strsplit(level_intervals, "\\]"))

#get the mean value of every interval
mean_level <- sapply(level_intervals, function(x) mean(as.numeric(unlist((strsplit(x, ","))))))
mean_level <- round(mean_level, digits = 0)

tc_list <- list()
for(i in c(1:n)){
  tc_list[[i]] <- tc_data[, c(1,2, ((i-1)*4+3):((i-1)*4+6), ncol(tc_data))]
  # each fraction has a different total protein intensity --> to compare between fractions we need to normalize by the total protein quantity in the respective fraction
  total_intensity <- sum(apply(tc_list[[i]][,c(3:6)], 1,  function(x) return(mean(x))))
  tc_list[[i]]$mean <- apply(tc_list[[i]][,c(3:6)], 1,  function(x) return(mean(x)/total_intensity))
}

for(i in c(1:length(tc_list))){
  df <- tc_list[[i]]
  
  summed_LFQ <-  sapply(levels(df$bin), function (x) sum(rowSums(df[which(df$bin == x), c(3:6)])))
  summed_LFQ <- as.data.frame(summed_LFQ)
  
  interval_means <- as.data.frame((strsplit(str_sub(levels(df$bin), 2, nchar(levels(df$bin))-1 ), ",")))
  interval_means <- lapply(interval_means,as.numeric)
  interval_means <- as.vector(unlist(lapply(interval_means, mean)))
  
  summed_LFQ$intervals <- interval_means

  p <- ggplot(summed_LFQ, aes(x = intervals, y = summed_LFQ, fill = intervals)) + geom_histogram(stat = "identity", show.legend = FALSE, binwidth = 5) + theme_classic(base_size = 14) +
    xlab("") + ylab("")  +  xlim(c(0, 100)) + ylim(c(0, 2.5e9 )) +  scale_fill_gradient(low  = 'lightblue', high = 'orange')
  
  ggplots[[i]] <- p 
}

intensity_figure <- ggarrange(plotlist = ggplots, ncol = 3, nrow = 2,  label.x = 0.4, font.label = list(size = 11), labels = conditions)
annotate_figure(intensity_figure, bottom = text_grob("Molecular Weight (kDa)"), left = text_grob("Summed Protein Intensities", rot = 90))
ggsave(filename = paste("./Plots/", filename, "_Protein_Separtion.pdf", sep =""), device = "pdf", width = 20, height = 10, units = "cm")


# library(rcartocolor)
# nColor <- 4
# scales::show_col(carto_pal(nColor, "Safe"))

# highlights: yellow: #ffc900 red:#a7171a
# greys: grey75 
# colour gradient: barplot(1:12, col = wes_palette("FantasticFox1", 12, type = "continuous"))
# c("#b0c7e0","#0057b7","#001e64"))


library("spatstat")
jpeg(paste("./Plots/", filename, "Protein_Separation_Density.pdf", sep =""),  width = 20, height = 20,units = "cm", res = 2000)
plot(density(tc_list[[4]]$MW, weights = tc_list[[4]]$mean), xlim = c(0,100), ylim = c(0,0.035), col =  "#ffc900", lwd = 10, 
     xlab = "", main = "", ylab = "", cex.axis = 2.5)
abline(v = weighted.median(tc_list[[4]]$MW, w =  tc_list[[4]]$mean), col = "#ffc900", lwd = 5, lty=2)
lines(density(tc_list[[3]]$MW, weights = tc_list[[3]]$mean), col = "#a7171a", lwd = 10)
abline(v = weighted.median(tc_list[[3]]$MW, w =  tc_list[[3]]$mean), col = "#a7171a", lwd = 5, lty=2)
lines(density(tc_list[[2]]$MW, weights = tc_list[[2]]$mean), col = "#0057b7", lwd = 10)
abline(v = weighted.median(tc_list[[2]]$MW, w =  tc_list[[2]]$mean), col = "#0057b7", lwd = 5, lty=2)
lines(density(tc_list[[1]]$MW, weights = tc_list[[1]]$mean), col = "gray40", lwd = 10)
abline(v = weighted.median(tc_list[[1]]$MW, w =  tc_list[[1]]$mean), col = "gray40", lwd = 5, lty=2)

dev.off()

 #-----------------------------------------------------------------------------------------------------------------------------------------------------
# Percentage of Proteins Present in Complex and Monomer
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# ! also here plot looks different then  for original analysis in Master Thesis even though data is exactly the same
# how many percent of protein are present as monomers and as complexes
percentage_complexes <- c()
for(i in c(1:nrow(tc_data))){

  MW <- tc_data[i,]$MW
  
  if(MW < 30){
    if(MW > 10){
     
      monomers <- sum(tc_data[i,c(15:18)])
      complexes <- sum(tc_data[i,c(3:10)])
    
    }else{
      
      monomers <- NA
      complexes <- NA

    }
    percentage_complexes <- c(percentage_complexes, complexes/(complexes + monomers))
  }else{
  
    percentage_complexes <- c(percentage_complexes, NA)
  }
}

tc_data$percentage_complexes <- percentage_complexes


ggplot(tc_data, aes(x = percentage_complexes, fill = MW)) + geom_histogram(fill = "grey") + theme_classic(base_size = 20) +
  xlab("Percentage Complex-Bound Form") + ylab("Number of Proteins") +ylim(c(0,100))
ggsave(filename = paste("./Plots/", filename, "_Percentage_Complex_Bound_Form_yeast.pdf", sep = ""),  device = "pdf", width = 12, height = 10, units = "cm")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Subset LiP Data to only contain the Proteins present in the TC data and safe
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# merge lip data with the molecular weight of the tryptic control data - like this we can get rid of peptides from proteins which have not been quantified in the tryptic control
lip_data <- merge(lip_data, tc_data[, c(1,2)])
lip_data <- unique(lip_data)


# export the normalized and filtered lip_data
write.table(lip_data, paste("./FilteredData/", filename, "_lip_data_filtered.tsv", sep = ""),  row.names = FALSE, quote = FALSE, sep = "\t")



