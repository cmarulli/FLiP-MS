library("reshape2")
library("ggplot2")
library("EnhancedVolcano")
source("/Users/cathy/Documents/PhD/SEC-F-LiP/Rfunctions/violin_plot.R")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# define variables
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bait <- "P32494"

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read in data
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data <- read.csv("./Input/20230117_nvolkmar_EX140_AP_Protein_MS1_Intensities.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(data) <- c("Protein",  "AP_untreated_1", "AP_untreated_2", "AP_untreated_3", "AP_untreated_4", 
                    "AP_HU_1", "AP_HU_2", "AP_HU_3", "AP_HU_4", 
                    "control_untreated_1", "control_untreated_2", "control_untreated_3", "control_untreated_4", 
                    "control_HU_1", "control_HU_2", "control_HU_3", "control_HU_4")

# get the high-confidence interactors from SAINT scoring
SAINT_HU <- read.table("./Input/NV140_Ada3_APMS_HU_SAINT_results.tab", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
SAINT_untreated <- read.table("./Input/NV140_Ada3_APMS_untreated_SAINT_results.tab", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

SAINT_HU_interactor <- SAINT_HU[which(SAINT_HU$SP == 1), ]
SAINT_untreated_interactor <- SAINT_untreated[which(SAINT_untreated$SP == 1), ]
all_interactor <- unique(c(SAINT_HU_interactor$prey, SAINT_untreated_interactor$prey))

# only keep the quantification of high confidence interactors
data <- data[which(data$Protein %in% all_interactor), ]
data[is.na(data)] <-  runif(sum(is.na(data)), min = 0, max = 200000)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SAGA/SLIK subunits
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# extract information about which proteins are part of the SAGA and SLIK complex
SAGA_complex_portal <- read.csv("./Input/SAGA_SLIK.csv", header = TRUE, stringsAsFactors = FALSE)
SAGA <- SAGA_complex_portal[1,"Expanded.participant.list"] 
SAGA <- strsplit(SAGA, split = c("\\|"))
SAGA <- sapply(SAGA, function(x) substr(x, start = 1, stop = 6))

SLIK <- SAGA_complex_portal[2,"Expanded.participant.list"] 
SLIK <- strsplit(SLIK, split = c("\\|"))
SLIK <- sapply(SLIK, function(x) substr(x, start = 1, stop = 6))

SAGA_SLIK <- unique(c(SAGA, SLIK))
SAGA_SLIK <- c(SAGA_SLIK, "Q12433", "P25649")

protein_names <- read.csv("../../../Databases/210713_Yeast_Uniprot_Names.csv", header = TRUE, stringsAsFactors = FALSE)
protein_names$name <- sapply(protein_names$Gene.names, function(x) strsplit(x, split = " ")[[1]][1])

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# noramlization the data by total area sums
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# intensities prior to normalization
violin_plot(data[, c(2:17)], my_title = "Protein Intensities before Normalization", my_ylab = "log10(Intensity)", my_filename = "ProteinIntensities_no_norm")

# Total Area Sums
TAS <- apply(data[, c(2:17)], 2, function(x) return(sum(x, na.rm = TRUE)))
TAS_ratios <- TAS/max(TAS)
data_norm <- sweep(as.matrix(data[, c(2:17)]),2,TAS_ratios,FUN="/")

data_norm <- data.frame(cbind("Protein" = data[,1], data_norm))
for(i in c(2:17)){
  data_norm[, i] <- as.numeric(data_norm[,i])
}

# intensities after normalization
violin_plot(data_norm[, c(2:17)], my_title = "Protein Intensities after Normalization", my_ylab = "log10(Intensity)", my_filename = "ProteinIntensities_norm")

# separate the dataframe by condition
AP_untreated <- data_norm[, c(2:5)]
AP_HU <- data_norm[, c(6:9)]
control_untreated <- data_norm[, c(10:13)]
control_HU <- data_norm[, c(14:17)]

data_SAGA <- data[which(data$Protein %in% SAGA_SLIK), ]

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# calculate the mean and standard deviation of the respective replicates and keep only proteins measured in at least triplicates
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
AP_untreated$Mean <- apply(AP_untreated, 1, function(x) mean(x, na.rm = TRUE))
AP_untreated$SD <- apply(AP_untreated, 1, function(x) sd(x, na.rm = TRUE))
AP_untreated$n <- apply(AP_untreated, 1, function(x) length(which(x[1:4]!=0)))
rownames(AP_untreated) <- data[,1]
AP_untreated <- AP_untreated[which(AP_untreated$n >=3), ]

AP_HU$Mean <- apply(AP_HU, 1, function(x) mean(x, na.rm = TRUE))
AP_HU$SD <- apply(AP_HU, 1, function(x) sd(x, na.rm = TRUE))
AP_HU$n <- apply(AP_HU, 1, function(x) length(which(x[1:4]!=0)))
rownames(AP_HU) <- data[,1]
AP_HU <- AP_HU[which(AP_HU$n >=3), ]

control_untreated$Mean <- apply(control_untreated, 1, function(x) mean(x, na.rm = TRUE))
control_untreated$SD <- apply(control_untreated, 1, function(x) sd(x, na.rm = TRUE))
control_untreated$n <- apply(control_untreated, 1, function(x) length(which(x[1:4]!=0)))
rownames(control_untreated) <- data[,1]
control_untreated <- control_untreated[which(control_untreated$n >=3), ]

control_HU$Mean <- apply(control_HU, 1, function(x) mean(x, na.rm = TRUE))
control_HU$SD <- apply(control_HU, 1, function(x) sd(x, na.rm = TRUE))
control_HU$n <- apply(control_HU, 1, function(x) length(which(x[1:4]!=0)))
rownames(control_HU) <- data[,1]
control_HU <- control_HU[which(control_HU$n >=3), ]

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# for each SAGA/SLIK component calculate the mean ratio of prey to bait and the respective propagated error in control and HU-sress
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prey_bait_untreated_mean <- c()
prey_bait_untreated_sd <- c()

prey_bait_HU_mean <- c()
prey_bait_HU_sd <- c()

#for(prey in SAGA_SLIK){
for(prey in rownames(AP_untreated)){    
  
  current_prey_bait_untreated_mean <- (AP_untreated[prey, "Mean"] / control_untreated[prey, "Mean"]) / (AP_untreated[bait, "Mean"] / control_untreated[bait, "Mean"])
  current_prey_bait_untreated_sd <- current_prey_bait_untreated_mean * sqrt(  (AP_untreated[prey, "SD"]/AP_untreated[prey, "Mean"])^2 + (control_untreated[prey, "SD"]/control_untreated[prey, "Mean"])^2 +
                                                                                (AP_untreated[bait, "SD"]/AP_untreated[bait, "Mean"])^2 + (control_untreated[bait, "SD"]/control_untreated[bait, "Mean"])^2)
  
  prey_bait_untreated_mean <- c(prey_bait_untreated_mean, current_prey_bait_untreated_mean)
  prey_bait_untreated_sd <- c(prey_bait_untreated_sd, current_prey_bait_untreated_sd)
  
  current_prey_bait_HU_mean <- (AP_HU[prey, "Mean"] / control_HU[prey, "Mean"]) / (AP_HU[bait, "Mean"] / control_HU[bait, "Mean"])
  current_prey_bait_HU_sd <- current_prey_bait_HU_mean * sqrt(  (AP_HU[prey, "SD"]/AP_HU[prey, "Mean"])^2 + (control_HU[prey, "SD"]/control_HU[prey, "Mean"])^2 +
                                                                  (AP_HU[bait, "SD"]/AP_HU[bait, "Mean"])^2 + (control_HU[bait, "SD"]/control_HU[bait, "Mean"])^2)
  
  prey_bait_HU_mean <- c(prey_bait_HU_mean, current_prey_bait_HU_mean)
  prey_bait_HU_sd <- c(prey_bait_HU_sd, current_prey_bait_HU_sd)
  
}

mean_sd_df <- data.frame(cbind("Protein" = rownames(AP_untreated), "Untreated_Mean" = prey_bait_untreated_mean, "HU_Mean" = prey_bait_HU_mean, 
                               "Untreated_SD" = prey_bait_untreated_sd, "HU_SD" = prey_bait_HU_sd))

for(i in c(2:5)){
  mean_sd_df[, i] <- as.numeric(mean_sd_df[,i])
}

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# calculate the p-value and L2FC
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_sd_df$t <- apply(mean_sd_df, 1, function(x) return( (as.numeric(x[3])-as.numeric(x[2])) / (sqrt((as.numeric(x[4])^2 + as.numeric(x[5])^2) / 2) * sqrt(2/3))))
mean_sd_df$L2FC <- apply(mean_sd_df, 1, function(x) return(log2(as.numeric(x[3])/as.numeric(x[2]))))
mean_sd_df$p <- apply(mean_sd_df, 1, function(x) return(pt(q = abs(as.numeric(x[6])), df = 4, lower.tail = FALSE)))

mean_sd_df <- merge(mean_sd_df, protein_names[, c("Entry", "name")], by.x = "Protein", by.y = "Entry")
rownames(mean_sd_df) <- mean_sd_df$name

keyvlas <- ifelse(
  mean_sd_df$Protein %in% SAGA_SLIK, '#0057B7', 'grey70')
names(keyvlas)[keyvlas == '#0057B7'] <- 'SAGA'
names(keyvlas)[keyvlas == 'grey70'] <- ''


EnhancedVolcano(mean_sd_df, x = "L2FC", y = "p", lab = rownames(mean_sd_df), ylab = bquote(~-log[10]~italic(p)), pCutoff = 0.06, FCcutoff = 1,
                ylim = c(0,1.8), xlim = c(-2, 2), subtitle = "", pointSize = 2, gridlines.major = FALSE, gridlines.minor = FALSE, title = "", colCustom = keyvlas, 
                colAlpha = c(ifelse(mean_sd_df$Protein %in% SAGA_SLIK, 1, 0.3)))

ggsave("./Plots/volcano_prey_bait_all.pdf",  device = "pdf", width = 12 , height = 16, units = "cm")

mean_sd_df_SAGA <- mean_sd_df[which(mean_sd_df$Protein %in% SAGA_SLIK), ]

EnhancedVolcano(mean_sd_df_SAGA, x = "L2FC", y = "p", lab = rownames(mean_sd_df_SAGA), ylab = bquote(~-log[10]~italic(p)), pCutoff = 0.06, FCcutoff = 1,
                ylim = c(0,1.5), xlim = c(-2, 2), subtitle = "", pointSize = 2, gridlines.major = FALSE, gridlines.minor = FALSE, title = "")
ggsave("./Plots/volcano_prey_bait_SAGA.pdf",  device = "pdf", width = 12 , height = 16, units = "cm")






