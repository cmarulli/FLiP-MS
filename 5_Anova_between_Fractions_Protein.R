library(dplyr)
library(reshape2)
library(protti)
library("EnhancedVolcano")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Additional functions
#-----------------------------------------------------------------------------------------------------------------------------------------------------
source("./Rfunctions/violin_plot.R")
source("./Rfunctions/fisher_GO_test.R")
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data
#-----------------------------------------------------------------------------------------------------------------------------------------------------
data <- read.delim("./FilteredData/FLiP_PBI_library_tc_data_filtered.tsv", header = TRUE, stringsAsFactors = FALSE, sep = ';')
data[data == 0] <- NA
colnames(data) <-  c("Protein", "MW", "100K1", "100K2", "100K3", "100K4",
                     "50K1", "50K2", "50K3", "50K4", 
                     "30K1", "30K2", "30K3", "30K4", 
                     "10K1", "10K2", "10K3", "10K4")  

fractions <- c("100K", "50K", "30K", "10K")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Normalization
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Violin Plot of protein intensities in different fractions
violin_plot(data[, c(3:18)], my_title = "Protein_Intensities")

# median normalize protein intensities
median_all <- median(as.matrix(data[, c(3:18)]), na.rm = TRUE)
data_norm <- as.data.frame(apply(data[, c(3:18)], 2, function(x) return(x/median(x, na.rm = TRUE) * median_all)))
data_norm <- cbind(data[,c(1:2)], data_norm)
colnames(data_norm) <- colnames(data)

violin_plot(data_norm[, c(3:18)], my_title = "Protein_Intensities_Normalized")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Differential abundance analysis
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# calculate the mean intensities per fraction
data_means <- data.frame(Protein = data$Protein)
for(frac in fractions){
  data_means <- cbind(data_means, data_norm %>% rowwise() %>% summarise("mean_{frac}" := mean(c_across(paste(frac, c(1:4), sep = '')), na.rm=TRUE)))
}

# calculate the standard deviation of the intensities per fraction
data_sds <- data.frame(Protein = data$Protein)
for(frac in fractions){
  data_sds <- cbind(data_sds, data_norm %>% rowwise() %>% summarise("sd_{frac}" := sd(c_across(paste(frac, c(1:4), sep = '')), na.rm=TRUE)))
}

# calculate the replicates quantified per fraction
data_reps <- data.frame(Protein = data$Protein)
for(frac in fractions){
  data_reps <- cbind(data_reps, data_norm %>% rowwise() %>% summarise("reps_{frac}" := length(which(is.na(c_across(paste(frac, c(1:4), sep = ""))) == FALSE))))
}

# perform an ANOVA to test difference in protein intensities between filter fractions
anova_prot <- function(data_means, data_sds, data_reps,  data, fractions){
  
  # iterate over each peptide and calculate the F statistic
  F_values <- c()
  p_values <- c()
  MS_errors <- c()
  MS_groups <- c()
  
  # ANOVA: calculate F value and corresponding probability 
  for (i in c(1:nrow(data_means))){
    
    # extract the means and the standard deviations of the peptide/protein ratio in all fractions for the current peptide
    means <- as.vector(as.numeric(data_means[i, c(2:5)]))
    sds <- as.vector(as.numeric(data_sds[i, c(2:5)]))
    rep <-  as.vector(as.numeric(data_reps[i, c(2:5)]))
    
    # if index != 0 then peptide was not measured in this fraction
    index <- sapply(means, function(x) length(which(is.na(x))))
    index <- which(index == 0)
    current_means <- means[index]
    current_sds <- sds[index]
    current_rep <- rep[index]
    
    # define the number of groups that we compare (only take the groups that have non-zero peptide intensities)
    nb_groups <- length(current_means)
    total_nb_observations <- sum(current_rep)
    
    # mean of means
    grand_mean <- mean(current_means)
    # between group standard error
    MS_group <- sum(sapply(current_means, function(x) (x - grand_mean)^2) * current_rep)/ (nb_groups-1)
    # within group standard error
    MS_error <-  sum(sapply(current_sds, function(x) x^2) * (current_rep-1)) / (total_nb_observations - nb_groups) 
    
    # calculate the F value and the corresponding probability to observe such an event or a more extreme event
    F_value <- MS_group/MS_error
    p_value <- pf(F_value, nb_groups -1 , total_nb_observations - nb_groups, lower.tail = FALSE)
    
    # store the variables
    F_values <- c(F_values, F_value)
    p_values <- c(p_values, p_value)
    MS_errors <- c(MS_errors, MS_error)
    MS_groups <- c(MS_groups, MS_group)
    
  }
  
  # correct for multiple testing
  q_values <- p.adjust(p_values, method = "BH")
  
  # store the summary of the results in anova_summary
  anova_summary <- data.frame()
  anova_summary <- cbind(data_means$Protein, F_values, p_values, q_values, MS_errors, MS_groups)
  colnames(anova_summary) <- c("Protein", "F", "p", "q", "MSError", "MSGroup")
  anova_summary <- merge(anova_summary, data_means)
  
  anova_summary <- as.data.frame(anova_summary)
  
  for(i in c(2:6)){
    anova_summary[, i] <- as.numeric(anova_summary[, i])
  }
  
  anova_summary <- anova_summary %>% rowwise() %>% mutate(L2FC = max(c(log2(mean_100K/mean_50K), log2(mean_100K/mean_30K), log2(mean_100K/mean_10K), 
                                                     log2(mean_50K/mean_30K), log2(mean_50K/mean_10K), log2(mean_30K/mean_10K)), na.rm = TRUE))
  
  anova_summary[which(is.infinite(anova_summary$L2FC)), "L2FC"] <- 10
  return(anova_summary)
}

anova_summary <- anova_prot(data_means, data_sds, data_reps, data_norm, fractions)
anova_summary <- anova_summary[-which(is.na(anova_summary$q)), ]

keyvals <- ifelse((anova_summary$q < 0.05 & abs(anova_summary$L2FC) > 1) , '#0057b7','grey')
names(keyvals)[keyvals == '#0057b7'] <- 'significant'
names(keyvals)[keyvals == 'grey'] <- 'non-significant'

EnhancedVolcano(anova_summary, x = "L2FC", y = "q", lab = "", pCutoff = 0.05, FCcutoff = 1, colCustom = keyvals,
                ylab = bquote(~-log[10]~italic(q)),
                ylim = c(0,16), xlim = c(-12, 12), 
                title = "", subtitle = "", pointSize = 0.5, gridlines.major = FALSE, gridlines.minor = FALSE)
ggsave("./Plots/Volcano_Proteins_anova_between_fractions.pdf",  device = "pdf", width = 12 , height = 16, units = "cm")