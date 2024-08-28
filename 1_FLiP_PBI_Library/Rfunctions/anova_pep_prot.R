anova_pep_prot <- function(mean_peptides_proteins, sd_peptides_proteins, nb_replicates, lip_data, conditions){
  
  # iterate over each peptide and calculate the F statistic
  F_values <- c()
  p_values <- c()
  MS_errors <- c()
  MS_groups <- c()
  
  # ANOVA: calculate F value and corresponding probability 
  for (i in c(1:nrow(mean_peptides_proteins))){
    
    # extract the means and the standard deviations of the peptide/protein ratio in all fractions for the current peptide
    means <- as.vector(as.numeric(mean_peptides_proteins[i, ]))
    sds <- as.vector(as.numeric(sd_peptides_proteins[i, ]))
    rep <- as.vector(as.numeric(nb_replicates[i, ]))
    
    # if index != 0 then peptide was not measured in this fraction
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
  anova_summary <- cbind(rownames(mean_peptides_proteins), F_values, p_values, q_values, MS_errors, MS_groups, mean_peptides_proteins)
  colnames(anova_summary) <- c("Peptide", "F", "p", "q", "MSError", "MSGroup", conditions)
  anova_summary <- merge(anova_summary, lip_data[, c("Protein","Peptide", "tryptic")])
  
  return(anova_summary)
  
}
