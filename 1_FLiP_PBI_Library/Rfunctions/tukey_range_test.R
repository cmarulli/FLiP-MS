

tukey_range_test <- function(anova_summary, n, nb_replicates, my_conditions){

  p_tukey <- data.frame()
  L2FC_anova <- data.frame()
  
  for (i in c(1:nrow(anova_summary))){
    means <- anova_summary[i, grep('Mean', colnames(anova_summary))]
    means[which(means == 0)] <- NA
    peptide <- anova_summary[i, 'Peptide']
    
    p_values <- c()
    t_values <- c()
    L2FC_values <- c()
    
    for (a in c(1:(n-1))){
      for (b in c((a+1):n)){
        
        # caluclate the q value
        t <- as.numeric((means[a] - means[b])/sqrt(anova_summary[i,5]/4))
        # calculate the log2 FC
        L2FC<- log2(as.numeric(means[a])/ as.numeric(means[b]))
        # get the probablity to observe a q-value as extream as the calculated one
        # nmeans: number of fractions = 4,  degrees of freedom (N)
        # number of replicates: as given in nb_replicates (R)
        # degree of freedom = R-N
        # lower.tail = FALSE --> probability to observe higher value
        # use the absolute value as some means are negative
        nb_rep <- sum(nb_replicates[peptide, ])
        
        p <- c(ptukey(abs(t), nmeans = 4, df = (nb_rep - n), lower.tail = FALSE))
        
        # store p and q
        p_values <- c(p_values, p)
        L2FC_values <- c(L2FC_values, L2FC)
      }
    }
    
    p_values <- p.adjust(p_values, method = "BH")
    
    # store p and q for all combinations
    p_tukey <- rbind(p_tukey, p_values)
    L2FC_anova <- rbind(L2FC_anova, L2FC_values)
  }
  
  combinations <- get_combinations(my_conditions) 
  
  # adjust column and row names of Tukey p and q values

  colnames(p_tukey) <- paste("q(", combinations, ")", sep = '')
  colnames(L2FC_anova) <- combinations

  rownames(p_tukey) <- anova_summary$Peptide
  rownames(L2FC_anova) <- anova_summary$Peptide
  
  
  
  # extract the most significant L2FC observed for each peptide
  maxL2FC <- c()
  most_sign_change <- c()
  for (i in c(1: nrow(L2FC_anova))){
    
    L2FCs <- as.numeric(L2FC_anova[i, ])
    p <- as.numeric(p_tukey[i, ])
    # set all the comparison to 0 intensities to 0
    L2FCs[which(is.infinite(L2FCs))] <- 0
    # extract the maximum L2FC for each peptide
    # maxL2FC <- c(maxL2FC, L2FCs[which(abs (L2FCs) == max(abs(L2FCs), na.rm = TRUE))[1]])
    # extract the most significant L2FC for each peptide
    maxL2FC <- c(maxL2FC, L2FCs[which(p == min(p, na.rm = TRUE))[1]])
    most_sign_change <- c(most_sign_change, combinations[which(p == min(p, na.rm = TRUE))[1]])
  }
  
  anova_summary$L2FC <- maxL2FC
  anova_summary$most_sig_change <- most_sign_change
  anova_summary <- cbind(anova_summary, p_tukey)
  
  # what is the most significiant difference in peptides intensity means observed and is the trend
  # difference_mean <- c()
  # effect_size <- c()
  # for (i in c(1:nrow(p_tukey))){
  #   
  #   # which comparison of means is significiantly different
  #   index <- which(as.numeric(p_tukey[i,]) < 0.05)
  #   
  #   if(length(index) == 0){
  #     effect_size <- c(effect_size, 0)
  #     difference_mean <- c(difference_mean, 'not_sign')
  #     
  #   }else{
  #     
  #     p_values <- p_tukey[i, index]
  #     max <- max(p_values)
  #     min <- min(p_values)
  #     
  #     if (abs(max) > abs(min)){
  #       effect_size <- c(effect_size, max)
  #     }else{
  #       effect_size <- c(effect_size, min)
  #     }
  #     
  #     # if t he q value is bigger than 0 assign TRUE, else assign FALSE
  #     p_values <- p_values > 0
  #     
  #     # is the trend of sinificiantly different means the same? e.g. all bigger fractions higher proteins intensities
  #     # if so, assign the label increase, decrease
  #     
  #     if(length(unique(as.vector(p_values))) == 1){
  #       if (p_values[1] == TRUE){
  #         difference_mean <- c(difference_mean, 'increased')
  #       }else{
  #         difference_mean <- c(difference_mean, 'decreased')
  #       }
  #     }else{
  #       difference_mean <- c(difference_mean, 'unclear')
  #     }
  #   }
  # }
  # 
  # anova_summary$effect_size <- effect_size
  # anova_summary$difference_mean <- difference_mean
  # 
  # # we expect semi-tryptic peptides to decrease in intensity in the complex form compared to the monomer form fully-tryptic peptides to increase
  # expected_at_interface <- c()
  # for (i in c(1:nrow(anova_summary))){
  #   
  #   tryp <- anova_summary[i, 'tryptic']
  #   diff <- anova_summary[i, 'difference_mean']
  #   
  #   if((tryp == "Specific" & diff == "increased") | (tryp != "Specific" & diff == "decreased")){
  #     expected_at_interface <- c(expected_at_interface, TRUE)
  #   } else if(diff == "unclear"){
  #     expected_at_interface <- c(expected_at_interface, TRUE)
  #   } else{
  #     expected_at_interface <- c(expected_at_interface, FALSE)
  #   }
  # }
  # 
  # anova_summary$expected_at_interface <- expected_at_interface
  
  return(list(anova_summary, p_tukey, L2FC_anova))
}
