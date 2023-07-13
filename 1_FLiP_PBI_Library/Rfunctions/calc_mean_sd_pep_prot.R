# function that calculates the mean peptide per protein value plus the corresponding standard deviation with error propagation


calc_mean_sd_pep_prot <- function(peptide_protein, lip_sub, tc_sub, n, column_names, filename){
  
  
  
  # store the mean peptides/protein intensities and the corresponding propagated error in mean_sd for each fraction
  mean_peptides_proteins <- data.frame(row.names = lip_sub$Peptide)
  sd_peptides_proteins <- data.frame(row.names = lip_sub$Peptide)
  nb_replicates <- data.frame(row.names = lip_sub$Peptide)
  
  # loop over each fraction/condition
  for (i in c(0:(n-1))){
    
    # get the peptides and proteins measured in one fraction
    pep_prot <- peptide_protein[, c((i*4 +1) : (i*4 + 4), ncol(lip_sub) -1 , ncol(peptide_protein))]
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
      current_pep_prot <- pep_prot[j,]
      index_replicates <- which(current_peptide[,c(1:4)]!=0)
      # store in how many replicates a peptide has been measured
      replicates <- c(replicates,  length(index_replicates))
      
      # if the peptide was measured at least in triplicates compute the mean and sd
      if(length(index_replicates) >= 3){
        
        mean_peptides <- mean(as.numeric(current_peptide[, index_replicates]))
        mean_proteins <- mean(as.numeric(proteins[current_peptide$Protein, index_replicates ]))
        
        # mean_pep_prot_current <- mean_peptides / mean_proteins
        mean_pep_prot_current <- mean(as.numeric(current_pep_prot[, index_replicates]))
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
  colnames(mean_peptides_proteins) <- column_names
  colnames(sd_peptides_proteins) <- column_names
  colnames(nb_replicates) <- column_names
  
  # replace 0 with NA for violin plots
  violin_df_mean_peptide_protein <- mean_peptides_proteins
  violin_df_mean_peptide_protein[violin_df_mean_peptide_protein == 0] <- NA 
  
  violin_df_sd_peptides_proteins <- sd_peptides_proteins
  violin_df_sd_peptides_proteins[violin_df_sd_peptides_proteins == 0] <- NA 
  
  violin_plot(violin_df_mean_peptide_protein, my_title = "Mean Peptide/Protein", my_ylab = "log10(mean(pep/prot))", my_filename = paste("mean_peptide_proteins_", filename, sep =""))
  violin_plot(violin_df_sd_peptides_proteins, my_title = "Standard Deviation Peptide/Protein", my_ylab = "log10(sd(pep/prot))", my_filename = paste("sd_peptide_proteins_", filename, sep =""))
  
  # keep only peptides that have been measured in at least 2 fractions --> hopefully once bound to complex and once as monomer
  mean_peptides_proteins[is.na(mean_peptides_proteins)] <- 0
  sd_peptides_proteins[is.na(sd_peptides_proteins)] <- 0
  
  index <- apply(mean_peptides_proteins, 1, function(x) length(which(x != 0 & !is.infinite(x) & !is.na(x))))
  mean_peptides_proteins <- mean_peptides_proteins[which(index >= 2), ]
  sd_peptides_proteins <- sd_peptides_proteins[which(index >= 2), ]
  nb_replicates <- nb_replicates[which(index >= 2), ] 
  
  return(list(mean_peptides_proteins, sd_peptides_proteins, nb_replicates))
  
}
