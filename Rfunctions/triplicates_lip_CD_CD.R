# function that returns dataframe of peptides/proteins that were measured in at least triplicates
# also returns the number of peptides/proteins per fraction
# distinct functions for lip/tc data
# requires function violin_plot

library("wesanderson")
pal <- wes_palette("Darjeeling1", 40, type = "continuous")


triplicates_lip <- function(data, n, conditions, filename){
  
  #exclude iRT peptides
  if(length(grep("iRT", data[,2])) > 0){
    data <- data[-grep("iRT", data[,2]), ]
  }
  
  # change type of intensity values to numerical
  for (i in c(3:(2+4*n))){
    data[,i] <- as.numeric(data[,i])
  }
  
  violin_plot(data[, c(3:(2+n*4))], my_title = "Peptide Intensities LiP", my_filename = paste(filename, "_Peptides", sep =""))
  
  #replcae NAs with 0
  data[is.na(data)] <- 0
  
  # safe metadata
  meta_data <- data[, - c(3:(2+n*4))]
  
  # create a list that contains all the dataframes
  peptides_triplicates_list <- list()
  for(i in c(1:n)){
    peptides_triplicates_list[[i]] <- data[, c(1,2, ((i-1)*4+3):((i-1)*4+6))]
  }
  
  # filter the peptides that have not been identified in all three replicates and plot how many peptides were identified in each condition
  nb_identified_peptides <- c()
  identified_peptides <- list()
  
  # assign intensity 0 to all peptides that have not been measured in at least triplicates
  for(i in 1:length(peptides_triplicates_list)){
    
    temp <- peptides_triplicates_list[[i]]
    index <- apply(temp[, c(3:6)], 1, function(x) length(which(x == 0)))
    temp[which(index > 1), c(3:6)] <- 0 
    
    # sum up peptide intensities of the same acetylated peptide with different charged states that have been identified in at least 3 runs
    lip_sub <- temp[,-2]
    lip_sub <- lip_sub %>% 
      group_by(mod_Peptide) %>% 
      summarise_all(funs(sum))
    
    temp <- merge(lip_sub, temp[, c(1,2)], by = "mod_Peptide" )
    temp <- unique(temp)
    temp <- temp[, c(1,6,2:5)]
    
    # add the coefficient of variation for each peptide
    CVs <- apply(temp[, c(3:6)], 1, function(x) sd(x) / mean(x))
    temp <- cbind(temp, CVs)
    colnames(temp) <- c(colnames(temp)[1:6], paste("CV_", colnames(temp)[3], sep = ""))
    
    peptides_triplicates_list[[i]] <- temp
    nb_identified_peptides <- c(nb_identified_peptides, length(which(index <= 0)))
    identified_peptides[[i]] <- temp[which(index <= 1), 1]
  }
  
  names(identified_peptides) <- conditions
  
  # save the complete filtered data set in the data variable
  data <- peptides_triplicates_list[[1]]
  for(i in c(2:length(peptides_triplicates_list))){
    data <- merge(data, peptides_triplicates_list[[i]])
  }
  
  data <- data[, c("mod_Peptide", "Protein", 
                   as.vector(sapply(conditions, function(x) return(paste(x, c(1:4), sep = "_")))), 
                   as.vector(sapply(conditions, function(x) return(paste("CV", x, "1", sep = "_")))))]

  # get rid of the peptides that have not been measured in at least triplicates in any of the fractions
  index <- apply(data[, c(3:(2+n*4))], 1, function(x) length(which(x == 0)))
  if (length(which(index == n*4) != 0)){
    data <- data[-which(index == n*4), ]
  }

  data <- merge(data, meta_data[-4])
  data <- unique(data)
  
  return(list(data, nb_identified_peptides, identified_peptides))
}




