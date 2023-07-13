triplicates_tc <- function(data, n, conditions, filename){
  
  
  # change type of intensity values to numerical
  for (i in c(2:(ncol(data)))){
    data[,i] <- as.numeric(data[,i])
  }
  
  #replcae NAs with 0
  data[is.na(data)] <- 0
  
  # create a list where each entry corresponds to the replicates of a single condition
  proteins_triplicates_list <- list()
  for(i in c(1:n)){
    proteins_triplicates_list[[i]] <- data[, c(1,2, ((i-1)*4+3):((i-1)*4+6))]
  }
  
  # filter the proteins that have not been identified in at least three replicates and plot how many proteins were identified in each condition
  nb_identified_proteins <- c()
  identified_proteins <- list()
  
  # assign intensity 0 to all peptides that have not been measured in at least triplicates
  for(i in 1:length(proteins_triplicates_list)){
    
    temp <- proteins_triplicates_list[[i]]
    index <- apply(temp[, 3:6], 1, function(x) length(which(x == 0)))
    temp[which(index > 1), c(3:6)] <- 0 
    
    proteins_triplicates_list[[i]] <- temp
    nb_identified_proteins <- c(nb_identified_proteins, length(which(index <= 0)))
    identified_proteins[[i]] <-  temp[which(index <= 1), 1]
  }
  
  names(identified_proteins) <- conditions
  
  # save the complete filtered data set in the data variable
  data <- proteins_triplicates_list[[1]]
  for(i in c(2:length(proteins_triplicates_list))){
    data <- merge(data, proteins_triplicates_list[[i]])
  }
  
  # get rid of the proteins that have not been measured in at least triplicates in any of the fractions
  index <- apply(data[, c(3:(2+n*4))], 1, function(x) length(which(x == 0)))
  if(length(which(index == n*4)) > 0){
    data <- data[-which(index == n*4), ]
  }
  
  return(list(data, nb_identified_proteins, identified_proteins))
  
}

