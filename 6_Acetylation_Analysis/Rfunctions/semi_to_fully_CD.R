semi_to_fully <- function(data){
  
  library("Biostrings")

  fastaFile <- readAAStringSet("20190927_uniprot_yeast_iRT_PK_MutL.fasta")
  library("stringi")
  
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df_aa_sequences <- data.frame(seq_name, sequence)
  colnames(df_aa_sequences) <- c("Protein", "Sequence")
  df_aa_sequences$Protein <- sapply(df_aa_sequences[,1], function(x) return(strsplit(x, split = "\\|")[[1]][2]))
  
  data <- merge(data, df_aa_sequences)
  
  fully_tryptic_peptides <- c()
  
  for (i in c(1:nrow(data))){
    
    # get the positions of the peptide in the protein sequence
    positions <- stri_locate(data[i,]$Sequence, regex = data[i,]$Peptide) 
    start <- positions[1]
    stop <- positions[2]
    
    if(data[i,]$tryptic == "Specific-N" || data[i,]$tryptic == "Unspecific"){
      
      start = start - 1
      while(start > 1){
        
        if(( substr(data[i,]$Sequence,start,start) == "K" || substr(data[i,]$Sequence,start,start) == "R") && substr(data[i,]$Sequence,(start+1),(start+1)) != "P"){
          start <- start + 1
          break()
        }else{
          start = start -1
        }
      }
      
      
    }
    
    if(data[i,]$tryptic == "Specific-C"|| data[i,]$tryptic == "Unspecific"){
      
      stop = positions[2] + 1
      while(stop < nchar(data[i, ]$Sequence)){
        
        if((substr(data[i,]$Sequence,stop,stop) == "K" || substr(data[i,]$Sequence,stop,stop) == "R") && substr(data[i,]$Sequence,(stop+1),(stop+1)) != "P"){
          break()
        }else{
          stop = stop + 1
        }
      }
    }
    
    fully_tryptic_peptides <- c(fully_tryptic_peptides, substr(data[i,]$Sequence, start, stop))
  }
  
  data$fully_tryptic <- fully_tryptic_peptides
  return(data)
}