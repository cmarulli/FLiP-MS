#function that annotates peptides as bein fully or semi-tryptic based on the output of the python script 20190322_find_peptide_in_fasta.ipynb
fully_semi_tryptic <- function(sequence_info){
  
  N_term_tryp <- which(sequence_info$AA_before == "-" | 
                         sequence_info$AA_before == "K" |
                         sequence_info$AA_before == "R" )
  
  # exclude peptides where a Proline comes prior to K or R beacuse it will not be recognized by trypsin
  N_term_P <- which(sequence_info$First_AA == "P" )
  N_term_tryp <- setdiff(N_term_tryp, N_term_P)
  
  sequence_info$trpN <- 0
  sequence_info$trpN[N_term_tryp] <- 1
  
  
  C_term_tryp <- which(sequence_info$Last_AA == "K" |
                         sequence_info$Last_AA == "R" |
                         sequence_info$AA_after == "-")
  
  C_term_P <-  which(sequence_info$AA_after == "P")
  C_term_tryp <- setdiff(C_term_tryp, C_term_P)
  
  sequence_info$trpC  <- 0
  sequence_info$trpC[C_term_tryp] <- 1
  
  sequence_info$tryptic <- apply(sequence_info[, c(ncol(sequence_info)-1, ncol(sequence_info))], 1, function(x) {
    
    if(x[1] == 1 & x[2] == 1){
      return("fully")
    }else if(x[1] == 0 & x[2] == 0){
      return("not")
    }else{
      return("semi")
    }
    
  })
  
  sequence_info <- unique(sequence_info)
  return(sequence_info)
}
