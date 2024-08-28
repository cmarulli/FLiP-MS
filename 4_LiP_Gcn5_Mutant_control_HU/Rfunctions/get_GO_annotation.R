library(httr)
library(CCprofiler)
get_GO_annotation <- function(protein_list){
  
  GO_identifiers <- character()
  GO_names <- character()
  proteins <- character()
  IPRs <- character()
  Start <- numeric()
  End <- numeric()
  type <- character()
  
  for (p in 1:length(protein_list) ){
    
    current_prot <- protein_list[p]
    
    querry <- paste("https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/", current_prot,
                    sep = "")
    interPro_data <- GET(querry)
    interPro_content <- content(interPro_data)
    
    interPro_results <- interPro_content$results
    
    for (domain in 1:length(interPro_results) ){
      protein <- interPro_results[[domain]]$proteins[[1]]$accession
      IPR <- interPro_results[[domain]]$metadata$accession
      
      if (length(interPro_results[[domain]]$metadata$go_terms) > 0){
        for (go_counter in 1:length(interPro_results[[domain]]$metadata$go_terms)){
          
          i <- 1
          j <- 1
          
          for (i in 1:length(interPro_results[[domain]]$proteins[[1]]$entry_protein_locations)){
            for (j in 1:length(interPro_results[[domain]]$proteins[[1]]$entry_protein_locations[[i]]$fragments)){
              
              
              
              Start <- rbind(Start, interPro_results[[domain]]$proteins[[1]]$entry_protein_locations[[i]]$fragments[[j]]$start)
              End <- rbind(End, interPro_results[[domain]]$proteins[[1]]$entry_protein_locations[[i]]$fragments[[j]]$end)
              
              GO_identifiers <- rbind(GO_identifiers, interPro_results[[domain]]$metadata$go_terms[[go_counter]]$identifier)
              GO_names <- rbind(GO_names, interPro_results[[domain]]$metadata$go_terms[[go_counter]]$name)
              
              if (length(interPro_results[[domain]]$metadata$go_terms[[go_counter]]$name) > 0){
                proteins <- rbind(proteins, protein)
                IPRs <- rbind(IPRs, IPR)
                type <- rbind(type, interPro_results[[domain]]$metadata$type)
              } # if there is a GO term
            } # for all entry locations
          } # for all fragments of the domain
        } # for all GO terms
      } # if there is GO metadata
    } # for all domains
  } # for all proteins
  
  output <- data.frame(proteins, IPRs, GO_names, GO_identifiers, Start, End, type, row.names = NULL, stringsAsFactors = FALSE)
  output$proteins <- toupper(output$proteins)
  return(output)
}
