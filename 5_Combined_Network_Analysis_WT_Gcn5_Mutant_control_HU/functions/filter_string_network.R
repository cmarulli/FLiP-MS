filter_string_network <- function(score_to_use = "combined_score",  score_cutoff = 700,
                                  string_input_file, mapping_table_file,
                                  string_output_filename){
  
  
  ##### load data
  string_input <- fread(string_input_file)
  mapping_table <- fread(mapping_table_file)
  
  ##### create mapping table for uniport IDs
  mapping_table$name <- paste("4932.", mapping_table$`Gene names  (ordered locus )`, sep ="")
  lookup <- setNames(mapping_table$Entry, mapping_table$name)
  
  #map
  string_input <- transform(string_input,  protein1_ID=lookup[ protein1], stringsAsFactors=FALSE)
  string_input <- transform(string_input,  protein2_ID=lookup[ protein2], stringsAsFactors=FALSE)
  
  
  ##### filter based on score
  index <- which(string_input[, eval(score_to_use)] >= score_cutoff)
  string_output <- string_input[index, ]
  
  
  #### filter NAs in mapping
  index  <- which(is.na(string_output$protein1_ID))
  string_output <- string_output[-index,] 
  index  <- which(is.na(string_output$protein2_ID))
  string_output <- string_output[-index,] 
  
  
  #### filter for self-loops
  index <- which(string_output$protein1_ID == string_output$protein2_ID)
  if (length(index) >0){
    string_output <- string_output[-index,]  
  }
  
  
  #### filter for AB = BA
  AB <- paste(string_output$protein1_ID, "_", string_output$protein2_ID, sep = "")
  BA <- paste(string_output$protein2_ID, "_", string_output$protein1_ID, sep = "")
  
  index <- which(AB == BA)
  if (length(index) >0){
    string_output <- string_output[-index,]  
  }
  
  
  #### reduce to columns of interest
  string_output <- string_output[, c("protein1_ID", "protein2_ID")]
  
  
  #### export
  write.table(string_output, paste("intermediate/", string_output_filename, ".tsv", sep = "" ), quote = F , sep = "\t", row.names = F)
  
}