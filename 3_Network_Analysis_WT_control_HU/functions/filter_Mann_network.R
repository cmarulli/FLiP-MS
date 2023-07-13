filter_Mann_network <- function(mann_input_file, mapping_table_file,
                                mann_output_filename){


  ##### load data
  mann_input <- fread(mann_input_file)
  mapping_table <- fread(mapping_table_file)
  
  
  ##### remove interactions of non-unique nodes
  index <-  which(grepl( ";", mann_input$source, fixed = TRUE) | grepl( ";", mann_input$target, fixed = TRUE)) # filter non-unique entries
  mann_input <- mann_input[-index]
  
  ##### create mapping table for uniport IDs
  mapping_table <- data.table(mapping_table$name, mapping_table$`UniProt ID`)
  index <-  which(grepl( ";", mapping_table$V1, fixed = TRUE) | grepl( ";", mapping_table$V2, fixed = TRUE)) # filter non-unique entries
  mapping_table <- mapping_table[-index]
  
  lookup <- setNames(mapping_table$V2, mapping_table$V1)
  
  #map
  mann_input <- transform(mann_input,  protein1_ID=lookup[ source], stringsAsFactors=FALSE)
  mann_input <- transform(mann_input,  protein2_ID=lookup[ target], stringsAsFactors=FALSE)
  
  
  #### filter for self-loops
  index <- which(mann_input$protein1_ID == mann_input$protein2_ID)
  if (length(index) >0){
    mann_input <- mann_input[-index,]  
  }
  
  
  
  
  
  #### filter for AB = BA
  AB <- paste(mann_input$protein1_ID, "_", mann_input$protein2_ID, sep = "")
  BA <- paste(mann_input$protein2_ID, "_", mann_input$protein1_ID, sep = "")
  
  index <- which(AB == BA)
  if (length(index) >0){
    mann_input <- mann_input[-index,]  
  }
  
  #### reduce to columns of interest
  mann_input <- mann_input[, c("protein1_ID", "protein2_ID")]
  
  
  #### export
  write.table(mann_input, paste("intermediate/", mann_output_filename, ".tsv", sep = "" ), quote = F , sep = "\t", row.names = F)
}
