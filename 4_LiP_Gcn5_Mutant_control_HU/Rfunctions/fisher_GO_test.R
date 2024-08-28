fisher_GO_test <- function(peptides_GO_df){
  
  GO_terms <- unique(peptides_GO_df$GO_names)
  p_values <- c()
  
  for (GO in GO_terms){
    
    peptides_GO_df$is_GO <- rep(FALSE, nrow(peptides_GO_df))
    peptides_GO_df[which(peptides_GO_df$GO_names == GO), ]$is_GO <- TRUE
    
    tab <- table(peptides_GO_df$is_GO, peptides_GO_df$sign)
    
    p_values <- c(p_values, fisher.test(tab, alternative = 'greater')$p.value) 
  }
  
  p_values <- p.adjust(p_values, method =  "BH")
  return(data.frame("GO_Terms" = GO_terms, "p_values" = p_values)) 
}
