run_pageRank <- function(dampening_factor = 0.85, pageRank_percentile_cutoff = 0.75, FiLiP_changes, filename){
  
  
  ##### load data
  load( file = "intermediate/full_network.rds")
  load( file = "intermediate/layout_full.rds")

  
  #### run page rank
  seeds <- rep(0, length(V(network_full)$name))
  index_seed <- which(V(network_full)$name %in% FiLiP_changes)
  seeds[index_seed] <- 1 
  
  pageRank_full <- page.rank(network_full, algo = "prpack",  personalized = seeds,
            directed = FALSE, damping = dampening_factor)
  
  boxplot_data <- data.frame(cbind(pageRank_full[[1]], rep("All Proteins", length(V(network_full)$name))))
  pageRank_pbody <- pageRank_full[[1]][names(pageRank_full[[1]]) %in% FiLiP_changes] 
  temp <- data.frame( cbind(pageRank_pbody, rep("FiLiP Changes Proteins", length(pageRank_pbody))) )
  colnames(temp) <- colnames(boxplot_data)
  
  boxplot_data <- rbind(boxplot_data, temp)
  boxplot_data[, 1] <- as.numeric(boxplot_data[, 1])
  
  p <- ggplot(boxplot_data, aes(x = X2, y = X1, fill = X2)) +            # Applying ggplot function
        geom_boxplot() +
        scale_fill_manual(values = c("grey75", "#0057b7")) +
        labs(y= "PageRank score") +
        ggtitle("Personalized pageRank score distribution") +
        theme_classic() +
        theme(axis.title.x=element_blank(), legend.title = element_blank())
  ggsave(paste("plots/pageRank_scores_", filename, ".pdf", sep = ''), p, width = 21, height = 29.7, units = "cm")
  
  
  #### plot full network colored by pageRank score
  V(network_full)$pageRank <-  pageRank_full[[1]]
  fine = 500 # this will adjust the resolving power.
  pal = colorRampPalette(c("white", "#0057b7","#001e64"))
  
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  graphCol  <- pal(fine)[as.numeric(cut(V(network_full)$pageRank, breaks = fine))]
  
  V(network_full)$size <- 1200
  V(network_full)$color <- graphCol
  V(network_full)$frame.color <- "grey75"
  V(network_full)[which(V(network_full)$name %in% FiLiP_changes)]$frame.color <- "#ffc900"
  
  pdf(paste("plots/full_network_colored_by_pageRank_", filename, ".pdf", sep = ''), pointsize = 2)
    plot(network_full, layout = layout_full, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)   # need to adjust limits with new data
    title("Full input network - colored by pageRank score", cex.main = 2)
  dev.off()
  
  
  #### get propagated network
  pageRank_score_cutoff <- quantile(pageRank_full[[1]], pageRank_percentile_cutoff)
  index_pagerank <- which(pageRank_full[[1]] >= pageRank_score_cutoff)
  
  propagted_proteins <- names(pageRank_full[[1]])[index_pagerank] 
  index_pagerank <- which(V(network_full)$name %in% propagted_proteins)
  
  network_propagated <- induced_subgraph(network_full, index_pagerank)
  layout_porpagated <- layout_full[index_pagerank, ]
  
  #### plot propagated network
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  pdf(paste("plots/pageRank_network_", filename, ".pdf", sep = ''),  pointsize = 2)
    plot(network_propagated, layout = layout_porpagated, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)   
    title("Network after pageRank", cex.main = 2)
  dev.off()
  
  #### export results
  save(network_propagated, file = paste("intermediate/pageRank_network_", filename, ".rds", sep = ''))
  save(layout_porpagated, file = paste("intermediate/layout_pageRank_", filename, ".rds", sep = ''))
  write.table(propagted_proteins, paste("intermediate/propagted_proteins_", filename, ".tsv", sep =''), quote = F , sep = "\t", row.names = F)
}