cluster_GO_analysis <- function(enrichment_q_cutoff = 0.05,
                                detected_proteins){
  

  
  
  
  ##### load data
  load(file = "final/clustered_network_final.rds")
  load( file = "final/layout_clustered_new.rds")
  load(file = "intermediate/pageRank_network_merged.rds")
  
  # proteins from network propagation
  network_proteins <- V(network_propagated)$name
  
  # merge
  background <- unique(c(detected_proteins, network_proteins))
  
  # annotate with uniprot info
  uniprot <- fetch_uniprot(background)
  
  
  #### define hit list
  is_Cluster_member <- rep(FALSE, length(background))
  
  membership <- unique(V(network_clustered)$Cluster_membership)
  
  for (i in 1:length(membership)){
    
    current_cluster <- membership[i]
    data_cluster <- data.frame("proteins" = background,
                               "cluster_membership" = is_Cluster_member)
    
    index_cluster <-  which(data_cluster$proteins %in%  V(network_clustered)$name[ which(V(network_clustered)$Cluster_membership == current_cluster) ])
    
    data_cluster$cluster_membership[index_cluster] <- TRUE
    
    data_cluster <- left_join(data_cluster, uniprot, by = c("proteins" = "id"))
    
    # GO analysis
    
    MF <- go_enrichment(data_cluster,
                          protein_id = proteins,
                          is_significant = cluster_membership,
                          go_annotations_uniprot = go_molecular_function,
                          plot = F)

    CC <- go_enrichment(data_cluster,
                          protein_id = proteins,
                          is_significant = cluster_membership,
                          go_annotations_uniprot = go_cellular_compartment,
                          plot = F)
    
    BP <- go_enrichment(data_cluster,
                  protein_id = proteins,
                  is_significant = cluster_membership,
                  go_annotations_uniprot = go_biological_process,
                  plot = F)
    
    EA_data <- data.frame( "term" = c(MF$term, CC$term, BP$term),
                           "q" = c(MF$adj_pval, CC$adj_pval, BP$adj_pval),
                           "number_of_proteins_in_cluster" = c(MF$n_significant_proteins_in_process, CC$n_significant_proteins_in_process, BP$n_significant_proteins_in_process),
                           "type" = c(rep("MF", length(MF$term)), rep("CC", length(CC$term)), rep("BP", length(BP$term))) )
    
    # limit number of letters of the terms (for plotting)
    EA_data$term <- sub("(?<=^.{34}).*","...", EA_data$term, perl = T)


    # EA_data <- data.frame( "term" = c(CC$term, BP$term),
    #                        "q" = c(CC$adj_pval, BP$adj_pval),
    #                        "number_of_proteins_in_cluster" = c(CC$n_significant_proteins_in_process, BP$n_significant_proteins_in_process),
    #                        "type" = c( rep("CC", length(CC$term)), rep("BP", length(BP$term))) )


    
    EA_data <- EA_data[EA_data$q <= enrichment_q_cutoff, ]
    
    
    
    
    # barplot
    EA_data$log_q <- -log10(EA_data$q)
    EA_data <- EA_data[order(EA_data$log_q, decreasing = FALSE), ]
    EA_data <- EA_data[-which(duplicated(EA_data$term)), ]
    
    EA_data$term <- factor(EA_data$term, levels = EA_data$term) # check if unique term: if not, keep smaller q value
    
    p <- ggplot(EA_data, aes(x = term, y = log_q, fill = type)) +
                      geom_bar(stat="identity") +
                      scale_fill_manual(values = c("#b0c7e0","#0057b7","#001e64")) +
                      theme_classic() +
                      ylab( "-log10(adj. p-value)") +
                      ggtitle(paste("GO enrichment cluster", as.character(current_cluster))) +
                      theme_classic(base_size = 22) +
                      theme( axis.title.y=element_blank(), legend.title = element_blank(),
                             axis.text.y = element_text(size = 14), legend.text = element_text(size = 14), 
                             axis.text.x = element_text(size = 14), plot.title = element_text(size = 14),
                             axis.title.x = element_text(size = 14)) +
                      coord_flip()
    ggsave(paste("plots/GO_Cluster_", as.character(current_cluster), ".pdf", sep = ""), p, width = 29.7, height = 21, units = "cm")
    
    
    # corresponding network
    upper <-  apply(layout_clustered_new, 2, max)
    lower <-  apply(layout_clustered_new, 2, min)
    
    network_current <- network_clustered
    index_current <- which(V(network_current)$Cluster_membership != current_cluster)
    V(network_current)$color[index_current] <- "grey75"
    V(network_current)$frame.color[index_current] <- "grey75"
    
    pdf(paste("plots/GO_Cluster_", as.character(current_cluster), "_network.pdf", sep = ""), pointsize = 2)
      plot(network_current, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0)  
      title(paste("Cluster", as.character(current_cluster)), cex.main = 2)
    dev.off()
    
    
    # save enrichment analysis data
    write.table(EA_data, paste("final/GO_Cluster_", as.character(current_cluster), ".tsv", sep = ""),
                quote = F , sep = "\t", row.names = F)
  
   }
}
