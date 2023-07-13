walktrap_clustering <- function(step_number = 6){
 
   ##### load data
  load(file = "intermediate/pageRank_network_M.rds")
  network_M <- network_propagated
  load(file = "intermediate/pageRank_network_WT.rds")
  network_WT <- network_propagated
  
  load(file = "intermediate/full_network.rds")
  load(file = "intermediate/layout_full.rds")
  
  propagted_proteins <- unique(c(V(network_M)$name, V(network_WT)$name)) 
  index_propagated <- which(V(network_full)$name %in% propagted_proteins)
  network_propagated <- induced_subgraph(network_full, index_propagated)
  layout_propagated <- layout_full[index_propagated, ]
  
  ##### run walktrap clustering
  walktrap_clusters <- cluster_walktrap(network_propagated, steps = step_number,
                                        merges = TRUE, modularity = TRUE, membership = TRUE)
  
  
  V(network_propagated)$Cluster_membership <- membership(walktrap_clusters)
  
  pdf("plots/Walktrap_cluster_size.pdf", pointsize = 2)
    hist(table(membership(walktrap_clusters)), breaks = 50)
  dev.off()
  
  ##### export results
  cluster_export_table <- as.matrix(membership(walktrap_clusters))
  write.table(cluster_export_table, "intermediate/walktrap_clusters.tsv", quote = F , sep = "\t", row.names = T)
  
  save(network_propagated, file = "intermediate/pageRank_network_merged.rds")
  save(layout_propagated, file = "intermediate/layout_pageRank.rds")
  
  
}