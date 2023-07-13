walktrap_clustering <- function(step_number = 6){
 
   ##### load data
  load(file = "intermediate/pageRank_network.rds")
  load( file = "intermediate/layout_pageRank.rds")
  
  ##### run walktrap clustering
  walktrap_clusters <- cluster_walktrap(network_porpagated, steps = step_number,
                                        merges = FALSE, modularity = TRUE, membership = TRUE)
  
  
  V(network_porpagated)$Cluster_membership <- membership(walktrap_clusters)
  
  pdf("plots/Walktrap_cluster_size.pdf", pointsize = 2)
    hist(table(membership(walktrap_clusters)), breaks = 50)
  dev.off()
  
  ##### export results
  cluster_export_table <- as.matrix(membership(walktrap_clusters))
  write.table(cluster_export_table, "intermediate/walktrap_clusters.tsv", quote = F , sep = "\t", row.names = T)
}