betweenness_clustering <- function(){
 
   ##### load data
  load(file = "intermediate/pageRank_network.rds")
  load(file = "intermediate/layout_pageRank.rds")
  
  ##### run clustering
  betweenness_clusters <- cluster_edge_betweenness(network_porpagated, directed = FALSE, 
                                                merges = TRUE, modularity = TRUE, membership = TRUE)
  
  
  
  
  V(network_porpagated)$Cluster_membership <- membership(betweenness_clusters)
  
  png("plots/Walktrap_cluster_size.png", width=720,height = 720)
    hist(table(membership(betweenness_clusters)), breaks = 50)
  dev.off()
  
  ##### export results
  cluster_export_table <- as.matrix(membership(betweenness_clusters))
  write.table(cluster_export_table, "intermediate/betweenness_clusters.tsv", quote = F , sep = "\t", row.names = T)
}