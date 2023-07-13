markov_clustering <- function(markov_expansion_parameter = 2, markov_inflation_coefficient = 2){

  ##### load data
  load(file = "intermediate/pageRank_network.rds")
  load(file = "intermediate/layout_pageRank.rds")
  
  ##### get adjacency matrix
  adj_matrix <- as_adjacency_matrix(network_porpagated, type = "both", edges = FALSE)
  
  ##### run markov clustering
  adj_matrix <- as.matrix(adj_matrix)
  
  markov_clusters <- mcl(adj_matrix, addLoops = TRUE, expansion = markov_expansion_parameter, inflation = markov_inflation_coefficient,
                          allow1 = TRUE, max.iter = 100, ESM = FALSE)
  
  
  V(network_porpagated)$Cluster_membership <- markov_clusters$Cluster
  
  png("plots/markov_cluster_size.png", width=720,height = 720)
    hist(table(markov_clusters$Cluster), breaks = 50)
  dev.off()
  
  ##### export results
  cluster_export_table <- data.frame(V(network_porpagated)$name, V(network_porpagated)$Cluster_membership)
  write.table(cluster_export_table, "intermediate/markov_clusters.tsv", quote = F , sep = "\t", row.names = T)
  
}