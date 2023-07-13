cluster_analysis2 <- function(cluster_method = "walktrap", cluster_min_size = 10, cluster_max_size = 500,
                             cluster_col_pal = "Dark2", cluster_q_value_cutoff = 0.01){

  ##### load data
  load(file = "intermediate/pageRank_network.rds")
  load( file = "intermediate/layout_pageRank.rds")
  load( file = "intermediate/layout_full.rds")
  load(file = "intermediate/p-body_proteins.rds")
  
  if(cluster_method == "walktrap"){
    cluster <- read.table("intermediate/walktrap_clusters.tsv", stringsAsFactors = FALSE)
    V(network_porpagated)$Cluster_membership <- cluster$V1 
    
  } else if (cluster_method == "markov"){
    cluster <- read.table("intermediate/markov_clusters.tsv", stringsAsFactors = FALSE)
    V(network_porpagated)$Cluster_membership <- cluster$V.network_porpagated..Cluster_membership 
    
  } else if (cluster_method == "betweenness"){
    cluster <- read.table("intermediate/betweenness_clusters.tsv", stringsAsFactors = FALSE)
    V(network_porpagated)$Cluster_membership <- cluster$V1
  }
  
  
  
  #### filter clusters based on size
  clustr_sizes <- table( V(network_porpagated)$Cluster_membership)
  index_cluster <- which(clustr_sizes >= cluster_min_size & clustr_sizes <= cluster_max_size)
  
  number_of_clusters <-length(index_cluster)
  
  clusters_to_use <- names(clustr_sizes[ index_cluster ])
  index_cluster <- which(V(network_porpagated)$Cluster_membership %in% clusters_to_use)
  
  network_clustered <- induced_subgraph(network_porpagated, index_cluster)
  layout_clustered <- layout_porpagated[index_cluster, ]
  
  
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  mycolors <- colorRampPalette(brewer.pal(8, cluster_col_pal))(number_of_clusters)
  
  membership <- unique(V(network_clustered)$Cluster_membership)
  lookup = setNames(mycolors, as.character(membership))
  
  df <- data.frame( as.character(V(network_clustered)$Cluster_membership) )
  colnames(df) <- "membership"
  df <- transform(df, membership=lookup[membership], stringsAsFactors=FALSE)
  
  
  V(network_clustered)$color <- df[,1]
  V(network_clustered)$size <- 450
  
  
  png("plots/clustered_network_size_filter.png", width=3600,height = 3600)
    plot(network_clustered, layout = layout_clustered, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
         col=mycolors)  
    title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 6)
  dev.off()
    
  
  #### change coordinates
  layout_clustered_new <- layout_with_kk(network_clustered,  )
  
  upper <- 1.1 * apply(layout_clustered_new, 2, max)
  lower <- 1.1 * apply(layout_clustered_new, 2, min)
  
  
  V(network_clustered)$size <- 30
  png("plots/clustered_network_size_filter_new_layout.png", width=3600,height = 3600)
    plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
         col=mycolors)  
    title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 6)
  dev.off()
  
  
  ##### Test clusters if they are enriched in our hits
  ### Fisher for enrichment of original hits
  fisher_p <- rep(0, length(membership))
  total_number <- length(V(network_clustered))
  
  V(network_clustered)$origial_hit <- rep(0, total_number)
  index_pody <- which(V(network_clustered)$name %in% p_body_proteins)
  V(network_clustered)$origial_hit[index_pody] <- 1
  
  
  for (i in 1:length(membership)){
    current_cluster <- membership[i]
    
    # define contingency table
    YY <- length(which( V(network_clustered)$origial_hit == 1 & V(network_clustered)$Cluster_membership == current_cluster ))   # in original hit list & in cluster
    NY <- length(which( V(network_clustered)$origial_hit == 0 & V(network_clustered)$Cluster_membership == current_cluster ))   # not in original hit list & in cluster
    YN <- length(which( V(network_clustered)$origial_hit == 1 & V(network_clustered)$Cluster_membership != current_cluster ))   # in original hit list & not in cluster
    NN <- length(which( V(network_clustered)$origial_hit == 0 & V(network_clustered)$Cluster_membership != current_cluster ))   # not in original hit list & not in cluster
    
    dat <- data.frame(
      "in_cluster" = c(YY, NY),
      "not_in_cluster" = c(YN, NN),
      row.names = c("hit_list", "not_in_hit_list"),
      stringsAsFactors = FALSE
    )
    
    # save p value
    fisher_p[i] <- fisher.test(dat)$p.value
  }
  
  # adjust for multiple testing
  fisher_q <- p.adjust(fisher_p, method = "BH") 
  
  
  ### Kolmogorov-Smirnov Test to see if higher pageRank score
  ks_p <- rep(0, length(membership))
  
  
  # i <- 9
  for (i in 1:length(membership)){
    current_cluster <- membership[i]
    scores_in_cluster <- V(network_clustered)$pageRank[ which(V(network_clustered)$Cluster_membership == current_cluster) ]
    scors_not_in_cluster <- V(network_clustered)$pageRank[ which(V(network_clustered)$Cluster_membership != current_cluster) ]
    
    ks_p[i] <- wilcox.test(scores_in_cluster, scors_not_in_cluster, alternative = "greater")$p.value
  }
  
  # adjust for multiple testing
  ks_q <- p.adjust(ks_p, method = "BH")
  
  
  cluster_q <- data.frame("cluster" = membership,
                          "Fisher_q" = fisher_q,
                          "wilcoxon_q" =  ks_q)
  
  
  
  ##### reduce to significant cluster and plot
  cluster_q_sig <- cluster_q[which(cluster_q$Fisher_q <= cluster_q_value_cutoff | cluster_q$wilcoxon_q <=  cluster_q_value_cutoff), ]
  index_sig <- which( V(network_clustered)$Cluster_membership %in% cluster_q_sig$cluster)
  
  network_sig <- induced_subgraph(network_clustered, index_sig)
  layout_sig <- layout_clustered_new[index_sig, ]
  
  
  
  png("plots/clustered_network_significant_new_layout.png", width=3600,height = 3600)
    plot(network_sig, layout = layout_sig, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
         col=mycolors)  
    title(paste(as.character(length(unique(cluster_q_sig$cluster))), " significant communities found by ",  cluster_method, " clustering"), cex.main = 6)
  dev.off()
  
  # new layout
  layout_clustered_sig_new <- layout_with_lgl(network_sig, )
  V(network_sig)$size <- 250
  
  
  upper <- apply(layout_clustered_sig_new, 2, max)
  lower <- apply(layout_clustered_sig_new, 2, min)
  
  
  png("plots/clustered_network_significant_final_layout.png", width=3600,height = 3600)
    plot(network_sig, layout = layout_clustered_sig_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
         col=mycolors)  
    title(paste(as.character(length(unique(cluster_q_sig$cluster))), " significant communities found by ",  cluster_method, " clustering"), cex.main = 6)
  dev.off()
  
  
  ##### export results
  write.table(cluster_q, "final/cluster_q_values.tsv", quote = F , sep = "\t", row.names = F)
  save(network_sig, file = "final/cluster_sig_network.rds")
  save(layout_clustered_sig_new, file = "final/layout_cluster_sig.rds")
}