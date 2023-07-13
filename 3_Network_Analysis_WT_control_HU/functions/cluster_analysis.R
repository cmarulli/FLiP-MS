cluster_analysis <- function(cluster_method = "walktrap", cluster_min_size = 10, cluster_max_size = 500,
                             cluster_col_pal = "Dark2", cluster_q_value_cutoff = 0.05, complex_portal, gene_names){

  ##### load data
  load(file = "intermediate/pageRank_network.rds")
  load( file = "intermediate/layout_pageRank.rds")
  load( file = "intermediate/layout_full.rds")
  load(file = "intermediate/FiLiP_changes.rds")
  
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
  cluster_sizes <- table(V(network_porpagated)$Cluster_membership)
  index_cluster <- which(cluster_sizes >= cluster_min_size & cluster_sizes <= cluster_max_size)
  
  number_of_clusters <-length(index_cluster)
  
  clusters_to_use <- names(cluster_sizes[index_cluster])
  index_cluster <- which(V(network_porpagated)$Cluster_membership %in% clusters_to_use)
  
  network_clustered <- induced_subgraph(network_porpagated, index_cluster)
  layout_clustered <- layout_porpagated[index_cluster, ]
  
  # count the number of FiLiP changes in each cluster and condition
  
  V(network_clustered)$markers_in_cluster <- rep(0, length(V(network_clustered)$name))
  for(c in clusters_to_use){
    nb_markes <- length(which(V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$FiLiP_changes == TRUE))
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$markers_in_cluster <- nb_markes
  }
  
  
  # remove the clusters that do not have a single FiLiP change in either conditions
  
  clusters_to_use <- unique(V(network_clustered)$Cluster_membership[which(V(network_clustered)$markers_in_cluster != 0)])
  index_cluster <- which(V(network_clustered)$Cluster_membership %in% clusters_to_use)
  network_clustered <- induced_subgraph(network_clustered, index_cluster)
  layout_clustered <- layout_clustered[index_cluster, ]
  number_of_clusters <- length(clusters_to_use)
  
  # add the name of the complex from the complex portal database
  
  V(network_clustered)$complex <- rep(NA, length(V(network_clustered)$name))
  df_cluster_name <- data.frame()
  for(c in clusters_to_use){
    # get the proteins in the cluster
    proteins_in_complex_propageted <- V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$name
    # get the corresponding complex names
    complexes_propagated <- as.vector(unlist(sapply(proteins_in_complex_propageted, function(x) return(complex_portal[grep(x, complex_portal$`Identifiers (and stoichiometry) of molecules in complex`), 2]))))
    
    # get the proteins with FiLiP marker changes in cluster
    proteins_in_complex_FiLiP <- V(network_clustered)[which(V(network_clustered)$Cluster_membership == c & V(network_clustered)$FiLiP_changes == TRUE)]$name
    # get the corresponding complex names
    complexes_FiliP <- as.vector(unlist(sapply(proteins_in_complex_FiLiP, function(x) return(complex_portal[grep(x, complex_portal$`Identifiers (and stoichiometry) of molecules in complex`), 2]))))
    
    # extract the complexes for which a FiLiP marker was detected
    unique_complexes <- unique(complexes_FiliP)
    complex_subunits <- as.vector(unlist(sapply(unique_complexes, function(x) return(str_count(as.character(complex_portal[grep(x, complex_portal$`Recommended name`), 5]), "\\|") +1))))
    
    # select the most abundant complex name
    # complex <- names(sort(table(complexes), decreasing = TRUE)[1])
    
    df_complex_summary <- data.frame(Complex = unique_complexes, Subunits = complex_subunits)
    
    # get the number of subunits changing in each complex
    df_complex_summary$detected_subunits_FiLiP <- sapply(df_complex_summary$Complex, function(x) return(length(grep(x, complexes_FiliP))))
    df_complex_summary$detected_subunits_propagated <- sapply(df_complex_summary$Complex, function(x) return(length(grep(x, complexes_propagated))))
    df_complex_summary$fraction_detected <- apply(df_complex_summary, 1, function(x) return(as.numeric(x[3])/as.numeric(x[2])))
    
    df_complex_summary <- df_complex_summary[order(df_complex_summary$detected_subunits_propagated, decreasing = TRUE), ]
    df_complex_summary <- df_complex_summary[order(df_complex_summary$detected_subunits_FiLiP, decreasing = TRUE), ]
    
    #complex <- paste(df_complex_summary[which(df_complex_summary$fraction_detected == max(df_complex_summary$fraction_detected)), "Complex"], collapse = ";")
    #complex <- df_complex_summary[which(df_complex_summary$detected_subunits_FiLiP == max(df_complex_summary$detected_subunits_FiLiP)), "Complex"]
                     
    # if there are more than 2 complexes with the highest number of FiLiP marker changes, chose the two that have the highest number of propagated proteins 
    #if(length(complex)>2){
     #complex <- df_complex_summary[which(df_complex_summary$Complex %in% complex & df_complex_summary$detected_subunits_propagated >= sort(df_complex_summary$detected_subunits_propagated, decreasing = TRUE)[2]), "Complex"]
    #}
    
    # merge into one name

     
    complex <- paste(unlist(apply(df_complex_summary[, c("Complex", "detected_subunits_FiLiP", "detected_subunits_propagated")], 1, function(x) return(paste(x, collapse = "|")))), collapse = ";")
    
    
      
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)[1]]$complex <- complex
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$complex_name <- complex
    
    df_cluster_name <- rbind(df_cluster_name, c(c, as.character(complex)))
  }
  
  df_cluster_name <- as.data.frame(df_cluster_name)
  colnames(df_cluster_name) <- c("Cluster", "ComplexName")
  write.csv(df_cluster_name, "./final/cluster_complex_names.csv", row.names = FALSE, quote = FALSE)
  
  # create a summary df with all the information
  df_summary <- data.frame(V(network_clustered)$name)
  colnames(df_summary) <- c("Protein")
  df_summary$Gene <- sapply(df_summary$Protein, function(x) return(gene_names[grep(x, gene_names$Entry), "Gene.names"]))
  df_summary$Cluster <- V(network_clustered)$Cluster_membership
  df_summary$Complex <- V(network_clustered)$complex_name
  
  df_summary$Detected <-  V(network_clustered)$Detected
  df_summary$Abundance_Change <-  V(network_clustered)$Abundance_changes
  df_summary$LiP_Change <-  V(network_clustered)$LiP_changes
  df_summary$FiLiP_Change <-  V(network_clustered)$FiLiP_changes
  df_summary$high_confidence_FiLiP <-  V(network_clustered)$high_confidence_FiLiP
  df_summary$Gcn5_target <- V(network_clustered)$Gcn5_acetylation_targets
  df_summary$SAGA_transcript <- V(network_clustered)$SAGA_targets
  
  df_summary[which(df_summary$Complex == "PRP19-associated complex"), "Complex"] <- "Splicesome LSM1-7 Pat1"
  write.table(df_summary, "./final/Complex_Summary_WT_control_HU.tsv", quote = FALSE, row.names = FALSE, sep = '\t')
  
  
  
  
  
  
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  mycolors <- colorRampPalette(brewer.pal(8, cluster_col_pal))(number_of_clusters)
  
  membership <- unique(V(network_clustered)$Cluster_membership)
  lookup = setNames(mycolors, as.character(membership))
  
  df <- data.frame(as.character(V(network_clustered)$Cluster_membership))
  colnames(df) <- "membership"
  df <- transform(df, membership=lookup[membership], stringsAsFactors=FALSE)
  
  
  #V(network_clustered)$color <- df[,1]
  V(network_clustered)$color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  V(network_clustered)$size <- 1200
  #V(network_clustered)$frame.color <- df[,1]
  V(network_clustered)$frame.color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  V(network_clustered)[which(V(network_clustered)$FiLiP_changes == TRUE)]$frame.color <- "black"
  
  pdf("plots/clustered_network_size_filter_cluster_label.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors,vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)   
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  pdf("plots/clustered_network_size_filter.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  
  #### change coordinates
  layout_clustered_new <- layout_nicely(network_clustered,  )
  
  upper <- 1.1 * apply(layout_clustered_new, 2, max)
  lower <- 1.1 * apply(layout_clustered_new, 2, min)
  
  
  V(network_clustered)$size <- 35
  pdf("plots/clustered_network_size_filter_new_layout.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  pdf("plots/clustered_network_size_filter_new_layout_labeled.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors,vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)    
  title(paste(as.character(number_of_clusters), " notable communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  
  ##### Test clusters if they are enriched in our hits
  ### Fisher for enrichment of original hits
  fisher_p <- rep(0, length(membership))
  total_number <- length(V(network_clustered))
  
  V(network_clustered)$origial_hit <- rep(0, total_number)
  index_pody <- which(V(network_clustered)$name %in% FiLiP_changes)
  V(network_clustered)$origial_hit[index_pody] <- 1
  
  
  for (i in 1:length(membership)){
    current_cluster <- membership[i]
    
    # define contingency table
    YY <- length(which( V(network_clustered)$origial_hit == 1 & V(network_clustered)$Cluster_membership == current_cluster ))   # in original hit list & in cluster
    NY <- length(which( V(network_clustered)$origial_hit == 0 & V(network_clustered)$Cluster_membership == current_cluster ))   # not in original hit list & in cluster
    YN <- length(which( V(network_clustered)$origial_hit == 1 & V(network_clustered)$Cluster_membership != current_cluster ))   # in original hit list & not in cluster
    NN <- length(which( V(network_clustered)$origial_hit == 0 & V(network_clustered)$Cluster_membership != current_cluster ))   # not in original hit list & not in cluster
    
    # dat <- data.frame(
    #   "in_cluster" = c(YY, NY),
    #   "not_in_cluster" = c(YN, NN),
    #   row.names = c("hit_list", "not_in_hit_list"),
    #   stringsAsFactors = FALSE
    # )
    # 
    
    dat <- data.frame(
      "hit_list" = c(YY, YN),
      "not_in_hit_list" = c(NY, NN),
      row.names = c("in_cluster", "not_in_cluster"),
      stringsAsFactors = FALSE
    )

    
    
    # save p value
    fisher_p[i] <- fisher.test(dat, alternative = "greater" )$p.value
  }
  
  # adjust for multiple testing
  #fisher_q <- p.adjust(fisher_p, method = "BH") 
  fisher_q <- fisher_p
  
  ### Kolmogorov-Smirnov Test to see if higher pageRank score
  ks_p <- rep(0, length(membership))
  
  
  # i <- 9
  for (i in 1:length(membership)){
    current_cluster <- membership[i]
    scores_in_cluster <- V(network_clustered)$pageRank[ which(V(network_clustered)$Cluster_membership == current_cluster) ]
    scors_not_in_cluster <- V(network_clustered)$pageRank[ which(V(network_clustered)$Cluster_membership != current_cluster) ]
    
    #ks_p[i] <- ks.test(scores_in_cluster, scors_not_in_cluster, alternative = "greater")$p.value
    ks_p[i] <- ks.test(scores_in_cluster, scors_not_in_cluster, alternative = "less")$p.value
  }
  
  # adjust for multiple testing
  #ks_q <- p.adjust(ks_p, method = "BH")
  ks_q <- ks_p
  
  cluster_q <- data.frame("cluster" = membership,
                          "Fisher_q" = fisher_p,
                          "Kolmogorov-Smirnov_q" =  ks_q)
  
  
  
  # keep only clusters that have at least one FiLiP marker
  
  clusters <- cluster_q$cluster
  keep_clusters <- c()
  for (c in 1:length(clusters) ){
    if(sum(V(network_clustered)$origial_hit[which(V(network_clustered)$Cluster_membership == clusters[c])]) > 0){
      keep_clusters <- c(keep_clusters, clusters[c]) 
    }
  } 
  
  cluster_q <- cluster_q[which(cluster_q$cluster %in% keep_clusters), ]
  
  ##### reduce to significant cluster and plot
  cluster_q_sig <- cluster_q[which(cluster_q$Fisher_q <= cluster_q_value_cutoff | cluster_q$Kolmogorov.Smirnov_q <=  cluster_q_value_cutoff), ]
  index_sig <- which(V(network_clustered)$Cluster_membership %in% cluster_q_sig$cluster)
  
  network_sig <- induced_subgraph(network_clustered, index_sig)
  layout_sig <- layout_clustered_new[index_sig, ]
  
  
  clusters_sig <- V(network_sig)$Cluster_membership
  list_of_clusters <- unique(clusters_sig)
  
  
  index_sig <- which(V(network_clustered)$Cluster_membership %in% list_of_clusters)
  network_sig <- induced_subgraph(network_clustered, index_sig)
  layout_sig <- layout_clustered_new[index_sig, ]
  
  
  
  pdf("plots/clustered_network_significant_new_layout.pdf", pointsize = 2)
  plot(network_sig, layout = layout_sig, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(length(list_of_clusters)), " significant communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  # new layout
  #layout_clustered_sig_new <- layout_with_lgl(network_sig, )
  #V(network_sig)$size <- 250
  layout_clustered_sig_new <- layout_nicely(network_sig, )
  V(network_sig)$size <- 35
  
  
  upper <- apply(layout_clustered_sig_new, 2, max)
  lower <- apply(layout_clustered_sig_new, 2, min)
  
  
  pdf("plots/clustered_network_significant_final_layout.pdf", pointsize = 2)
  plot(network_sig, layout = layout_clustered_sig_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(length(list_of_clusters)), " significant communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
  pdf("plots/clustered_network_significant_final_layout_label.pdf", pointsize = 2)
  plot(network_sig, layout = layout_clustered_sig_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors,  vertex.label = V(network_sig)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)    
  title(paste(as.character(length(list_of_clusters)), " significant communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
  
 
  
  # old layout
  index_sig <- which(V(network_porpagated)$Cluster_membership %in% list_of_clusters)
  network_sig <- induced_subgraph(network_porpagated, index_sig)
  layout_sig <- layout_porpagated[index_sig, ]

  V(network_sig)$color <- lookup[as.character(V(network_sig)$Cluster_membership)]
  V(network_sig)$frame.color <- lookup[as.character(V(network_sig)$Cluster_membership)]
  V(network_sig)$size <- 1200
  V(network_sig)[which(V(network_sig)$FiLiP_changes == TRUE)]$frame.color <- "#ffc900"
  
    
    
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
    
    
  pdf("plots/clustered_network_significant_original_layout.pdf", pointsize = 2)
  plot(network_sig, layout = layout_sig, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)
  title(paste(as.character(length(list_of_clusters)), " significant communities found by ",  cluster_method, " clustering"), cex.main = 2)
  dev.off()
    
  
  ###### plot LiP anf FLiP changes
  
  
  V(network_sig)$color <- "grey80"
  V(network_sig)$frame.color <- "grey80"
  V(network_sig)$size <- 1200
  V(network_sig)[which(V(network_sig)$FiLiP_changes == TRUE)]$frame.color <- "#ffc900"
  V(network_sig)[which(V(network_sig)$LiP_changes == TRUE)]$color <- "#0057b7"
      
  
  
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  
  pdf("plots/clustered_network_significant_FLiP_and_LiP_changes.pdf", pointsize = 2)
  plot(network_sig, layout = layout_sig, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)
  title("FLiP changes (yellow) and LiP changes (blue)", cex.main = 2)
  dev.off()
  
  
  ###### plot SAGA target and abundance changes
  
  
  V(network_sig)$color <- "grey80"
  V(network_sig)$frame.color <- "grey80"
  V(network_sig)$size <- 1200
  V(network_sig)[which(V(network_sig)$SAGA_targets == TRUE)]$frame.color <- "#ffc900"
  V(network_sig)[which(V(network_sig)$Abundance_changes == TRUE)]$color <- "#0057b7"
          
        
        
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
        
        
  pdf("plots/clustered_network_significant_SAGA_target_and_Abundance_changes.pdf", pointsize = 2)
  plot(network_sig, layout = layout_sig, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)
  title("SAGA transcription targets (yellow) and abundance changes (blue)", cex.main = 2)
  dev.off()
  
  
  
  ###### plot Gcn5 acetylation target and FiLiP changes
  
  
  V(network_sig)$color <- "grey80"
  V(network_sig)$frame.color <- "grey80"
  V(network_sig)$size <- 1200
  V(network_sig)[which(V(network_sig)$FiLiP_changes == TRUE)]$frame.color <- "#ffc900"
  V(network_sig)[which(V(network_sig)$Gcn5_acetylation_targets == TRUE)]$color <- "#0057b7"
          
        
        
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
        
        
  pdf("plots/clustered_network_significant_FiliP_change_and_Gcn5_acetylation_target.pdf", pointsize = 2)
  plot(network_sig, layout = layout_sig, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)
  title("FiLiP change (yellow) and Gcn5 acetylation targets (blue)", cex.main = 2)
  dev.off()
        
        
  
  
  
  ##### export results
  write.table(cluster_q, "final/cluster_q_values.tsv", quote = F , sep = "\t", row.names = F)
  save(network_sig, file = "final/cluster_sig_network.rds")
  save(layout_clustered_sig_new, file = "final/layout_cluster_sig.rds")
  save(layout_sig, file = "final/layout_sig_original.rds")
}
