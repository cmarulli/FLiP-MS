cluster_analysis <- function(cluster_min_size = 10, cluster_max_size = 500,
                             cluster_col_pal = "Dark2", cluster_q_value_cutoff = 0.05, complex_portal){

  ##### load data
  load(file = "intermediate/pageRank_network_merged.rds")
  load( file = "intermediate/layout_pageRank.rds")
  load( file = "intermediate/layout_full.rds")
  
  cluster <- read.table("intermediate/walktrap_clusters.tsv", stringsAsFactors = FALSE)
  V(network_propagated)$Cluster_membership <- cluster$V1 
  
  #### filter clusters based on size
  cluster_sizes <- table(V(network_propagated)$Cluster_membership)
  index_cluster <- which(cluster_sizes >= cluster_min_size & cluster_sizes <= cluster_max_size)
  number_of_clusters <- length(index_cluster)
  clusters_to_use <- names(cluster_sizes[index_cluster])
  index_cluster <- which(V(network_propagated)$Cluster_membership %in% clusters_to_use)
  
  network_clustered <- induced_subgraph(network_propagated, index_cluster)
  layout_clustered <- layout_propagated[index_cluster, ]
  
  
  # count the number of FiLiP changes in each cluster and condition
  
  V(network_clustered)$markers_in_cluster_WT <- rep(0, length(V(network_clustered)$name))
  for(c in clusters_to_use){
    nb_markes <- length(which(V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$FiLiP_changes_WT == TRUE))
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$markers_in_cluster_WT <- nb_markes
  }
  
  V(network_clustered)$markers_in_cluster_M <- rep(0, length(V(network_clustered)$name))
  for(c in clusters_to_use){
    nb_markes <- length(which(V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$FiLiP_changes_M == TRUE))
    V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$markers_in_cluster_M <- nb_markes
  }

  
  
  # remove the clusters that do not have a single FiLiP change in either conditions
  
  clusters_to_use <- unique(V(network_clustered)$Cluster_membership[which(V(network_clustered)$markers_in_cluster_WT != 0 | V(network_clustered)$markers_in_cluster_M != 0)])
  index_cluster <- which(V(network_clustered)$Cluster_membership %in% clusters_to_use)
  network_clustered <- induced_subgraph(network_clustered, index_cluster)
  layout_clustered <- layout_clustered[index_cluster, ]
  number_of_clusters <- length(clusters_to_use)
  
  
  # add the name of the complex from the complex portal database
  # 
  # V(network_clustered)$complex <- rep(NA, length(V(network_clustered)$name))
  # df_cluster_name <- data.frame()
  # for(c in clusters_to_use){
  #   # get the proteins in the cluster
  #   proteins_in_complex <- V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$name
  #   # get the corresponding complex names
  #   complexes <- as.vector(unlist(sapply(proteins_in_complex, function(x) return(complex_portal[grep(x, complex_portal$`Identifiers (and stoichiometry) of molecules in complex`), 2]))))
  #   # select the most abundant compes name
  #   complex <- names(sort(table(complexes), decreasing = TRUE)[1])
  # 
  #   V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$complex <- complex
  #   df_cluster_name <- rbind(df_cluster_name, c(c, as.character(complex)))
  # }
  # 
  # 
  # colnames(df_cluster_name) <- c("Cluster", "ComplexName")
  # write.csv(df_cluster_name, "./final/cluster_complex_names.csv", row.names = FALSE, quote = FALSE)
  
  # add the name of the complex from the complex portal database
  
  V(network_clustered)$complex <- rep(NA, length(V(network_clustered)$name))
  df_cluster_name <- data.frame()
  for(c in clusters_to_use){
    # get the proteins in the cluster
    proteins_in_complex_propageted <- V(network_clustered)[which(V(network_clustered)$Cluster_membership == c)]$name
    # get the corresponding complex names
    complexes_propagated <- as.vector(unlist(sapply(proteins_in_complex_propageted, function(x) return(complex_portal[grep(x, complex_portal$`Identifiers (and stoichiometry) of molecules in complex`), 2]))))
    
    # get the proteins with FiLiP marker changes in WT or mutant in cluster
    proteins_in_complex_FiLiP <- V(network_clustered)[which(V(network_clustered)$Cluster_membership == c & (V(network_clustered)$FiLiP_changes_WT == TRUE | V(network_clustered)$FiLiP_changes_M == TRUE))]$name
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
  
  df_summary$Detected_WT <-  V(network_clustered)$Detected_WT
  df_summary$Detected_M <-  V(network_clustered)$Detected_M
  df_summary$Abundance_Change_WT <-  V(network_clustered)$Abundance_changes_WT
  df_summary$Abundance_Change_M <-  V(network_clustered)$Abundance_changes_M
  df_summary$LiP_Change_WT <-  V(network_clustered)$LiP_changes_WT
  df_summary$LiP_Change_M <-  V(network_clustered)$LiP_changes_M
  df_summary$FiLiP_Change_WT <-  V(network_clustered)$FiLiP_changes_WT
  df_summary$FiLiP_Change_M <-  V(network_clustered)$FiLiP_changes_M
  df_summary$FiLiP_Change_WT <-  V(network_clustered)$FiLiP_changes_WT
  df_summary$FiLiP_Change_M <-  V(network_clustered)$FiLiP_changes_M
  df_summary$high_confidence_FiLiP_WT <-  V(network_clustered)$high_confidence_FiLiP_WT
  df_summary$high_confidence_FiLiP_M <-  V(network_clustered)$high_confidence_FiLiP_M
  df_summary$Gcn5_target <- V(network_clustered)$Gcn5_acetylation_targets
  df_summary$SAGA_transcript <- V(network_clustered)$SAGA_targets
  

  write.table(df_summary, "./final/Complex_Summary_WT_Mutant_control_HU.tsv", quote = FALSE, row.names = FALSE, sep = '\t')
  
  
  
  
  
  
  # get the difference in FLiP hits in WT vs mutant for all clusters
  
  V(network_clustered)$difference_markers <- sign(V(network_clustered)$markers_in_cluster_WT - V(network_clustered)$markers_in_cluster_M)
  
  fine = 1000
  
  pal_WT = colorRampPalette(c("white", "#0057b7","#001e64"))
  graphCol_WT  <- pal_WT(fine)[as.numeric(cut(V(network_clustered)$markers_in_cluster_WT, breaks = fine))]
  
  pal_M = colorRampPalette(c("white", "red","#a7171a"))
  graphCol_M  <- pal_M(fine)[as.numeric(cut(V(network_clustered)$markers_in_cluster_M, breaks = fine))]
  
  pal_WT_M = colorRampPalette(c( "#FFB6B6", "white","#99BCE2"))
  graphCol_WT_M  <- pal_WT_M(3)[as.numeric(cut(V(network_clustered)$difference_markers, breaks = 3))]
  
  
  
  
  
  
  
  #### plot network with clusters
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  mycolors <- colorRampPalette(brewer.pal(8, cluster_col_pal))(number_of_clusters)
  
  membership <- unique(V(network_clustered)$Cluster_membership)
  lookup = setNames(mycolors, as.character(membership))
  
  df <- data.frame(as.character(V(network_clustered)$Cluster_membership))
  colnames(df) <- "membership"
  df <- transform(df, membership=lookup[membership], stringsAsFactors=FALSE)
  
  V(network_clustered)$color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  V(network_clustered)$size <- 1200
  V(network_clustered)$frame.color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
        
  
  pdf("plots/clustered_network_size_filter.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors, vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)  
  title(paste(as.character(number_of_clusters), " notable communities found by walktrap clustering"), cex.main = 2)
  dev.off()
  
  
  
  
  
  #### change coordinates
  layout_clustered_new <- layout_nicely(network_clustered,  )
  
  upper <- 1.1 * apply(layout_clustered_new, 2, max)
  lower <- 1.1 * apply(layout_clustered_new, 2, min)
  
  
  V(network_clustered)$size <- 35
  pdf("plots/clustered_network_size_filter_new_layout.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors, vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)  
  title(paste(as.character(number_of_clusters), " notable communities found by walktrap clustering"), cex.main = 2)
  dev.off()
  
  V(network_clustered)$size[which(V(network_clustered)$Gcn5_acetylation_targets == TRUE)] <- 50
  V(network_clustered)$color[which(V(network_clustered)$Gcn5_acetylation_targets == TRUE)] <- "#ffc900"
    
  pdf("plots/clustered_network_size_filter_new_layout_Gcn5_targets.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors, vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)  
  title(paste(as.character(number_of_clusters), " notable communities found by walktrap clustering with Gcn5 Targets"), cex.main = 2)
  dev.off()
  
  
  
  V(network_clustered)$frame.color[which(V(network_clustered)$FiLiP_changes_WT == TRUE)] <- "black"
    
  pdf("plots/clustered_network_size_filter_new_layout_Gcn5_targets_WT_FliP_markers.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(number_of_clusters), " notable communities found by walktrap clustering with Gcn5 Targets with WT FLiP Markers"), cex.main = 2)
  dev.off()
  
  V(network_clustered)$size <- 35
  V(network_clustered)$color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
    
  pdf("plots/clustered_network_size_filter_new_layout_WT_FliP_markers.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(number_of_clusters), " notable communities found by walktrap clustering with WT FLiP Markers"), cex.main = 2)
  dev.off()
 
  V(network_clustered)$frame.color[which(V(network_clustered)$FiLiP_changes_WT == TRUE)] <- "#ffc900"
  
  pdf("plots/clustered_network_size_filter_new_layout_WT_FliP_markers_yellow.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0,
       col=mycolors)  
  title(paste(as.character(number_of_clusters), " notable communities found by walktrap clustering with WT FLiP Markers"), cex.main = 2)
  dev.off()
  
  
  # plot WT hits
  
  V(network_clustered)$size <- 35
  V(network_clustered)$color <- graphCol_WT
  V(network_clustered)$frame.color <- "grey75"
  V(network_clustered)$frame.color[which(V(network_clustered)$FiLiP_changes_WT == TRUE)] <- "black"
      
  pdf("plots/clustered_network_WT.pdf", pointsize = 2)
    plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0)  
    title("WT FiLiP Hits", cex.main = 2)
  dev.off()
  
  # plot mutant hits
  
  V(network_clustered)$size <- 35
  V(network_clustered)$color <- graphCol_M
  V(network_clustered)$frame.color <- "grey75"
  V(network_clustered)$frame.color[which(V(network_clustered)$FiLiP_changes_M == TRUE)] <- "black"
      
  pdf("plots/clustered_network_M.pdf", pointsize = 2)
    plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0)  
    title("Mutant FiLiP Hits", cex.main = 2)
  dev.off()
  
  
  # plot wild type minus mutant hits
  
  V(network_clustered)$size <- 35
  V(network_clustered)$color <- graphCol_WT_M
  V(network_clustered)$frame.color <- "grey75"
  V(network_clustered)$frame.color[which(V(network_clustered)$FiLiP_changes_WT == TRUE)] <- "#0057b7"
  V(network_clustered)$frame.color[which(V(network_clustered)$FiLiP_changes_M == TRUE)] <- "#a7171a"
  V(network_clustered)$frame.color[which(V(network_clustered)$FiLiP_changes_WT == TRUE & V(network_clustered)$FiLiP_changes_M == TRUE) ] <- "black"
      
  
 
  pdf("plots/clustered_network_WT_vs_M.pdf", pointsize = 2)
   plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0)  
   title("Difference WT vs Mutant", cex.main = 2)
  dev.off()
 
  pdf("plots/clustered_network_WT_vs_M_labelled.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0, 
       vertex.label = V(network_clustered)$complex, vertex.label.family="Helvetica", vertex.label.color= "black", vertex.label.cex=0.5)  
  title("Difference WT vs Mutant", cex.main = 2)
  dev.off()
  
  
  
  # plot wild type minus mutant hits
  
  V(network_clustered)$size[which(V(network_clustered)$Gcn5_acetylation_targets == TRUE)] <- 50
  V(network_clustered)$color[which(V(network_clustered)$Gcn5_acetylation_targets == TRUE)] <- "#ffc900"

  pdf("plots/clustered_network_WT_vs_M_acetylation.pdf", pointsize = 2)
  plot(network_clustered, layout = layout_clustered_new, rescale = FALSE , ylim=c(lower[2],upper[2]), xlim=c(lower[1],upper[1]), asp = 0)  
  title("Difference WT vs Mutant with Gcn5 Targets", cex.main = 2)
  dev.off()
        
  
 
  
  V(network_clustered)$color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
  V(network_clustered)$frame.color <- lookup[as.character(V(network_clustered)$Cluster_membership)]
        
  
  ##### export results
  save(layout_clustered_new, file = "final/layout_clustered_new.rds")
  save(network_clustered, file = "final/clustered_network_final.rds")
  
  
  
   
  
  
}

