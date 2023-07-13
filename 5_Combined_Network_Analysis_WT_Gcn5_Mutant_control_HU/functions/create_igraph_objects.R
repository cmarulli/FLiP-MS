

create_igraph_objects <- function(network_edgelist_file,
                                  FiLiP_changes_WT, FiLiP_changes_M, LiP_changes_WT, LiP_changes_M,  
                                  abundance_changes_WT, abundance_changes_M,
                                  detected_WT, detected_M, 
                                  SAGA_targets, Gcn5_acetylation_targets, high_confidence_markers_WT, high_confidence_markers_M){

  ##### load data
  network_edgelist <- as.matrix(fread(network_edgelist_file))
  
  #### create igraph object
  network_full <- graph_from_edgelist(network_edgelist, directed = FALSE)
  
  #### simplify
  network_full <- simplify(network_full, remove.multiple = TRUE, remove.loops = TRUE)
  
  #### plot the network as force directed network
  #layout_full <- layout_with_drl(network_full,  )
  layout_full <- layout_with_graphopt(network_full,  )
  
  
  #### define network plot properties
  # V(network_full)$size <- 150
  V(network_full)$size <- 600
  V(network_full)$color <- "grey75"
  V(network_full)$frame.color <- "grey75"
  V(network_full)$label  <- ""
  E(network_full)$color <- "grey90"
  
  
  #### get limits for plot
  upper <- 1.1 * apply(layout_full, 2, max)
  lower <- 1.1 * apply(layout_full, 2, min)
  
  pdf("plots/full_input_network.pdf", pointsize = 2)
    plot(network_full, layout = layout_full, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)   # need to adjust limits with new data
    title("Full input network", cex.main = 2)
  dev.off()
  
  
  #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # annotate the network with information from the experiment
  #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  ##### add info of abundance changes WT 
  index_changes <- which(V(network_full)$name %in% abundance_changes_WT)
  V(network_full)$Abundance_changes_WT <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Abundance_changes_WT[index_changes] <- TRUE
  
  ##### add info of abundance changes M 
  index_changes <- which(V(network_full)$name %in% abundance_changes_M)
  V(network_full)$Abundance_changes_M <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Abundance_changes_M[index_changes] <- TRUE
  
  ##### add info of LiP changes WT
  index_changes <- which(V(network_full)$name %in% LiP_changes_WT)
  V(network_full)$LiP_changes_WT <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$LiP_changes_WT[index_changes] <- TRUE
  
  ##### add info of LiP changes M
  index_changes <- which(V(network_full)$name %in% LiP_changes_M)
  V(network_full)$LiP_changes_M <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$LiP_changes_M[index_changes] <- TRUE
  
  ##### add info from FiLiP marker changes WT
  index_changes <- which(V(network_full)$name %in% FiLiP_changes_WT)
  V(network_full)$FiLiP_changes_WT <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$FiLiP_changes_WT[index_changes] <- TRUE
  
  ##### add info from FiLiP marker changes M
  index_changes <- which(V(network_full)$name %in% FiLiP_changes_M)
  V(network_full)$FiLiP_changes_M <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$FiLiP_changes_M[index_changes] <- TRUE
  
  ##### add info about high confidence FiLiP change WT
  index_changes <- which(V(network_full)$name %in% high_confidence_markers_WT)
  V(network_full)$high_confidence_FiLiP_WT <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$high_confidence_FiLiP_WT[index_changes] <- TRUE
  
  ##### add info about high confidence FiLiP change M
  index_changes <- which(V(network_full)$name %in% high_confidence_markers_M)
  V(network_full)$high_confidence_FiLiP_M <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$high_confidence_FiLiP_M[index_changes] <- TRUE
  
  
  ##### add info about protein detection WT
  index_changes <- which(V(network_full)$name %in% detected_WT)
  V(network_full)$Detected_WT <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Detected_WT[index_changes] <- TRUE
  
  ##### add info about protein detection M
  index_changes <- which(V(network_full)$name %in% detected_M)
  V(network_full)$Detected_M <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Detected_M[index_changes] <- TRUE
  
  
  ##### add info of being a SAGA transcription target
  index_changes <- which(V(network_full)$name %in% SAGA_targets)
  V(network_full)$SAGA_targets <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$SAGA_targets[index_changes] <- TRUE
  
  
  ##### add info of being a Gcn5 acetylation target
  index_changes <- which(V(network_full)$name %in% Gcn5_acetylation_targets)
  V(network_full)$Gcn5_acetylation_targets <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Gcn5_acetylation_targets[index_changes] <- TRUE
  
  #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  
  V(network_full)$frame.color[which(V(network_full)$FiLiP_changes_WT == TRUE)] <- "#0057b7"
  V(network_full)$frame.color[which(V(network_full)$FiLiP_changes_M == TRUE)] <- "#a7171a" 
  V(network_full)$frame.color[which(V(network_full)$FiLiP_changes_WT == TRUE & V(network_full)$FiLiP_changes_M == TRUE)] <- "#ffc900"   
  
  pdf("plots/full_input_network_FiLiP_changes.pdf", pointsize = 2)
    plot(network_full, layout = layout_full, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)
    title("Full input network - FiLiP Changes highlighted", cex.main = 2)
  dev.off()
  
  #### save network and layout
  save(network_full, file = "intermediate/full_network.rds")
  save(layout_full, file = "intermediate/layout_full.rds")
}

