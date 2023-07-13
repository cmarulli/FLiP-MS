create_igraph_objects <- function(p_cutoff = 0.05, filter_column = "q", add_non_connected = FALSE,
                                  network_edgelist_file, detected, FiLiP_changes, LiP_changes, abundance_changes, SAGA_targets, Gcn5_acetylation_targets, high_confidence_markers){

  ##### load data
  network_edgelist <- as.matrix(fread(network_edgelist_file))
  
  #### create igraph object
  network_full <- graph_from_edgelist(network_edgelist, directed = FALSE)
  
  #### simplify
  network_full <- simplify(network_full, remove.multiple = TRUE, remove.loops = TRUE)
  
  
  #### add not-connected vertices from p-body data
  # if (add_non_connected){
  #   index_pbody <- which(FiLiP_changes %in% V(network_full)$name)
  #   missing_v <- FiLiP_changes[-index_pbody]
  #   
  #   network_full <- add_vertices(network_full, length(missing_v), name = missing_v)
  #   
  # }
  # 
  
  #### plot the network as force directed network
  # layout_full <- layout_with_drl(network_full,  )
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
  
  
  
  ##### add info from experiment

  ##### add info of protein detection
  
  index_changes <- which(V(network_full)$name %in% detected)
  V(network_full)$Detected <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Detected[index_changes] <- TRUE
  
  
  ##### add info of abundance changes
  
  index_changes <- which(V(network_full)$name %in% abundance_changes)
  V(network_full)$Abundance_changes <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Abundance_changes[index_changes] <- TRUE
  
  ##### add info of LiP changes
  
  index_changes <- which(V(network_full)$name %in% LiP_changes)
  V(network_full)$LiP_changes <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$LiP_changes[index_changes] <- TRUE
  
  
  ##### add info of being a SAGA transcription target
  
  index_changes <- which(V(network_full)$name %in% SAGA_targets)
  V(network_full)$SAGA_targets <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$SAGA_targets[index_changes] <- TRUE
  
  
  ##### add info of being a Gcn5 acetylation target
  
  index_changes <- which(V(network_full)$name %in% Gcn5_acetylation_targets)
  V(network_full)$Gcn5_acetylation_targets <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$Gcn5_acetylation_targets[index_changes] <- TRUE
  
  ##### add info of being a high confidence FLiP marker
  
  index_changes <- which(V(network_full)$name %in% high_confidence_markers)
  V(network_full)$high_confidence_FiLiP <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$high_confidence_FiLiP[index_changes] <- TRUE
  
  
  ##### add info from FLiP marker changes
  
  
  index_changes <- which(V(network_full)$name %in% FiLiP_changes)
  V(network_full)$FiLiP_changes <- rep("FALSE", length(V(network_full)$name))
  V(network_full)$FiLiP_changes[index_changes] <- TRUE
  V(network_full)$size[index_changes] <- 1200
  V(network_full)$frame.color[index_changes] <- "#ffc900"
  V(network_full)$size[index_changes] <- 1200
  
  
  pdf("plots/full_input_network_FiLiP_changes.pdf", pointsize = 2)
    plot(network_full, layout = layout_full, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)
    title("Full input network - FiLiP Changes highlighted", cex.main = 2)
  dev.off()
  
  
  #### plot p-body part only
  network_pbody <- induced_subgraph(network_full, index_changes)
  layout_pbody <- layout_full[index_changes, ]
  
  pdf("plots/FiLiP_changes_network.pdf", pointsize = 2)
    plot(network_pbody, layout = layout_pbody, rescale = FALSE ,ylim=c(lower[2],upper[2]),xlim=c(lower[1],upper[1]), asp = 0)   
    title("FiLiP Changes input network", cex.main = 2)
  dev.off()
  
  
  
  #### save network and layout
  save(network_full, file = "intermediate/full_network.rds")
  save(layout_full, file = "intermediate/layout_full.rds")
  save(FiLiP_changes, file = "intermediate/FiLiP_changes.rds")

}