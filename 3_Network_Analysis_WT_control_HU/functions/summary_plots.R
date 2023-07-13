summary_plots <- function(){
  ##### load data
  load(file = "intermediate/FiLiP_changes.rds")
  propagted_proteins <- fread("intermediate/propagted_proteins.tsv")
  load(file = "final/cluster_sig_network.rds")
  load(file = "intermediate/full_network.rds")
  
  
  
  #### plot 1: number of proteins at different steps
  number_of_input_hits <- length(intersect(FiLiP_changes, names(V(network_full))))
  number_of_propagted_proteins <- length(unique(propagted_proteins$x))
  number_of_proteins_in_sig_custers <- length(unique(V(network_sig)$name))
  
  plot_1_data <- data.frame("category" = c("PBI Marker Change",
                                           "Propagated Network",
                                           "Final Clusters"),
                            "number_of_proteins" = c(number_of_input_hits,
                                                     number_of_propagted_proteins,
                                                     number_of_proteins_in_sig_custers))
  
  plot_1_data$category <- factor(plot_1_data$category, levels = plot_1_data$category)
  
  p <- ggplot(plot_1_data, aes(x = category, y = number_of_proteins)) +
    geom_bar(stat="identity", fill = c("#0057b7", "#0057b7", "#ffc900")) +
    geom_text(aes(label=number_of_proteins), vjust=1.6, color="white", size=8)+
    theme_classic() +
    xlab( "Number of proteins in ...") +
    ggtitle("Number of proteins at key steps") +
    theme( axis.title.y=element_blank(), legend.title = element_blank(),
           axis.title.x = element_text(size =16), text = element_text(size = 20))
  ggsave("plots/number_of_proteins_at_key_steps.pdf", p, width = 21, height = 21, units = "cm")
  
  
  ### plot overview_ cluster size
  plot_2_data <- as.data.frame(table(V(network_sig)$Cluster_membership))
  
    p <- ggplot(plot_2_data, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = plot_2_data$Var1) +
    geom_text(aes(label=Freq), vjust=1.6, color="white", size=8)+
    theme_classic() +
    xlab( "Cluster") +
    ylab("Number of proteins per cluster") +
    theme( axis.title.y=element_blank(), legend.position = "none",
           axis.title.x = element_text(size =16), text = element_text(size = 20), axis.text = element_text(size = 24))
  ggsave("plots/number_of_proteins_per_cluster.pdf", p, width = 60, height = 21, units = "cm")
  
  
}
