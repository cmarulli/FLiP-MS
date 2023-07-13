library(protti)
library(ggplot2)

data <- read.csv("./OutputData/LiP_WT_control_HU_anova_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")

get_barcode <- function(data, protein_ID, protein_name = "", q_cut, L2FC_cut){
  
  my_prot <- data[which(data$Protein == protein_ID), c("Peptide", "Protein", "Position", "q", "L2FC", "has_FLiP")]
  my_prot$Position <- as.numeric(my_prot$Position)
  my_prot$stop <- as.numeric(apply(my_prot, 1, function(x) return(as.numeric(x[3]) + nchar(x[1]))))
  my_prot$length <- as.numeric(rep( nchar(data[which(data$Protein == protein_ID), "Sequence"][1]), nrow(my_prot)))
  
  barcode_plot(my_prot, start_position = Position, end_position = stop, protein_length = length, facet = Protein, cutoffs = c(L2FC = 0, q = 0.05))
  ggsave(paste("./Plots/Barcode_", protein_name, ".pdf", sep = ''),  device = "pdf", width = 10, height = 2)
}

get_barcode(data, "P25644", "Pat1", 0.05, 0)
get_barcode(data, "P39517", "Dhh1", 0.05, 0)
get_barcode(data, "P22147", "Xrn1", 0.05, 0)
get_barcode(data, "P53550", "Dcp2", 0.05, 0)
