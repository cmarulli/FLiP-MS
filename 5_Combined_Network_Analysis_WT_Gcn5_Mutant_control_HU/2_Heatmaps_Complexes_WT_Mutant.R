library(ztable)
library(reshape2)
library(magrittr)
library(ggplot2)

data_input <- read.table("./final/Complex_Summary_WT_Mutant_control_HU.tsv", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
data <- data_input[which(data_input$FiLiP_Change_WT | data_input$FiLiP_Change_M),  c("Complex", "Gene", "FiLiP_Change_WT", "FiLiP_Change_M")]
data$diff <- apply(data[, c(3,4)], 1, function(x) return(x[1]-x[2]) )
data$Gene <- factor(data$Gene, levels = data$Gene)

complex_diff <- data.frame()
for(c in unique(data$Complex)){
  
  # number of marker changes in WT and M
  WT <- sum(as.numeric(data[which(data$Complex == c), "FiLiP_Change_WT"]))
  M <- sum(as.numeric(data[which(data$Complex == c), "FiLiP_Change_M"]))
  
  # log2 fold change difference of makers between WT and mutant
  diff <- log2(WT/M)

  # logical telling whether markers are only lost in one condition or also gained
  same_markers <- if(diff > 0){
    all(data[which(data$Complex == c), "diff"] >= 0)
  }else if(diff < 0){
    all(data[which(data$Complex == c), "diff"] <= 0)
  }else if(diff == 0){
    all(data[which(data$Complex == c), "diff"] == 0)
  }
  
  # number of gcn5 targets in the complex
  gcn5_target <- sum(as.numeric(data_input[which(data_input$Complex == c), "Gcn5_target"]))
  complex_diff <- rbind(complex_diff, c(c, WT, M, diff, same_markers, gcn5_target))
  
}

colnames(complex_diff) <- c("Complex", "WT", "Mutant", "Complex_Diff", "Same_Markers", "Gcn5_target")
complex_diff$Complex_Diff <- as.numeric(complex_diff$Complex_Diff)
complex_diff <- complex_diff[order(complex_diff$Complex, decreasing = TRUE), ]
complex_diff <- complex_diff[order(complex_diff$WT), ]
complex_diff <- complex_diff[order(complex_diff$M), ]
complex_diff <- complex_diff[order(complex_diff$Same_Markers), ]
complex_diff <- complex_diff[order(complex_diff$Complex_Diff), ]
complex_diff$Complex <- factor(complex_diff$Complex, levels = complex_diff$Complex)

complex_diff[which(is.infinite(complex_diff$Complex_Diff)), "Complex_Diff"] <- sign(complex_diff[which(is.infinite(complex_diff$Complex_Diff)), "Complex_Diff"]) *2

complex_diff <- melt(complex_diff, id.vars = "Complex")
complex_diff[which(complex_diff$value == TRUE), "value"] <- 0
complex_diff[which(complex_diff$value == FALSE), "value"] <- 1
complex_diff$value <- as.numeric(complex_diff$value)

complex_diff$scaled_value <- rep(0, nrow(complex_diff))
complex_diff[which(complex_diff$variable == "WT"), "scaled_value"] <- complex_diff[which(complex_diff$variable == "WT"), "value"]/10
complex_diff[which(complex_diff$variable == "Mutant"), "scaled_value"] <- -complex_diff[which(complex_diff$variable == "Mutant"), "value"]/7
complex_diff[which(complex_diff$variable == "Complex_Diff"), "scaled_value"] <- complex_diff[which(complex_diff$variable == "Complex_Diff"), "value"]/2
complex_diff[which(complex_diff$variable == "Complex_Diff"), "value"] <- rep(NA, nrow(complex_diff[which(complex_diff$variable == "Complex_Diff"), ]))

ggplot(complex_diff, aes(x = variable, y = Complex, fill = scaled_value)) +
  geom_tile(show.legend = FALSE) + scale_fill_gradient2(low = "#ffc900", mid ="grey90", high = "#0057b7")  + theme_minimal(base_size = 15) + xlab("") + ylab("") +
  geom_text(aes(label = round(value, 1)), ) + theme(axis.text.x=element_blank(), axis.ticks.x = element_blank())
ggsave('./plots/Heatmap_summary_complexes_WT_Mutant_control_HU.pdf', device = "pdf", width = 120 , height = 40, units = "cm")