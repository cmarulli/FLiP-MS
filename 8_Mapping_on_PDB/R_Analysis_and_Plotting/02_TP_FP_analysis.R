library(UpSetR)
library(ggplot2)

source("../../../Rfunctions/barplot_numbers_ggplot.R")
source("../../../GitHub_updated/1_FLiP_PBI_Library/Rfunctions/get_GO_annotation.R")
source("../../../Rfunctions/fisher_GO_test.R")

data_mapped <- read.csv("../Output/Average_distance_best.csv", header = TRUE, stringsAsFactors = FALSE)
data_mapped$sign <- sapply(data_mapped$q, function(x) return(ifelse(x < 0.05, TRUE, FALSE)))

length(unique(data_mapped$protein))

TP <- data_mapped[which(data_mapped$sign == TRUE & data_mapped$labels == 1), ]
FP <- data_mapped[which(data_mapped$sign == TRUE & data_mapped$labels == -1), ]

length(unique(data_mapped$protein))
length(unique(TP$protein))

protein_list <- list()
protein_list[["TP"]] <- unique(TP$protein)
protein_list[["FP"]] <- unique(FP$protein)

jpeg("./Plots/Upset_TP_FP_proteins.jpg", width = 1500, height = 1500, res = 400)
upset(fromList(protein_list))
dev.off()

overlap_proteins <- intersect(protein_list[[1]], protein_list[[2]])

TP_only <- TP[-which(TP$protein %in% overlap_proteins), ]
FP_only <- FP[-which(FP$protein %in% overlap_proteins), ]
FP_only <- FP_only[-which(FP_only$protein == "P14832"), ]

TP_overlap <- TP[which(TP$protein %in% overlap_proteins), ]
FP_overlap <-  FP[which(FP$protein %in% overlap_proteins), ]

overlap <- rbind(TP_overlap, FP_overlap)
only <- rbind(TP_only, FP_only)

sort(table(overlap$protein))

barplot_numbers_ggplot(c("TP", "FP", "TPO", "FPO"), c(nrow(TP_only), nrow(FP_only), nrow(TP_overlap), nrow(FP_overlap)), file_location = "./Plots/Number_peptides_TP_FP_class.pdf", y_lim = 1800)

data_mapped$type <- rep(NA, nrow(data_mapped))
data_mapped[which(data_mapped$peptide %in% TP_only$peptide), "type"] <- "TP_only"
data_mapped[which(data_mapped$peptide %in% FP_only$peptide), "type"] <- "FP_only"
data_mapped[which(data_mapped$peptide %in% TP_overlap$peptide), "type"] <- "TP_overlap"
data_mapped[which(data_mapped$peptide %in% FP_overlap$peptide), "type"] <- "FP_overlap"

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# GO enrichment analysis of InterPro domains in the different subgroups
#-----------------------------------------------------------------------------------------------------------------------------------------------------
GO_annotation <- read.csv("../../1_FLiP_PBI_Library/OutputData/FLiP_PBI_Library_lip_data_GO.csv", sep = ";")
data_mapped_GO <- merge(data_mapped, GO_annotation, by.x = 'peptide', by.y = "Peptide")
data_mapped_GO <- data_mapped_GO[which(data_mapped_GO$GO_type %in% c("domain", "repeat")), ]

# FP_only
data_mapped_GO$sign <- rep(FALSE, nrow(data_mapped_GO))
data_mapped_GO[which(data_mapped_GO$type  == "FP_only"), "sign"] <- TRUE
data_mapped_GO_sub <- data_mapped_GO[-which(data_mapped_GO$type %in% c("TP_only", "TP_overlap", "FP_overlap")), ]
p_values_GO_terms_domain <- fisher_GO_test(data_mapped_GO_sub)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(data_mapped_GO_sub$GO_names == x & data_mapped_GO_sub$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/FP_only_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 21 , height = 16, units = "cm")

# TP_only
data_mapped_GO$sign <- rep(FALSE, nrow(data_mapped_GO))
data_mapped_GO[which(data_mapped_GO$type  == "TP_only"), "sign"] <- TRUE
data_mapped_GO_sub <- data_mapped_GO[-which(data_mapped_GO$type %in% c("FP_only", "TP_overlap", "FP_overlap")), ]
p_values_GO_terms_domain <- fisher_GO_test(data_mapped_GO_sub)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(data_mapped_GO_sub$GO_names == x & data_mapped_GO_sub$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/TP_only_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 21 , height = 16, units = "cm")

# TP_overlap
data_mapped_GO$sign <- rep(FALSE, nrow(data_mapped_GO))
data_mapped_GO[which(data_mapped_GO$type  == "TP_overlap"), "sign"] <- TRUE
data_mapped_GO_sub <- data_mapped_GO[-which(data_mapped_GO$type %in% c("TP_only", "FP_only", "FP_overlap")), ]
p_values_GO_terms_domain <- fisher_GO_test(data_mapped_GO_sub)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(data_mapped_GO_sub$GO_names == x & data_mapped_GO_sub$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/TP_overlap_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 21 , height = 16, units = "cm")

# FP_overlap
data_mapped_GO$sign <- rep(FALSE, nrow(data_mapped_GO))
data_mapped_GO[which(data_mapped_GO$type  == "FP_overlap"), "sign"] <- TRUE
data_mapped_GO_sub <- data_mapped_GO[-which(data_mapped_GO$type %in% c("TP_only", "FP_only", "TP_overlap")), ]
p_values_GO_terms_domain <- fisher_GO_test(data_mapped_GO_sub)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(data_mapped_GO_sub$GO_names == x & data_mapped_GO_sub$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/FP_overlap_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 21 , height = 16, units = "cm")


# TP_all
data_mapped_GO$sign <- rep(FALSE, nrow(data_mapped_GO))
data_mapped_GO[which(data_mapped_GO$type %in% c("TP_only", "TP_overlap")), "sign"] <- TRUE
data_mapped_GO_sub <- data_mapped_GO[-which(data_mapped_GO$type %in% c("FP_only", "TP_overlap")), ]
p_values_GO_terms_domain <- fisher_GO_test(data_mapped_GO_sub)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(data_mapped_GO_sub$GO_names == x & data_mapped_GO_sub$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/TP_all_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 21 , height = 16, units = "cm")

# FP_all
data_mapped_GO$sign <- rep(FALSE, nrow(data_mapped_GO))
data_mapped_GO[which(data_mapped_GO$type %in% c("FP_only", "FP_overlap")), "sign"] <- TRUE
data_mapped_GO_sub <- data_mapped_GO[-which(data_mapped_GO$type %in% c("TP_only", "TP_overlap")), ]
p_values_GO_terms_domain <- fisher_GO_test(data_mapped_GO_sub)
p_values_GO_terms_domain <- p_values_GO_terms_domain[which(p_values_GO_terms_domain$p_values < 0.05), ]
p_values_GO_terms_domain$count <- sapply(p_values_GO_terms_domain$GO_Terms, function(x) length(which(data_mapped_GO_sub$GO_names == x & data_mapped_GO_sub$sign == TRUE)))
p_values_GO_terms_domain <- p_values_GO_terms_domain[order(p_values_GO_terms_domain$count, decreasing = FALSE), ]
p_values_GO_terms_domain$GO_Terms <- factor(p_values_GO_terms_domain$GO_Terms, levels = p_values_GO_terms_domain$GO_Terms)

ggplot(p_values_GO_terms_domain, aes(x= count, y = GO_Terms, fill = p_values)) + geom_bar(stat = "identity")  +theme_classic(base_size = 16) +
  xlab("Significant Peptide Count") + ylab("")
ggsave(filename = paste("./Plots/FP_all_GO_enrichment_domain.jpg", sep =""),  device = "jpg", width = 21 , height = 16, units = "cm")
