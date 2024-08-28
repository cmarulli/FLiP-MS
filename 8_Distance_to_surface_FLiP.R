library(ggplot2)
library(tidyr)
library(ggpubr)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FLiP dataset: are PK cleavage sites of semi-tryptic peptides that significantly change located closer to the surface than the ones that do no change? 
# surface is generated with a rolling ball alrgithm with different ball radii, the bigger the radius, the lower the chance to classify a burried site as being located at the surface
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
distances_st_1.5 <- read.csv("./Distance_to_surface/semi-tryptic/Flip_distance_surface.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_ft_1.5 <- read.csv("./Distance_to_surface/fully-tryptic/flip_distance_surface_fully.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_1.5 <- rbind(distances_st_1.5, distances_ft_1.5)

distances_st_4 <- read.csv("./Distance_to_surface/semi-tryptic/flip_distance_surface_4.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_ft_4 <- read.csv("./Distance_to_surface/fully-tryptic/flip_distance_surface_fully_4.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_4 <- rbind(distances_st_4, distances_ft_4)

distances_st_5.5 <- read.csv("./Distance_to_surface/semi-tryptic/flip_distance_surface_5_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_ft_5.5 <- read.csv("./Distance_to_surface/fully-tryptic/flip_distance_surface_fully_5_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_5.5 <- rbind(distances_st_5.5, distances_ft_5.5)

distances_st_7 <- read.csv("./Distance_to_surface/semi-tryptic/flip_distance_surface_7.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_ft_7 <- read.csv("./Distance_to_surface/fully-tryptic/flip_distance_surface_fully_7.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_7 <- rbind(distances_st_7, distances_ft_7)

distances_st_8.5 <- read.csv("./Distance_to_surface/semi-tryptic/flip_distance_surface_8_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_ft_8.5 <- read.csv("./Distance_to_surface/fully-tryptic/flip_distance_surface_fully_8_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_8.5 <- rbind(distances_st_8.5, distances_ft_8.5)

distances_st_10 <- read.csv("./Distance_to_surface/semi-tryptic/flip_distance_surface_10.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_ft_10 <- read.csv("./Distance_to_surface/fully-tryptic/flip_distance_surface_fully_10.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_10 <- rbind(distances_st_10, distances_ft_10)

FLiP_lib <- read.csv('./OutputData/FLiP_PBI_Library_anova_results.csv', header = TRUE, stringsAsFactors = FALSE, sep = ';')
FLiP_lib$sign <- sapply(FLiP_lib$q, function(x) return(ifelse(x < 0.05, TRUE, FALSE)))
FLiP_lib <- merge(FLiP_lib[, c("Peptide", "tryptic", "sign")], distances_1.5[, c("Peptide", "DistanceSurface.CA")])
FLiP_lib <- merge(FLiP_lib, distances_4[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
FLiP_lib <- merge(FLiP_lib, distances_5.5[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
FLiP_lib <- merge(FLiP_lib, distances_7[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
FLiP_lib <- merge(FLiP_lib, distances_8.5[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
FLiP_lib <- merge(FLiP_lib, distances_10[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")

colnames(FLiP_lib) <- c("Peptide", "tryptic", "sign", "1.5", "4", "5.5", "7", "8.5", "10")

FLiP_lib <- pivot_longer(FLiP_lib, cols = c("1.5", "4", "5.5", "7", "8.5", "10"), names_to = "Radius", values_to = "DistanceSurface.CA")
FLiP_lib$Radius <- factor(FLiP_lib$Radius, levels = c("1.5", "4", "5.5", "7", "8.5", "10"))

my_comparisons <- list(c("TRUE", "FALSE"))

ggplot(FLiP_lib, aes(x = Radius, y = DistanceSurface.CA, color = sign, fill = sign)) +  xlab("Radius (Å)") + scale_color_manual(values = c("#2680FF", "#FF8000")) + scale_fill_manual(values = c("#2680FF", "#FF8000")) +
  geom_boxplot(alpha = 0.4, outlier.colour="black", outlier.shape=16, outlier.size=1, notch=FALSE, show.legend = FALSE, facet.by = "Radius") + theme_classic(base_size = 22) + ylab("Distance to Surface (Å)") +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), method = "t.test")
ggsave("Plots/Distance_to_Surface_FLiP.pdf", device = "pdf", width = 25 , height = 12, units = "cm")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# WT HU dataset: are FLiP changes located closer to the surface than LiP changes?
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
distances_wt_st <- read.csv("./Distance_to_surface/semi-tryptic/wt_distance_surface.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_4_st <- read.csv("./Distance_to_surface/semi-tryptic/wt_distance_surface_4.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_5_5_st <- read.csv("./Distance_to_surface/semi-tryptic/wt_distance_surface_5_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_7_st <- read.csv("./Distance_to_surface/semi-tryptic/wt_distance_surface_7.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_8_5_st <- read.csv("./Distance_to_surface/semi-tryptic/wt_distance_surface_8_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_10_st <- read.csv("./Distance_to_surface/semi-tryptic/wt_distance_surface_10.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")

distances_wt_ft <- read.csv("./Distance_to_surface/fully-tryptic/wt_distance_surface_fully.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_4_ft <- read.csv("./Distance_to_surface/fully-tryptic/wt_distance_surface_fully_4.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_5_5_ft <- read.csv("./Distance_to_surface/fully-tryptic/wt_distance_surface_fully_5_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_7_ft <- read.csv("./Distance_to_surface/fully-tryptic/wt_distance_surface_fully_7.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_8_5_ft <- read.csv("./Distance_to_surface/fully-tryptic/wt_distance_surface_fully_8_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_wt_10_ft <- read.csv("./Distance_to_surface/fully-tryptic/wt_distance_surface_fully_10.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")

distances_wt <- rbind(distances_wt_ft, distances_wt_st)
distances_wt_4 <- rbind(distances_wt_4_ft, distances_wt_4_st)
distances_wt_5.5 <- rbind(distances_wt_5_5_ft, distances_wt_5_5_st)
distances_wt_7 <- rbind(distances_wt_7_ft, distances_wt_7_st)
distances_wt_8.5 <- rbind(distances_wt_8_5_ft, distances_wt_8_5_st)
distances_wt_10 <- rbind(distances_wt_10_ft, distances_wt_10_st)

WT <- read.csv("../2_LiP_WT_control_HU/OutputData/LiP_WT_control_HU_anova_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ';')
WT$sign <- rep(FALSE, nrow(WT))
WT[which(WT$q > 0.05 & WT$tryptic == "Specific"), "sign"] <- "Fully"
WT[which(WT$q > 0.05 & WT$tryptic != "Specific"), "sign"] <- "Semi"
WT[which(WT$q <= 0.05 & WT$tryptic != "Specific"), "sign"] <- "LiP Semi"
WT[which(WT$q <= 0.05 & WT$has_FLiP == TRUE  & WT$tryptic != "Specific"), "sign"] <- "FLiP Semi"
WT[which(WT$q <= 0.05 & WT$tryptic == "Specific"), "sign"] <- "LiP Fully"
WT[which(WT$q <= 0.05 & WT$has_FLiP == TRUE  & WT$tryptic == "Specific"), "sign"] <- "FLiP Fully"

WT <- merge(WT[, c("Peptide", "sign")], distances_wt[, c("Peptide", "DistanceSurface.CA")])
WT <- merge(WT, distances_wt_4[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
WT <- merge(WT, distances_wt_5.5[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
WT <- merge(WT, distances_wt_7[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
WT <- merge(WT, distances_wt_8.5[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
WT <- merge(WT, distances_wt_10[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")

colnames(WT) <- c("Peptide", "sign", "1.5", "4", "5.5", "7", "8.5", "10")

WT <- pivot_longer(WT, cols = c("1.5", "4", "5.5", "7", "8.5", "10"), names_to = "Radius", values_to = "DistanceSurface.CA")
WT$Radius <- factor(WT$Radius, levels = c("1.5", "4", "5.5","7", "8.5", "10"))
WT$sign <- factor(WT$sign,  levels= c("Semi", "Fully", "LiP Semi", "LiP Fully", "FLiP Semi", "FLiP Fully"))


ggplot(WT, aes(x = Radius, y= DistanceSurface.CA, color = sign, fill = sign)) +  xlab("Radius (Å)") + scale_color_manual(values = c("lightgrey", "grey", "lightblue", "orange",  "#2680FF", "#FF8000")) + 
  scale_fill_manual(values = c("lightgrey", "grey", "lightblue", "orange",  "#2680FF", "#FF8000")) +
  geom_boxplot(alpha = 0.4, outlier.colour="black", outlier.shape=16, outlier.size=1, notch=FALSE, facet.by = "Radius") + theme_classic(base_size = 22) + ylab("Distance to Surface (Å)") +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), method = "anova")
ggsave("Plots/Distance_to_Surface_FLiP_vs_LiP_wt.pdf", device = "pdf", width = 25 , height = 12, units = "cm")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Mutant HU dataset: are FLiP changes located closer to the surface than LiP changes?
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
distances_m_st <- read.csv("./Distance_to_surface/semi-tryptic/mut_distance_surface.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_4_st <- read.csv("./Distance_to_surface/semi-tryptic/mut_distance_surface_4.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_5_5_st <- read.csv("./Distance_to_surface/semi-tryptic/mut_distance_surface_5_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_7_st <- read.csv("./Distance_to_surface/semi-tryptic/mut_distance_surface_7.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_8_5_st <- read.csv("./Distance_to_surface/semi-tryptic/mut_distance_surface_8_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_10_st <- read.csv("./Distance_to_surface/semi-tryptic/mut_distance_surface_10.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")

distances_m_ft <- read.csv("./Distance_to_surface/fully-tryptic/mut_distance_surface_fully.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_4_ft <- read.csv("./Distance_to_surface/fully-tryptic/mut_distance_surface_fully_4.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_5_5_ft <- read.csv("./Distance_to_surface/fully-tryptic/mut_distance_surface_fully_5_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_7_ft <- read.csv("./Distance_to_surface/fully-tryptic/mut_distance_surface_fully_7.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_8_5_ft <- read.csv("./Distance_to_surface/fully-tryptic/mut_distance_surface_fully_8_5.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
distances_m_10_ft <- read.csv("./Distance_to_surface/fully-tryptic/mut_distance_surface_fully_10.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")

distances_m <- rbind(distances_m_ft, distances_m_st)
distances_m_4 <- rbind(distances_m_4_ft, distances_m_4_st)
distances_m_5.5 <- rbind(distances_m_5_5_ft, distances_m_5_5_st)
distances_m_7 <- rbind(distances_m_7_ft, distances_m_7_st)
distances_m_8.5 <- rbind(distances_m_8_5_ft, distances_m_8_5_st)
distances_m_10 <- rbind(distances_m_10_ft, distances_m_10_st)

M <- read.csv("../4_LiP_Gcn5_Mutant_control_HU/OutputData/LiP_Gcn5_Mutant_control_HU_anova_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ';')
M$sign <- rep(FALSE, nrow(M))
M[which(M$q > 0.05 & M$tryptic == "Specific"), "sign"] <- "Fully"
M[which(M$q > 0.05 & M$tryptic != "Specific"), "sign"] <- "Semi"
M[which(M$q <= 0.05 & M$tryptic != "Specific"), "sign"] <- "LiP Semi"
M[which(M$q <= 0.05 & M$has_FLiP == TRUE  & M$tryptic != "Specific"), "sign"] <- "FLiP Semi"
M[which(M$q <= 0.05 & M$tryptic == "Specific"), "sign"] <- "LiP Fully"
M[which(M$q <= 0.05 & M$has_FLiP == TRUE  & M$tryptic == "Specific"), "sign"] <- "FLiP Fully"

M <- merge(M[, c("Peptide", "sign")], distances_m[, c("Peptide", "DistanceSurface.CA")])
M <- merge(M, distances_m_4[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
M <- merge(M, distances_m_5.5[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
M <- merge(M, distances_m_7[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
M <- merge(M, distances_m_8.5[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")
M <- merge(M, distances_m_10[, c("Peptide", "DistanceSurface.CA")], by.x = "Peptide", by.y = "Peptide")

colnames(M) <- c("Peptide", "sign", "1.5", "4", "5.5", "7", "8.5", "10")

M <- pivot_longer(M, cols = c("1.5", "4", "5.5", "7", "8.5", "10"), names_to = "Radius", values_to = "DistanceSurface.CA")
M$Radius <- factor(M$Radius, levels = c("1.5", "4", "5.5","7", "8.5", "10"))
M$sign <- factor(M$sign,  levels= c("Semi", "Fully", "LiP Semi", "LiP Fully", "FLiP Semi", "FLiP Fully"))

ggplot(M, aes(x = Radius, y= DistanceSurface.CA, color = sign, fill = sign)) +  xlab("Radius (Å)") + scale_color_manual(values = c("lightgrey", "grey", "lightblue", "orange",  "#2680FF", "#FF8000")) + 
  scale_fill_manual(values = c("lightgrey", "grey", "lightblue", "orange",  "#2680FF", "#FF8000")) +
  geom_boxplot(alpha = 0.4, outlier.colour="black", outlier.shape=16, outlier.size=1, notch=FALSE, facet.by = "Radius") + theme_classic(base_size = 22) + ylab("Distance to Surface (Å)") +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), method = "anova")
ggsave("Plots/Distance_to_Surface_FLiP_vs_LiP_M.pdf", device = "pdf", width = 25 , height = 12, units = "cm")