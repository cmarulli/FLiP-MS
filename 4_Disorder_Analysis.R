library(data.table)

# color palette
colfunc<-colorRampPalette(c("#b0c7e0","#0057b7","#001e64"))

# disorder information from UniProt
disorder_annotation <- read.delim("../Files/uniprot-proteome_UP000002311+reviewed_yes.tab", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
disorder_annotation <- disorder_annotation[grep("Disorder", disorder_annotation$Region), ]

# extract position information about disorder
disorder_annotation$Position <- sapply(disorder_annotation$Region, function(x){
  
  a <- strsplit(x, split = ";")[[1]]
  a <- a[grep("REGION", a)]
  a <- as.vector(unlist(sapply(a, function(x) return(strsplit(x, " ")[[1]]))))
  a <- a[grep("\\..", a)]
  return(paste(a, collapse = ";"))
  
  } 
)

disorder_annotation <- setDT(disorder_annotation)[, .(Position = as.vector(unlist(strsplit(Position, ";")))), by = "Entry"]
disorder_annotation$Start <- unlist(sapply(disorder_annotation$Position, function(x) strsplit(x, "\\..")[[1]][1]))
disorder_annotation$End <- unlist(sapply(disorder_annotation$Position, function(x) strsplit(x, "\\..")[[1]][2]))
disorder_annotation <- disorder_annotation[, c(1,3,4)]
colnames(disorder_annotation) <-  c("Protein", "Disorder_Start","Disorder_End")
write.csv(disorder_annotation, "../Files/disorder_annotation.csv", row.names = FALSE, quote = FALSE)

# read in FLiP library and add end position information
FLiP_lib <- read.delim("./OutputData/FLiP_PBI_Library_anova_results.csv", header = TRUE, stringsAsFactors = FALSE, sep =";")
FLiP_lib <- FLiP_lib[, c("Protein", "Peptide", "q", "tryptic", "L2FC", "Position", "fully_tryptic")]
FLiP_lib$Position_End <- apply(FLiP_lib,1, function(x) return(as.numeric(x[6]) + nchar(x[2])))
FLiP_lib_merged <- merge(FLiP_lib, disorder_annotation)              

# plot the peptide length distribution of all peptides in the FLiP Library - separate by semi, fully, all
FLiP_lib_length <- FLiP_lib[, c("Peptide", "tryp tic")]
FLiP_lib_length$length <- sapply(FLiP_lib_length$Peptide, function(x) return(nchar(x)))

library("spatstat")
jpeg(paste("./Plots/Peptides_length_distribution.jpeg", sep =""),  width = 20, height = 20,units = "cm", res = 2000)
plot(density(FLiP_lib_length$length), xlim = c(0,30), ylim = c(0,0.1), col =  "#ffc900", lwd = 5, 
     xlab = "", main = "", ylab = "", cex.axis = 2.5)
abline(v = median(FLiP_lib_length$length), col = "#ffc900", lwd = 5, lty=2)
lines(density(FLiP_lib_length[which(FLiP_lib_length$tryptic == "Specific"), "length"]), col = "#a7171a", lwd = 5)
abline(v = median(FLiP_lib_length[which(FLiP_lib_length$tryptic == "Specific"), "length"]), col = "#a7171a", lwd = 5, lty=2)
lines(density(FLiP_lib_length[which(FLiP_lib_length$tryptic != "Specific"), "length"]), col = "#0057b7", lwd = 5)
abline(v = median(FLiP_lib_length[which(FLiP_lib_length$tryptic != "Specific"), "length"]), col = "#0057b7", lwd = 5, lty=2)
dev.off()


# Keep only the disorder annotations for which the location in the protein overlaps with the detected peptide 
peptides <- unique(FLiP_lib_merged$Peptide)
disordered_peptides <- c()
for(i in c(1:length(peptides))){
  
  current_peptide <- FLiP_lib_merged[which(FLiP_lib_merged$Peptide == peptides[i]), c("Position", "Position_End", "Disorder_Start", "Disorder_End")]
  intersect <- as.vector(apply(current_peptide, 1, function(x) length(intersect(c(x[1]:x[2]), c(x[3]:x[4])))))
  
  if(sum(intersect) > 1){
    disordered_peptides <- c(disordered_peptides, peptides[i])
  }
}

FLiP_lib$disordered <- rep(FALSE, nrow(FLiP_lib))
FLiP_lib[which(FLiP_lib$Peptide %in% disordered_peptides), "disordered"] <- TRUE

FLiP_lib$sign <- rep(FALSE, nrow(FLiP_lib))
FLiP_lib[which(FLiP_lib$q < 0.05), "sign"] <- TRUE

# test whether FLiP marker peptides are more disordered than non-marker peptides
p_value_FLiP_markers_disordered<- fisher.test(FLiP_lib$disordered, FLiP_lib$sign, alternative = "greater")

percentage_disorder_FLiP <- length(which(FLiP_lib$disordered == TRUE))/nrow(FLiP_lib)*100
percentage_disorder_FLiP_not_sign <- length(which(FLiP_lib$disordered == TRUE & FLiP_lib$sign == FALSE))/length(which(FLiP_lib$sign == FALSE))*100
percentage_disorder_FLiP_sign <- length(which(FLiP_lib$disordered == TRUE & FLiP_lib$sign == TRUE))/length(which(FLiP_lib$sign == TRUE))*100

df <- data.frame(name = c(c("All", "Significant", "Not Signficant", "All", "Significant", "Not Signficant")), 
                 percentage = c(100-percentage_disorder_FLiP, 100-percentage_disorder_FLiP_sign, 100-percentage_disorder_FLiP_not_sign, percentage_disorder_FLiP, percentage_disorder_FLiP_sign, percentage_disorder_FLiP_not_sign))
df$Disorder <- c("Not Disordered", "Not Disordered", "Not Disordered", "Disordered", "Disordered" ,"Disordered")
df$Disorder <- as.factor(df$Disorder)
df$Disorder <- relevel(df$Disorder, 'Not Disordered')

ggplot(df, aes(x = name, y = percentage, fill = Disorder)) + geom_bar( stat="identity", show.legend = TRUE, ) + ylim(c(0, 100)) +theme_classic(base_size = 22) + xlab("") +
  ylab ("% Disordered Peptides") + scale_fill_manual(values = c("grey75", "#0057b7")) + geom_text(aes(label = round(percentage, 1)), vjust = -0.5, size = 8)

ggsave("./Plots/FLiP_Disorder.pdf", device = "pdf", width = 20 , height = 16, units = "cm", useDingbats=FALSE)

barplot_numbers_ggplot(c("All", "Significant", "Not Signficant"), number = c(percentage_disorder_FLiP, percentage_disorder_FLiP_sign, percentage_disorder_FLiP_not_sign), 
                       file_location = "./Plots/FLiP_Disorder.jpg", y_lim = c(0, 100))
