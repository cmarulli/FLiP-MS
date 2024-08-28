library(hash)
source("../../Rfunctions/barplot_numbers_ggplot.R")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Read in the data
#-----------------------------------------------------------------------------------------------------------------------------------------------------
# read in FLiP library
FLiP_lib <- read.csv("./final_FLiP_lib/FiLiP_lib_marker_confidence.csv", header = TRUE, stringsAsFactors = FALSE)
# read in complex subunits
complex_portal <- read.delim('../Files/complex_portal_559292.tsv', header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# read in UniProt "Mass" annotation
uniprot <- read.delim("../Files/uniprotkb_proteome_UP000002311_2023_12_20.tsv", header = TRUE, stringsAsFactors = FALSE, sep = '\t')

# read in protein intensities detected in a full lysate tryptic control
TC <- read.delim("../2_LiP_WT_control_HU/FilteredData/LiP_WT_control_HU_tc_data_filtered.tsv", header = TRUE, stringsAsFactors = FALSE, sep = ';')
TC <- TC[-grep("iRT", TC$Protein), ]
TC$Mean <- apply(TC[, c(2:5)], 1, function(x) return(mean(x, na.rm = TRUE)))
TC$MW <- sapply(TC$Protein, function(x) return(uniprot[which(uniprot$Entry == x), "Mass"]))


#-----------------------------------------------------------------------------------------------------------------------------------------------------
# plot MW density of proteins detected in a yeast lysate 
#-----------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(TC, aes(x=MW, weights = Mean)) +  geom_density(color = "#0057b7", size = 1.5) + theme_classic( base_size = 26) + xlim(0, 4e5) + xlab("Molecular Weight (kDa)") + ylab("Density")
ggsave("Plots/Density_MolecularWeigth_detected_proteins_weighted.pdf", device = "pdf", width = 26 , height = 16, units = "cm")

#-----------------------------------------------------------------------------------------------------------------------------------------------------
# Plot number of protein complex subunits for proteins with and without FLiP markers
#-----------------------------------------------------------------------------------------------------------------------------------------------------
FLiP_proteins <- unique(FLiP_lib$Protein)

complex_subunits_FLiP <- c()
complex_subunits_all <- c()

for(i in c(1:nrow(complex_portal))){
  
  proteins <- unlist(strsplit(complex_portal[i, "Identifiers..and.stoichiometry..of.molecules.in.complex"], "\\|"))
  if(length(grep("C", proteins))> 0){
    proteins <- proteins[-grep("C", proteins)]
  }
  stoich <- sapply(proteins, function(x) return(strsplit(x, "\\(")[[1]][2]))
  stoich <- sapply(stoich, function(x) return(strsplit(x, "\\)")[[1]][1]))
  
  # if stoichiometry is 0 assign 1, because this means that it is simply not known
  stoich[stoich==0] <- 1
  proteins <- sapply(proteins, function(x) return(strsplit(x, "\\(")[[1]][1]))
  
  if(length(intersect(proteins, FLiP_proteins) > 0)){
    complex_subunits_FLiP <- c(complex_subunits_FLiP, sum(as.numeric(stoich)))
  }
  
  complex_subunits_all <- c(complex_subunits_all, sum(as.numeric(stoich)))
}

complex_portal$subunits <- complex_subunits_all
complex_portal$has_FLiP <- rep(FALSE, nrow(complex_portal)) 
complex_portal[grep(paste(FLiP_proteins, collapse = "|"), complex_portal$Expanded.participant.list), "has_FLiP"] <- TRUE

barplot_numbers_ggplot(names(table(complex_subunits_FLiP)[-1]), table(complex_subunits_FLiP)[-1], file_location = "./Plots/FLiP_Marker_Complex_Subunits", y_lim = 78, my_width = 25, my_basesize = 16, 
                       my_geom_text_size = 6, xlab = "Number of Complex Subunits", ylab = "Number of Complexes with FLiP Marker")

barplot_numbers_ggplot(names(table(complex_subunits_all)[-1]), table(complex_subunits_all)[-1], file_location = "./Plots/All_Complex_Subunits", y_lim = 300, my_width = 25, my_basesize = 16, 
                       my_geom_text_size = 6, xlab = "Number of Complex Subunits", ylab = "Number of Complexes with FLiP Marker")

ggplot(complex_portal, aes(x = subunits, fill = has_FLiP, color = has_FLiP)) + geom_histogram(binwidth = 1) + theme_classic() + 
  scale_color_manual(values=c("grey", "#0057b7")) + scale_fill_manual(values=c("grey", "#0057b7")) + 
  xlab("Number of Subunits") + ylab("Number of Protein Complexes") + xlim(0, 60)
ggsave("./Plots/Number_of_subunits.pdf", device = "pdf", width = 20 , height = 14, units = "cm", useDingbats=FALSE)