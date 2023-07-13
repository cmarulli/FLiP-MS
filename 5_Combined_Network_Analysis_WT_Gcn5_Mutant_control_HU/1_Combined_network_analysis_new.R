library('data.table')
library('ggplot2')
library("RColorBrewer")
library("plyr")
library("dplyr")
library("stringr")
library("purrr")
library("igraph")
library("MCL")
library("protti")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# source all functions
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
files.sources = list.files("functions", full.names=TRUE, pattern="*.R$")
sapply(files.sources, source)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read in information from wt and mutant datasets
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LiP changes
LiP_changes_WT <- read.csv("../2_LiP_WT_control_HU/OutputData/LiP_WT_control_HU_anova_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
LiP_changes_M <-  read.csv("../4_LiP_Gcn5_Mutant_control_HU/OutputData/LiP_Gcn5_Mutant_control_HU_anova_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")

# exclude peptides that cannot be located to specific position on proteins
LiP_changes_WT <- LiP_changes_WT[-which(grepl(',', LiP_changes_WT$Position)), ]
LiP_changes_M <- LiP_changes_M[-which(grepl(',', LiP_changes_M$Position)), ]

# only keep the peptides in protein regions covered in both datasets

# which peptides in the WT dataset overlap with peptides in the mutant dataset
overlapping_WT <- c()
for(i in c(1:nrow(LiP_changes_WT))){
  
  pep_pos_WT <- c(LiP_changes_WT[i, "Position"]:LiP_changes_WT[i, "Stop_Position"])
  pep_pos_M <- LiP_changes_M[which(LiP_changes_M$Protein == LiP_changes_WT[i, "Protein"]), ]
  pep_pos_M <- apply(pep_pos_M, 1, function(x) return(c(x[12]:x[15])))
  
  intersect <- unlist(lapply(pep_pos_M, function(x) return(length(intersect(pep_pos_WT, x)) > 0)))
  
  if(TRUE %in% intersect){
    overlapping_WT <- c(overlapping_WT, TRUE)
  }else{
    overlapping_WT <- c(overlapping_WT, FALSE)
  }
}

# which peptides in the M dataset overlap with peptides in the WT dataset
overlapping_M <- c()
for(i in c(1:nrow(LiP_changes_M))){
  
  pep_pos_M <- c(LiP_changes_M[i, "Position"]:LiP_changes_M[i, "Stop_Position"])
  pep_pos_WT <- LiP_changes_WT[which(LiP_changes_WT$Protein == LiP_changes_M[i, "Protein"]), ]
  pep_pos_WT <- apply(pep_pos_WT, 1, function(x) return(c(x[12]:x[15])))
  
  intersect <- unlist(lapply(pep_pos_WT, function(x) return(length(intersect(pep_pos_M, x)) > 0)))
  
  if(TRUE %in% intersect){
    overlapping_M <- c(overlapping_M, TRUE)
  }else{
    overlapping_M <- c(overlapping_M, FALSE)
  }
}

# keep only peptides for which the corresponding region is also covered in the other dataset
LiP_changes_WT <- LiP_changes_WT[which(overlapping_WT == TRUE), ]
LiP_changes_M <- LiP_changes_M[which(overlapping_M == TRUE), ]

FiLiP_changes_WT <- LiP_changes_WT[which(LiP_changes_WT$has_FLiP == TRUE), ]
FiLiP_changes_M <- LiP_changes_M[which(LiP_changes_M$has_FLiP == TRUE), ]

detected_WT <- unique(LiP_changes_WT[ , "Protein"])
detected_M <- unique(LiP_changes_M[ , "Protein"])

# store all detected proteins
detected_proteins <- c(unique(LiP_changes_WT$Protein, LiP_changes_M$Protein))

# subset to significant proteins
LiP_changes_WT <- unique(LiP_changes_WT[which(LiP_changes_WT$q < 0.05), "Protein"])
LiP_changes_M <- unique(LiP_changes_M[which(LiP_changes_M$q < 0.05), "Protein"])

# subset to significant high confidence marker proteins 
high_confidence_markers_WT <- unique(FiLiP_changes_WT[which(FiLiP_changes_WT$q < 0.05 & FiLiP_changes_WT$confidence =='high'), "Protein"])
high_confidence_markers_M <- unique(FiLiP_changes_M[which(FiLiP_changes_M$q < 0.05 & FiLiP_changes_M$confidence =='high'), "Protein"])

# subset to significant proteins
FiLiP_changes_WT <- unique(FiLiP_changes_WT[which(FiLiP_changes_WT$q < 0.05), "Protein"])
FiLiP_changes_M <- unique(FiLiP_changes_M[which(FiLiP_changes_M$q < 0.05), "Protein"])

# Abundance changes
abundance_changes_WT <- read.csv("../2_LiP_WT_control_HU/OutputData/LiP_WT_control_HU_Protein_Abundance_Changes.csv", header = TRUE, stringsAsFactors = FALSE)
abundance_changes_M <- read.csv("../4_LiP_Gcn5_Mutant_control_HU/OutputData/LiP_Gcn5_Mutant_control_HU_Protein_Abundance_Changes.csv", header = TRUE, stringsAsFactors = FALSE)

# only keep the proteins regions covered in both datasets
abundance_changes_WT <- abundance_changes_WT[which(abundance_changes_WT$Protein %in% abundance_changes_M$Protein), ]
abundance_changes_M <- abundance_changes_M[which(abundance_changes_M$Protein %in% abundance_changes_WT$Protein), ]

# subset to significant proteins
abundance_changes_WT <- unique(abundance_changes_WT[which(abundance_changes_WT$q < 0.05 & abs(abundance_changes_WT$L2FC) > 0.5), "Protein"])
abundance_changes_M <- unique(abundance_changes_M[which(abundance_changes_M$q < 0.05 & abs(abundance_changes_M$L2FC) > 0.5), "Protein"])

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read in information about SAGA and Gcn5 targets
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SAGA_targets <- read.csv("../Files/SAGA_transcription_targets.txt", header = TRUE, sep = '\t')
gene_names <- read.csv("../Files/210713_Yeast_Uniprot_Names.csv", header = TRUE, stringsAsFactors = FALSE)
gene_names$Gene.names <- sapply(gene_names$Gene.names, function(x) return(strsplit(x, " ")[[1]][1]))
SAGA_targets$Protein <- sapply(SAGA_targets$Gene.Name, function(x) return(gene_names[grep(paste('^', x, '$', sep = ""), gene_names$Gene.names), "Entry"]))
SAGA_targets <- unlist(unique(SAGA_targets$Protein))

# https://doi.org/10.1074/mcp.M114.043141
Gcn5_acetylation_targets <- read.csv("../Files/GCN5_acetylation_targets.csv", header = TRUE)
Gcn5_acetylation_targets <- unlist(unique(sapply(c(Gcn5_acetylation_targets$GENE, "NUP60", "UME6"), function(x) return(gene_names[grep(paste('^', x, '$', sep = ""), gene_names$Gene.names), "Entry"]))))

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 0) house keeping: create empty folders for output
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings
delete_existing_files <- TRUE
initialize_folders(delete_existing_files)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1) filter input network
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings
network <- "complexPortal"     
network_output_filename <-  "network_interactions"

# complexPortal
complexPortal_input_file <- "../Files/complex_portal_559292.tsv"
complexPortal_mapping_table_file <- "../Files/uniprot-proteome_mapping.tab"
complexPortal_output_filename <-  "complex_portal_interactions"

filter_complexPortal_network(complexPortal_input_file, complexPortal_mapping_table_file,
                             network_output_filename)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2) create i-graph objects
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings
network_edgelist_file <- paste("intermediate/", network_output_filename, ".tsv", sep = "" )

# run the function
create_igraph_objects(network_edgelist_file, 
                      FiLiP_changes_WT, FiLiP_changes_M, LiP_changes_WT, LiP_changes_M,  
                      abundance_changes_WT, abundance_changes_WT, 
                      detected_WT, detected_M, SAGA_targets, Gcn5_acetylation_targets, 
                      high_confidence_markers_WT, high_confidence_markers_M)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3) network propagation
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
run_pageRank(dampening_factor = 0.9, pageRank_percentile_cutoff = 0.6, FiLiP_changes = FiLiP_changes_WT, filename = "WT")
run_pageRank(dampening_factor = 0.9, pageRank_percentile_cutoff = 0.6, FiLiP_changes = FiLiP_changes_M, filename = "M")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4) Clustering
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings walktrap
walktrap_clustering(step_number  = 6)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5) Filter Clusters and Plot
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cluster_min_size <- 4
cluster_max_size <- 500
cluster_col_pal = "Dark2"
complex_portal <- fread("../Files/complex_portal_559292.tsv")

# run the function
cluster_analysis(cluster_min_size, cluster_max_size,
                 cluster_col_pal, cluster_q_value_cutoff, complex_portal)
