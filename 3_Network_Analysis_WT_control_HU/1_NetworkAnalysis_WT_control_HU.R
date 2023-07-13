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

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# source all functions
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
files.sources = list.files("functions", full.names=TRUE, pattern="*.R$")
sapply(files.sources, source)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 0) house keeping: create empty folders for output
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
delete_existing_files <- TRUE
initialize_folders(delete_existing_files)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1) filter input network
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings
network <- "complexPortal"     
network_output_filename <-  "network_interactions"

# complexPortal
complexPortal_input_file <- "../Files/complex_portal_559292.tsv"
complexPortal_mapping_table_file <- "../Files/uniprot-proteome_mapping.tab"
complexPortal_output_filename <-  "complex_portal_interactions"

#    run the function
if (network == "string" ){
  filter_string_network(string_score_to_use, string_score_cutoff, 
                        string_input_file, string_mapping_table_file,
                        network_output_filename)
} else if (network =="complexPortal"){
  filter_complexPortal_network(complexPortal_input_file, complexPortal_mapping_table_file,
                               network_output_filename)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2) create i-graph objects
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings
p_cutoff <- 0.05
filter_column <- "q"  
add_non_connected <- TRUE
network_edgelist_file <- paste("intermediate/", network_output_filename, ".tsv", sep = "" )

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3) read in experimental and annotation data
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
LiP_changes <- read.csv("../2_LiP_WT_control_HU/OutputData/LiP_WT_control_HU_anova_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";")
FiLiP_changes <- LiP_changes[which(LiP_changes$has_FLiP == TRUE), ]
# detected proteins
detected <- unique(LiP_changes[ , "Protein"])

LiP_changes <- unique(LiP_changes[which(LiP_changes$q < 0.05 & LiP_changes$has_FLiP == FALSE), "Protein"])
high_confidence_markers <- unique(FiLiP_changes[which(FiLiP_changes$q < 0.05 & FiLiP_changes$confidence =='high'), "Protein"])
FiLiP_changes <- unique(FiLiP_changes[which(FiLiP_changes$q < 0.05), "Protein"])

abundance_changes <- read.csv("../2_LiP_WT_control_HU/OutputData/LiP_WT_control_HU_Protein_Abundance_Changes.csv", header = TRUE, stringsAsFactors = FALSE)
abundance_changes <- unique(abundance_changes[which(abundance_changes$q < 0.05 & abs(abundance_changes$L2FC) > 0.5 ), "Protein"])

# read in annotation about SAGA transcription targets and Gcn5 acetylation targets
SAGA_targets <- read.csv("../Files/SAGA_transcription_targets.txt", header = TRUE, sep = '\t')
gene_names <- read.csv("../Files/210713_Yeast_Uniprot_Names.csv", header = TRUE, stringsAsFactors = FALSE)
gene_names$Gene.names <- sapply(gene_names$Gene.names, function(x) return(strsplit(x, " ")[[1]][1]))
SAGA_targets$Protein <- sapply(SAGA_targets$Gene.Name, function(x) return(gene_names[grep(paste('^', x, '$', sep = ""), gene_names$Gene.names), "Entry"]))
SAGA_targets <- unlist(unique(SAGA_targets$Protein))

# https://doi.org/10.1074/mcp.M114.043141
Gcn5_acetylation_targets <- read.csv("../Files/GCN5_acetylation_targets.csv", header = TRUE)
Gcn5_acetylation_targets <- unlist(unique(sapply(c(Gcn5_acetylation_targets$GENE, "NUP60", "UME6"), function(x) return(gene_names[grep(paste('^', x, '$', sep = ""), gene_names$Gene.names), "Entry"]))))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4) create the annotated network
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
create_igraph_objects(p_cutoff, filter_column, add_non_connected,
                      network_edgelist_file, detected, FiLiP_changes, LiP_changes, abundance_changes, 
                      SAGA_targets, Gcn5_acetylation_targets, high_confidence_markers)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 5) Netowork propagation using pagerank
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings
dampening_factor <- 0.90
pageRank_percentile_cutoff <- 0.60

# run the function
run_pageRank(dampening_factor, pageRank_percentile_cutoff)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 6) Clustering of the network
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# settings walktrap
cluster_method <- "walktrap"
walktrap_step_number <- 4
walktrap_clustering(walktrap_step_number)

cluster_min_size <- 4
cluster_max_size <- 500
cluster_col_pal = "Dark2"
cluster_q_value_cutoff <- 1

# read in information about the protein complex names from the complex portal database
complex_portal <- fread("../Files/complex_portal_559292.tsv")

# run the function
cluster_analysis(cluster_method, cluster_min_size, cluster_max_size,
                 cluster_col_pal, cluster_q_value_cutoff, complex_portal, gene_names)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 7) Global statistics
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# run the function
summary_plots()
save.image(file='final/settings_used.RData')