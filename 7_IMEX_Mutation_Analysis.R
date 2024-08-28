library('data.table')
library('ggplot2')
library("RColorBrewer")
library("dplyr")
library("stringr")
library("purrr")
library('protti')

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### settings
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(123)
species <- "559292 - Saccharomyces cerevisiae"
sinlge_AA_mutation_only <- FALSE
exclude_protein_RNA <- FALSE
exclude_actin <- TRUE # extreme outlier in yeast with a lot of mutations
FLiP_file <- "../1_FLiP_PBI_Library/OutputData/FLiP_PBI_Library_anova_results.csv"
AA_tolerance <- 0 # mutation can be X AAs to the left or right of a peptide
peptide_or_parent <- "parent"  #c("peptide", "parent") 
number_of_control_runs <- 10 

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# read in an format data
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#### load data
mutation_data <- fread("../Files/mutations.tsv")
FLiP_library <- read.csv(FLiP_file, header = TRUE, stringsAsFactors = FALSE, sep = ";")

# reduce to proteins which have at least one marker
marker_proteins <- unique(FLiP_library[which(FLiP_library$q < 0.05), "Protein"])
FLiP_library <- FLiP_library[which(FLiP_library$Protein %in% marker_proteins), ]

# add protein length
FLiP_library$length <- sapply(FLiP_library$Sequence, function(x) return(nchar(x)+1))

#### filter for species
mutation_data <- mutation_data[mutation_data$`Affected protein organism` == species,]

#### get type of mutation & filter for single AA point mutations
orig_seq_len <- nchar(mutation_data$`Original sequence`)
new_seq_len <- nchar(mutation_data$`Resulting sequence`)

if (sinlge_AA_mutation_only){
  mutation_data <- mutation_data[orig_seq_len == 1 & new_seq_len == 1, ]
}

#### filter out non protein interactions
if (exclude_protein_RNA){
  index <-  which(  str_count( mutation_data$`Interaction participants`, "stable complex") >= 1|
                      str_count( mutation_data$`Interaction participants`, "protein") >= 2 |
                      str_count( mutation_data$`Interaction participants`, "protein") ==1 &  str_count( mutation_data$`Interaction participants`, "small molecule") == 0  & str_count(mutation_data$`Interaction participants`, "rna\\(MI") == 0)
  
} else {
  index <-  which(  str_count( mutation_data$`Interaction participants`, "stable complex") >= 1|
                      str_count( mutation_data$`Interaction participants`, "protein") >= 2 |
                      str_count( mutation_data$`Interaction participants`, "protein") ==1 &  str_count( mutation_data$`Interaction participants`, "small molecule") == 0)
}

mutation_data <- mutation_data[index, ]

#### parse mutant data columns a bit - get UniprotID, start and end position in separate columns
uniprot_ID <- character(nrow(mutation_data))
start <- numeric(nrow(mutation_data))
end <- numeric(nrow(mutation_data))

uniprot_list <- strsplit(mutation_data$`Feature short label`, ":")
seq_list <- strsplit(mutation_data$`Feature range(s)`, "-")

for (i in 1:nrow(mutation_data)){
  uniprot_ID[i] <- uniprot_list[[i]][1]
  start[i] <- seq_list[[i]][1]
  end[i] <- seq_list[[i]][2]
}

mutation_data <- cbind(mutation_data, uniprot_ID, start, end)
mutation_data$start <- as.numeric(mutation_data$start)
mutation_data$end <- as.numeric(mutation_data$end)

### filter mutations - are they in a detected protein?
mutation_protein <- unique(mutation_data$uniprot_ID)
mutation_data <- mutation_data[mutation_data$uniprot_ID %in% unique(FLiP_library$Protein), ]

# get start and end positions
if (peptide_or_parent == "peptide"){
  FLiP_library$Position <- as.numeric(FLiP_library$Position)
  FLiP_library$Position_end <- FLiP_library$Position + nchar(FLiP_library$Peptide)
  
} else if(peptide_or_parent == "parent"){
  Position <- numeric()
  i <- 1
  for (i in 1:nrow(FLiP_library)){
    current_pos <-  as.numeric(gregexpr(pattern = FLiP_library$fully_tryptic[i], text = FLiP_library$Sequence[i])[[1]])
    Position <- rbind(Position, current_pos)
  }
  
  FLiP_library$Position <- NULL
  FLiP_library <- cbind(FLiP_library, Position)
  FLiP_library$Position_end <- FLiP_library$Position + nchar(FLiP_library$fully_tryptic)
}

# get rid of mutations at the same position
mutation_data$positional_mutant <- paste(mutation_data$uniprot_ID, "_", as.character(mutation_data$start), sep = "")
mutation_data <- mutation_data[-which(duplicated(mutation_data$positional_mutant)), ]

# summarize mutation that are very close to each other as a single mutation (otherwise they are overrepresented)
mutation_data_reduced <- data.frame()

for(p in unique(mutation_data$uniprot_ID)){
  
  # get the mutation sites of the current protein
  data <- mutation_data[which(mutation_data$uniprot_ID == p), ]
  data <- data[order(data$start), ]
  # extend the start and end position by two
  data$start_ext <- sapply(data$start, function(x) return(x-2))
  data$end_ext <- sapply(data$end, function(x) return(x+2))
  
  # if the extented start and end positions overlap, assign the same index to the mutation sites
  if(nrow(data) > 2){
    
    current_index <- 1
    index <- c(current_index)
    
    for(i in c(1:(nrow(data) - 1))){
      
      if(length(intersect(c(as.numeric(data[i, "start_ext"]) : as.numeric(data[i, "end_ext"])), c(as.numeric(data[i+1, "start_ext"]) : as.numeric(data[i+1, "end_ext"])))) == 0){
        current_index <- current_index+1 
        index <- c(index, current_index)
      }else{
        index <- c(index, current_index)
      }
    }
    
    data$index <- index
    
    # reduce the data to unique indexes spanning multiple closely located mutation sites
    new_data <- data %>%  group_by(index) %>% summarize(start = min(start), end = max(end))
    new_data <- unique(merge(new_data, data[, c("index", "uniprot_ID")]))
    new_data <- new_data[, c("uniprot_ID", "start", "end")]
    
    mutation_data_reduced <- rbind(mutation_data_reduced, new_data)
    
  }else{
    new_data <- data[, c("uniprot_ID", "start", "end")]
    mutation_data_reduced <- rbind(mutation_data_reduced, new_data)
  }
}

mutation_data <- mutation_data_reduced

# store the current mutation data table as control
mutation_data_control <- mutation_data

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# For each mutation check whether it was detected in the FLiP dataset and store the lowest q value
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mutation_detected <- rep(FALSE, nrow(mutation_data)) # vector: is there a peptide overlapping with the mutation
lowest_FLiP_q <- rep(1, nrow(mutation_data)) # vector: lowest q values of the overlapping peptides
lowest_FLiP_q_peptide <- rep("", nrow(mutation_data))

FLiP_library_total <- FLiP_library # needed for sampling random proteins (further down)
FLiP_library <- FLiP_library[FLiP_library$Protein %in% mutation_protein, ]

for (m in 1:nrow(mutation_data)){
  FLiP_peptides <- FLiP_library[FLiP_library$Protein == mutation_data$uniprot_ID[m], ]
  mut_covered <- rep(FALSE, nrow(FLiP_peptides))
  
  for(p in 1:nrow(FLiP_peptides)){ # check each peptide of the protein if it overlaps
    start_overlap <- mutation_data$start[m] >= FLiP_peptides$Position[p] - AA_tolerance & mutation_data$start[m] <= FLiP_peptides$Position_end[p] + AA_tolerance # LHS
    end_overlap <- mutation_data$end[m] >= FLiP_peptides$Position[p] - AA_tolerance & mutation_data$end[m] <= FLiP_peptides$Position_end[p] + AA_tolerance       # RHS
    mut_covered[p] <- start_overlap+end_overlap > 0 # if the peptide overlaps to the left or right, it is overlapping
    if (sum(mut_covered) > 0){
      mutation_detected[m] <- TRUE
    } 
  }
  
  if (mutation_detected[m] == TRUE){ # if there are several peptides overlapping, take the one with the lowest q-value of all overlapping peptides to mark the mutation
    lowest_FLiP_q[m]<- min(FLiP_peptides$q[mut_covered]) 
    lowest_FLiP_q_peptide[m] <- FLiP_peptides$fully_tryptic[which(FLiP_peptides$q == lowest_FLiP_q[m])]
  }
}

mutation_data <- cbind(mutation_data, lowest_FLiP_q)
mutation_data <- cbind(mutation_data, lowest_FLiP_q_peptide)
mutation_data <- mutation_data[mutation_detected == TRUE, ] # reduce to mutations detected in dataset

# get rid of mutations that are covered by the same peptide
overlapping_peptides <- unique(mutation_data$lowest_FLiP_q_peptide)
mutation_data_reduced <- data.frame()
for(p in overlapping_peptides){
  mutation_data_reduced <- rbind(mutation_data_reduced, mutation_data[grep(p, mutation_data$lowest_FLiP_q_peptide)[1], ])
}
mutation_data <- mutation_data_reduced

export_total_number_mutations <- nrow(mutation_data)
export_actin_number_mutations <- nrow(mutation_data[which(mutation_data$uniprot_ID == "P60010"), ])

### remove actin due to its large number of mutations (extreme outlier)
if (exclude_actin & "P60010" %in% FLiP_library$Protein){
  mutation_data <- mutation_data[-which(mutation_data$uniprot_ID == "P60010"), ]
  FLiP_library <- FLiP_library[-which(FLiP_library$Protein  == "P60010"), ]
}

### calculate how many proteins have a FLiP marker covering them
percentage_M_with_FLiP <- sum(mutation_data$lowest_FLiP_q <= 0.05) / length(mutation_data$lowest_FLiP_q)

export_number_mutations_covered <- sum(mutation_data$lowest_FLiP_q <= 0.05)
export_number_mutations_not_covered <- sum(mutation_data$lowest_FLiP_q > 0.05)
export_share_covered <- percentage_M_with_FLiP
mutations_per_protein <- table(mutation_data$uniprot_ID)
write.table(mutations_per_protein, "./OutputData/Number_of_mutations_per_protein.csv", sep = ",", row.names = F)

plot(table(mutations_per_protein), col = "blue", xlab = "number of mutations per protein", ylab = "number of proteins", 
     ylim = c(0, max(table(mutations_per_protein)+5)),  cex.main=2, cex.lab=2,cex.axis=2 )

plot( table(FLiP_library$Protein), col = "blue", xlab = "number of peptides per protein", ylab = "number of proteins", 
     ylim = c(0, max(table(FLiP_library$Protein)+5)),  cex.main=2, cex.lab=2,cex.axis=2 )

# some plots
for(p in unique(mutation_data$uniprot_ID)){

  my_prot <- FLiP_library[which(FLiP_library$Protein == p & FLiP_library$q < 0.05), c("Protein", "Position", "Position_end", "length")]
  my_prot$q <- rep(0.001, nrow(my_prot))
  my_prot$L2FC <- rep(0, nrow(my_prot))

  my_mutation <- mutation_data[which(mutation_data$uniprot_ID == p), c("start", "end")]
  my_mutation$end <- sapply(my_mutation$end, function(x) return(x+1))
  colnames(my_mutation) <- c("Position", "Position_end")
  my_mutation$Protein <- rep(p, nrow(my_mutation))
  my_mutation$length <- rep(my_prot$length[1], nrow(my_mutation))
  my_mutation$q <- rep(1, nrow(my_mutation))
  my_mutation$L2FC <- rep(0, nrow(my_mutation))

  my_prot <- rbind(my_prot, my_mutation)

  barcode_plot(my_prot, start_position = Position, end_position = Position_end, protein_length = length, facet = Protein, cutoffs = c(q = 0.05, L2FC = 0))
  ggsave(paste("./Plots/BarcodeIMEX/", p, ".pdf", sep = ''),  device = "pdf", width = 10, height = 2)

}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Control: Random draw of mutations in the same proteins
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mutation_data_control$CoveredRegions <- sapply(mutation_data_control$uniprot_ID, function(x) return(nchar(FLiP_library[which(FLiP_library$Protein == x), "Sequence"])[1]))

percentage_M_with_FLiP_control <- c()
nb_with_FLiP_control <- c()
nb_without_FLiP_control <- c()
set.seed(123)

# for number_of_control_runs random draws
for(i in c(1:number_of_control_runs)){
  
  # randomly sample mutation sites in the protein
  current_mutation_data_control <- mutation_data_control
  current_mutation_data_control$RandomMutation <- sapply(current_mutation_data_control$CoveredRegions, function(x) return(sample(x, 1)))
  print(current_mutation_data_control$RandomMutation)
  # exclude mutations for which the protein is not part of the FLiP dataset
  if(length(which(is.na(current_mutation_data_control$RandomMutation))) > 0){
    current_mutation_data_control <- current_mutation_data_control[-which(is.na(current_mutation_data_control$RandomMutation)), ]
  }
  
  mutation_detected_control <- rep(FALSE, nrow(current_mutation_data_control)) # vector: is there a peptide overlapping with the mutation
  lowest_FLiP_q_control <- rep(1, nrow(current_mutation_data_control)) # vector: lowest q values of the overlapping peptides
  lowest_FLiP_q_peptide_control <- rep("", nrow(current_mutation_data_control))
  
  for (m in 1:nrow(current_mutation_data_control)){
    
    FLiP_peptides <- FLiP_library[FLiP_library$Protein == current_mutation_data_control$uniprot_ID[m], ]
    mut_covered <- rep(FALSE, nrow(FLiP_peptides))
    
    for(p in 1:nrow(FLiP_peptides)){ # check each peptide of the protein if it overlaps
      start_overlap <- current_mutation_data_control$RandomMutation[m] >= FLiP_peptides$Position[p] - AA_tolerance & current_mutation_data_control$RandomMutation[m] <= FLiP_peptides$Position_end[p] + AA_tolerance # LHS
      mut_covered[p] <- start_overlap > 0 #if the peptide overlaps to the left or right, it is overlapping
      if(sum(mut_covered) > 0){
        mutation_detected_control[m] <- TRUE
      } 
    }
    
    if (mutation_detected_control[m] == TRUE){ # if there are several peptides overlapping, take the one with the lowest q-value of all overlapping peptides to mark the mutation
      lowest_FLiP_q_control[m]<- min(FLiP_peptides$q[mut_covered]) 
      lowest_FLiP_q_peptide_control[m] <- FLiP_peptides$fully_tryptic[which(FLiP_peptides$q == lowest_FLiP_q_control[m])[1]]
    }
  }
  
  current_mutation_data_control$lowest_FLiP_q_control <- lowest_FLiP_q_control
  current_mutation_data_control$lowest_FLiP_q_peptide_control <- lowest_FLiP_q_peptide_control
  current_mutation_data_control <- current_mutation_data_control[mutation_detected_control == TRUE, ] # reduce to mutations detected in dataset
  
  # get rid of mutations that are covered by the same peptide
  overlapping_peptides <- unique(current_mutation_data_control$lowest_FLiP_q_peptide_control)
  mutation_data_reduced <- data.frame()
  for(p in overlapping_peptides){
    mutation_data_reduced <- rbind(mutation_data_reduced, current_mutation_data_control[grep(p, current_mutation_data_control$lowest_FLiP_q_peptide_control)[1], ])
  }
  
  current_mutation_data_control <- mutation_data_reduced
  
  export_total_number_mutations <- nrow(current_mutation_data_control)
  export_actin_number_mutations <- nrow(current_mutation_data_control[which(current_mutation_data_control$uniprot_ID == "P60010"), ])
  
  ### remove actin due to its large number of mutations (extreme outlier)
  if (exclude_actin){
    if(length(which(current_mutation_data_control$uniprot_ID == "P60010")) > 0){
      current_mutation_data_control <- current_mutation_data_control[-which(current_mutation_data_control$uniprot_ID == "P60010"), ]
      #FLiP_library <- FLiP_library[-which(FLiP_library$Protein  == "P60010"), ]
    }
  }
  
  ### calculate how many proteins have a FLiP marker covering them
  sum(current_mutation_data_control$lowest_FLiP_q_control <= 0.05)
  
  percentage_M_with_FLiP_control <- c( percentage_M_with_FLiP_control, sum(current_mutation_data_control$lowest_FLiP_q_control <= 0.05) / length(current_mutation_data_control$lowest_FLiP_q_control))
 
  nb_with_FLiP_control <- c(nb_with_FLiP_control, sum(current_mutation_data_control$lowest_FLiP_q_control <= 0.05))
  nb_without_FLiP_control <- c(nb_without_FLiP_control, sum(current_mutation_data_control$lowest_FLiP_q_control > 0.05))
}

percentage_M_with_FLiP_control
nb_with_FLiP_control
nb_without_FLiP_control

test <- t.test(percentage_M_with_FLiP_control, mu = percentage_M_with_FLiP, alternative = "less") 
print(test)

# some plots
for(p in unique(current_mutation_data_control$uniprot_ID)){
  
  my_prot <- FLiP_library[which(FLiP_library$Protein == p & FLiP_library$q < 0.05), c("Protein", "Position", "Position_end", "length")]
  my_prot$q <- rep(0.001, nrow(my_prot))
  my_prot$L2FC <- rep(0, nrow(my_prot))
  
  my_mutation <- current_mutation_data_control[which(current_mutation_data_control$uniprot_ID == p), "RandomMutation"]
  end <- unlist(sapply(my_mutation, function(x) return(x+1+AA_tolerance)))
  my_mutation <- data.frame("Position" = my_mutation, "Position_end" = end)
  colnames(my_mutation) <- c("Position", "Position_end")
  my_mutation$Protein <- rep(p, nrow(my_mutation))
  my_mutation$length <- rep(my_prot$length[1], nrow(my_mutation))
  my_mutation$q <- rep(1, nrow(my_mutation))
  my_mutation$L2FC <- rep(0, nrow(my_mutation))
  
  my_prot <- rbind(my_prot, my_mutation)
  
  barcode_plot(my_prot, start_position = Position, end_position = Position_end, protein_length = length, facet = Protein, cutoffs = c(q = 0.05, L2FC = 0))
  ggsave(paste("./Plots/BarcodeIMEXRandom//", p, ".pdf", sep = ''),  device = "pdf", width = 10, height = 2)
}


df_percentages <- data.frame("Draw" = c("Random", "IMEX"), 
                          "Percentage" = c(mean(percentage_M_with_FLiP_control)*100, percentage_M_with_FLiP*100), 
                          "SD" = c(sd(percentage_M_with_FLiP_control)*100, 0))
df_percentages <- df_percentages[order(df_percentages$Percentage), ]
df_percentages$Draw <- factor(df_percentages$Draw, levels = df_percentages$Draw)

ggplot(df_percentages) +
  geom_bar( aes(x=Draw, y=Percentage), stat="identity", fill="#0057b7", alpha=0.7) +
  geom_errorbar( aes(x=Draw, ymin=Percentage-SD, ymax=Percentage+SD), width=0.4, colour="#ffc900", alpha=0.9, size=1.3) + ylim(0,100) +
  theme_classic(base_size = 26) + geom_text(aes(x = Draw, y = Percentage, label = round(Percentage, 1)), vjust = -1.5, size = 8) + xlab("")

ggsave(paste("./Plots/IMEX_overlap_", AA_tolerance, "_AA_tolerance.jpg", sep =""), device = "jpg", width = 16 , height = 14, units = "cm")
