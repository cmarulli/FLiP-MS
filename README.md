# FLiP-MS
This repository contains most of the code used for the data analysis of the manuscript "Probing protein interactome dynamics using an experimental library of protein complex interfaces". It contains all the input files, intermediate and final files as well as plots and databases needed for the analysis. If you would like to run the code, download the whole repository and store all folders in one folder, relatives paths to other folders in this repository are used in the different analyses. In each directory, run the scripts in increasing order. 


### 1_FLiP_PBI_Library
This directory provides all the resources needed to create the FLiP protein binding interface (PBI) library. 

### 2_LiP_WT_control_HU
This directory analyzes the LiP-MS data performed on wild type cells under hydroxyurea(HU)-induced DNA replication stress. It identifies all protein regions that undergo structural changes under this condition and uses the 1_FLiP_PBI_Library to identify proteins that change protein-protein interactions (PPIs).

### 3_Network_Analysis_WT_control_HU
This directory projects the PBI marker changes identified in 2_LiP_WT_control_HU on a protein-protein interaction network, extracts the part of the network closely connected to the hits and clusters the network. This provides information about protein complexes likely to underg changes in assembly state under HU-stress in wild type cells.

### 4_LiP_Gcn5_Mutant_control_HU
This directory analyzes the LiP-MS data performed on Gcn5 catalytic dead mutant cells under hydroxyurea(HU)-induced DNA replication stress. It identifies all protein regions that undergo structural changes under this condition and uses the 1_FLiP_PBI_Library to identify proteins that change protein-protein interactions (PPIs).

### 5_Combined_Network_Analysis_WT_Gcn5_Mutant_control_HU
This directory projects the PBI marker changes identified 2_LiP_WT_control_HU and 4_LiP_Gcn5_Mutant_control_HU (only protein regions detected in both datasets considered), projects them on a protein-protein interaction network, extracts the part of the network closely connected to the hits and clusters the network. This provides information about differences in protein complexes likely to underg changes in assembly state under HU-stress in wild type cells versus Gcn5 mutant cells. 

### 6_Acetylation_Analysis
This directory identifies differences in acetylation of peptides in wild type and/or Gcn5 mutant cells under hydroxyurea(HU)-induced DNA replication stress.

### Files
This directory contain databases required in the analyses above. 
