# FLiP-MS
This repository contains the code used for the data analysis of the manuscript "Global profiling of protein-protein interaction markers for probing complex protein dynamics". 

### 1_FLiP_PBI_Library
This directory provides all the resources needed to create the FLiP protein binding interface (PBI) library. 

### 2_LiP_WT_control_HU
This directory analyzes the LiP-MS data performed on wild type cells under hydroxyurea(HU)-induced DNA replication stress. It identifies all protein regions that undergo structural changes under this condition and uses the FLiP library to identify proteins that change protein-protein interactions (PPIs).

### 3_Network_Analysis_WT_control_HU
This directory projects the PBI marker changes identified in 2_LiP_WT_control_HU on a protein-protein interaction network, extracts the part of the network closely connected to the hits and clusters the network. This provides information about protein complexes likely to undergo changes in assembly state under HU-stress in wild type cells.

### 4_LiP_Gcn5_Mutant_control_HU
This directory analyzes the LiP-MS data performed on Gcn5 catalytic dead mutant cells under hydroxyurea(HU)-induced DNA replication stress. It identifies all protein regions that undergo structural changes under this condition and uses the FLiP Library to identify proteins that change protein-protein interactions (PPIs).

### 5_Combined_Network_Analysis_WT_Gcn5_Mutant_control_HU
This directory projects the PBI marker changes identified 2_LiP_WT_control_HU and 4_LiP_Gcn5_Mutant_control_HU (only protein regions detected in both datasets considered), projects them on a protein-protein interaction network, extracts the part of the network closely connected to the hits and clusters the network. This provides information about differences in protein complexes likely to underg changes in assembly state under HU-stress in wild type cells versus Gcn5 mutant cells. 

### 6_Acetylation_Analysis
This directory identifies differences in acetylation of peptides in wild type and/or Gcn5 mutant cells under hydroxyurea(HU)-induced DNA replication stress.

### 7_AP_MS
This directory identifies differences in interactors of Ada3 from AP-MS experiments under HU-stress compared to control conditions. 

### 8_Mapping_on_PDB
This directory maps FLiP maker peptides onto mulitmeric PDB structures and calculated the distance to the PBI. A ROC anlaysis is done assessing how good marker peptides located to interfaces. 

### Files
This directory contain databases required in the analyses above. 
