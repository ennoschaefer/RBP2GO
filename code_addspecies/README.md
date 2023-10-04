In this folder, you can find code that was used to add additional species. The code serves as an orientation as for the addition of a new species many steps need to be done manually. However, in [create_table_XX_dataset.Rmd](create_table_XX_dataset.Rmd) one can find a guide for adding a new species.

For orientation purposes, you can find a detailed description for every column of the *table_XX_Dataset.rds* below. It is also listed which databank the information from the column comes from. The *table_XX_Non_Listed_Proteins.rds* columns are identical, except for the column identifying the proteins as RBPs or non-RBPs. Each row represents a protein, all rows constitute the proteome that is considered by RBP2GO.
  - Entry name: Short ID-like label, less consistent then the Uniprot_ID
  - Uniprot_ID: unique ID for each protein given by Uniprot. Most consistent ID, therefore this column is the main identifier to match data from different datasets to the same protein
  - Protein_Name (string): Name of the protein
  - RBP2GO_Score (float): */*
  - Nb_Datasets (int): Amount of studies identifying RBPs
  - Listing_Count (int): */* 
  - AVG10_Int_Listing_Count (float): */* 
  - Mass_kDa (float): Mass in kDa
  - Length_AA (int): amount of AAs
  - pI (float): Isoelectric point
  - Gene_Name (string): Name of the respective gene
  - Alias_Names (string): Other names of the respective gene
  
Each protein can be only marked as one out of the two folllowing
  - Only with poly(A) enrichment ("X"): All studies enriched protein-coding (poly-A) mRNAs to identify this protein as a RBP
  - Only without poly(A) enrichment ("X"): Not a single study enriched protein-coding (poly-A) mRNAs to identify this protein as a RBP
 
  *STUDIES* (at least 1, "X" marks those proteins identified as RBPs by a certain study, column named after publication)
 
  - GO_BP (string): "Gene Ontology - biological process terms", listed, entries separated  by ";" 
  - GO_MF (string): "Gene Ontology - molecular function" terms, listed, entries separated  by ";"
  - GO_CC (string): "Gene Ontology - molecular function" terms, listed, entries separated  by ";"
  - String_ID (string): ID given by String database
  - String_PPI (string): Proteins that are interacted with (listed, entries seperated by ";")
  - String_PPI_Scores (int): Score of the interaction (listed, entries separated  by ";")
  - GO_BP_Tree (string): Parent terms of the GO_BP terms 
  - GO_MF_Tree (string): Parent terms of the GO_MF terms
  - GO_CC_Tree (string): Parent terms of the GO_CC terms 
  
-------------------
Homolog proteins for the protein at hand are only investigated for the species that are covered by RBP2GO.
  - Nb_Homologs (int): How many other species have a homolog to this respective protein
  - Nb_RBP_Homologs (int): How many of the homologs are RBPs
  - List_Homologs (string): Uniprot_IDs of all homolog proteins
  - List_RBP_Homologs (string): Uniprot_IDs of homolog proteins that are RBPs 
-------------------
  
  - InterPro (string): domains contained in the respective protein
  - has_RBD (No/Yes): is any domain RNA binding 
  - RBDs_nb (int): amount of RBDs
  - RBDs_content_fraction (float): */*
  - RBD_list (string): InterPro IDs of the RBDs
  - has_fam_ID (No/Yes): */*
  - Fam_ID_list: InterPro IDs of the family IDs
  - Domain_score (int): */*
  - Composite_Score (float): */*
  - RBP_status (string): RBPs are marked with "RBP"
  - nb_RBDs_individual (string): InterPro IDs of the RBDs
  - IDR_content_fraction (float): */*
  
  
