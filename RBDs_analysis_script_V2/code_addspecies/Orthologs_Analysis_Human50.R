### RBP2GO - Analysis of the protein orthologs ###

# Load the information that is available for all species
setwd("~/Desktop/RBP2GO/RDA_data_RBP2GO")

load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_HS_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_AT_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_SC_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_EC_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_SaE_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_TB_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_LD_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_LM_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_PF_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_DM_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_CE_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_DaR_alias.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/df_MM_alias.rda")

load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_HS_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_AT_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_SC_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_EC_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_SaE_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_TB_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_LD_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_LM_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_PF_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_DM_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_CE_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_DaR_Dataset.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_MM_Dataset.rda")

load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_HS_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_AT_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_SC_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_EC_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_SaE_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_TB_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_LD_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_LM_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_PF_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_DM_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_CE_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_DaR_Non_Listed_Proteins.rda")
load("~/Desktop/RBP2GO/RDA_data_RBP2GO/table_MM_Non_Listed_Proteins.rda")

# Load the 50% orthologs table
setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs")
Human_UniRef50 <- read.table("UniRef50 human.txt", header=T, sep = "\t", check.names = F, quote = "\"", comment.char = "", stringsAsFactors = F, colClasses = c(rep('character',9)))
# Get the name of the listed organisms UniRef50
n_row_human_UniRef50 <- dim(Human_UniRef50)[1]
organisms_list <- ""
for (i in 1:n_row_human_UniRef50) {
  organisms <- Human_UniRef50[i,"Organisms"]
  organisms <- unlist(strsplit(organisms,"; "))
  organisms_list <- unique(c(organisms_list,organisms))
}
organisms_list_50 <- organisms_list[organisms_list != ""]

# Seek for the 13 RBP2GO species in Human_UniRef50
organisms_list_50[grep("Homo sapiens \\(Human\\)",organisms_list_50,ignore.case = T)] # OK
organisms_list_50[grep("Mus musculus \\(Mouse\\)",organisms_list_50,ignore.case = T)] # OK
organisms_list_50[grep("Danio rerio \\(Zebrafish\\) \\(Brachydanio rerio\\)",organisms_list_50,ignore.case = T)] # OK
organisms_list_50[grep("Drosophila melanogaster \\(Fruit fly\\)",organisms_list_50,ignore.case = T)] # OK
organisms_list_50[grep("Caenorhabditis elegans",organisms_list_50,ignore.case = T)]
organisms_list_50[grep("Saccharomyces cerevisiae",organisms_list_50,ignore.case = T)]
organisms_list_50[grep("Arabidopsis thaliana",organisms_list_50,ignore.case = T)]
organisms_list_50[grep("Trypanosoma brucei",organisms_list_50,ignore.case = T)] # OK
organisms_list_50[grep("Leishmania donovani",organisms_list_50,ignore.case = T)]
organisms_list_50[grep("Leishmania mexicana",organisms_list_50,ignore.case = T)]
organisms_list_50[grep("Plasmodium falciparum",organisms_list_50,ignore.case = T)]
organisms_list_50[grep("Salmonella Typhimurium",organisms_list_50,ignore.case = T)]
organisms_list_50[grep("Escherichia coli",organisms_list_50,ignore.case = T)] # OK

# Rename first column
colnames(Human_UniRef50)[1] <- "Yourlist"

# Make a table for each protein to list the orthologs information
protein_info <- Human_UniRef50[1,]
# Create a dataframe with the information
df_protein <- data.frame(UniProt_ID = protein_info$Yourlist,
                         Cluster_ID = protein_info$`Cluster ID`,
                         Cluster_Name = protein_info$`Cluster name`,
                         Cluster_Size = protein_info$Size,
                         Cluster_Members = unlist(strsplit(protein_info$`Cluster members`,"; ")),
                         #Organisms = unlist(strsplit(protein_info$Organisms,"; ")),
                         Length = protein_info$Length,
                         Identity = protein_info$Identity,
                         stringsAsFactors = F)
df_Human_UniRef50 <- df_protein

for (i in 2:dim(Human_UniRef50)[1]) {
  #for (i in 20081:20083) {  
  protein_info <- Human_UniRef50[i,]
  nb_protein <- length(unlist(strsplit(protein_info$Yourlist,",")))
  if (nb_protein == 1) {
    df_protein <- data.frame(UniProt_ID = protein_info$Yourlist,
                             Cluster_ID = protein_info$`Cluster ID`,
                             Cluster_Name = protein_info$`Cluster name`,
                             Cluster_Size = protein_info$Size,
                             Cluster_Members = unlist(strsplit(protein_info$`Cluster members`,"; ")),
                             #Organisms = unlist(strsplit(protein_info$Organisms,"; ")),
                             Length = protein_info$Length,
                             Identity = protein_info$Identity,
                             stringsAsFactors = F)
    df_Human_UniRef50 <- rbind(df_Human_UniRef50,df_protein)
  } else {
    for (j in 1:nb_protein) {
      df_protein <- data.frame(UniProt_ID = unlist(strsplit(protein_info$Yourlist,","))[j],
                               Cluster_ID = protein_info$`Cluster ID`,
                               Cluster_Name = protein_info$`Cluster name`,
                               Cluster_Size = protein_info$Size,
                               Cluster_Members = unlist(strsplit(protein_info$`Cluster members`,"; ")),
                               #Organisms = unlist(strsplit(protein_info$Organisms,"; ")),
                               Length = protein_info$Length,
                               Identity = protein_info$Identity,
                               stringsAsFactors = F)
      df_Human_UniRef50 <- rbind(df_Human_UniRef50,df_protein)
    }
  }
}

setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
save(df_Human_UniRef50,file = "df_Human_UniRef50.rda")










# Select only the interesting information (only the 13 species of the RBP2GO database)
# Return to the working directory
setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
load("df_Human_UniRef50.rda")

# total number of proteins
dim(table_HS_Dataset)[1]+dim(table_HS_Non_Listed_Proteins)[1]
# total number of RBPs
dim(table_HS_Dataset)[1]
# total number of non-RBPs
dim(table_HS_Non_Listed_Proteins)[1]

# Remove proteins ortholog with themselves
df_Human_UniRef50 <- df_Human_UniRef50[!df_Human_UniRef50$UniProt_ID == df_Human_UniRef50$Cluster_Members,]

# total number of proteins with orthologs information
length(unique(df_Human_UniRef50$UniProt_ID))
# total number of RBPs with orthologs information
length(unique(df_Human_UniRef50[df_Human_UniRef50$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]$UniProt_ID))
# total number of non-RBPs with orthologs information
length(unique(df_Human_UniRef50[df_Human_UniRef50$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]$UniProt_ID))

# Check for proteins from same species and other species - Homo sapiens
df_Human_UniRef50_HS <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_HS_alias$alias,]
df_Human_UniRef50_HS$Organisms <- "Homo sapiens"
# Get the information whether the ortholog is RBP-listed or not
df_Human_UniRef50_HS_RBP <- df_Human_UniRef50_HS[df_Human_UniRef50_HS$Cluster_Members %in% table_HS_Dataset$Uniprot_ID,]
df_Human_UniRef50_HS_RBP$Times_Listed <- apply(df_Human_UniRef50_HS_RBP,1,function(x) {
  table_HS_Dataset[table_HS_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_HS_Non_RBP <- df_Human_UniRef50_HS[!df_Human_UniRef50_HS$Cluster_Members %in% table_HS_Dataset$Uniprot_ID,]
df_Human_UniRef50_HS_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Mus musculus
df_Human_UniRef50_MM <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_MM_alias$alias,]
df_Human_UniRef50_MM$Organisms <- "Mus musculus"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_MM_RBP <- df_Human_UniRef50_MM[df_Human_UniRef50_MM$Cluster_Members %in% table_MM_Dataset$Uniprot_ID,]
df_Human_UniRef50_MM_RBP$Times_Listed <- apply(df_Human_UniRef50_MM_RBP,1,function(x) {
  table_MM_Dataset[table_MM_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_MM_Non_RBP <- df_Human_UniRef50_MM[!df_Human_UniRef50_MM$Cluster_Members %in% table_MM_Dataset$Uniprot_ID,]
df_Human_UniRef50_MM_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Danio rerio
df_Human_UniRef50_DaR <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_DaR_alias$alias,]
df_Human_UniRef50_DaR$Organisms <- "Danio rerio"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_DaR_RBP <- df_Human_UniRef50_DaR[df_Human_UniRef50_DaR$Cluster_Members %in% table_DaR_Dataset$Uniprot_ID,]
df_Human_UniRef50_DaR_RBP$Times_Listed <- apply(df_Human_UniRef50_DaR_RBP,1,function(x) {
  table_DaR_Dataset[table_DaR_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_DaR_Non_RBP <- df_Human_UniRef50_DaR[!df_Human_UniRef50_DaR$Cluster_Members %in% table_DaR_Dataset$Uniprot_ID,]
df_Human_UniRef50_DaR_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Drosophila melanogaster
df_Human_UniRef50_DM <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_DM_alias$alias,]
df_Human_UniRef50_DM$Organisms <- "Drosophila melanogaster"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_DM_RBP <- df_Human_UniRef50_DM[df_Human_UniRef50_DM$Cluster_Members %in% table_DM_Dataset$Uniprot_ID,]
df_Human_UniRef50_DM_RBP$Times_Listed <- apply(df_Human_UniRef50_DM_RBP,1,function(x) {
  table_DM_Dataset[table_DM_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_DM_Non_RBP <- df_Human_UniRef50_DM[!df_Human_UniRef50_DM$Cluster_Members %in% table_DM_Dataset$Uniprot_ID,]
df_Human_UniRef50_DM_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Caenorhabditis elegans
df_Human_UniRef50_CE <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_CE_alias$alias,]
df_Human_UniRef50_CE$Organisms <- "Caenorhabditis elegans"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_CE_RBP <- df_Human_UniRef50_CE[df_Human_UniRef50_CE$Cluster_Members %in% table_CE_Dataset$Uniprot_ID,]
df_Human_UniRef50_CE_RBP$Times_Listed <- apply(df_Human_UniRef50_CE_RBP,1,function(x) {
  table_CE_Dataset[table_CE_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_CE_Non_RBP <- df_Human_UniRef50_CE[!df_Human_UniRef50_CE$Cluster_Members %in% table_CE_Dataset$Uniprot_ID,]
df_Human_UniRef50_CE_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Saccharomyces cerevisiae
df_Human_UniRef50_SC <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_SC_alias$alias,]
df_Human_UniRef50_SC$Organisms <- "Saccharomyces cerevisiae"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_SC_RBP <- df_Human_UniRef50_SC[df_Human_UniRef50_SC$Cluster_Members %in% table_SC_Dataset$Uniprot_ID,]
df_Human_UniRef50_SC_RBP$Times_Listed <- apply(df_Human_UniRef50_SC_RBP,1,function(x) {
  table_SC_Dataset[table_SC_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_SC_Non_RBP <- df_Human_UniRef50_SC[!df_Human_UniRef50_SC$Cluster_Members %in% table_SC_Dataset$Uniprot_ID,]
df_Human_UniRef50_SC_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Arabidopsis thaliana
df_Human_UniRef50_AT <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_AT_alias$alias,]
df_Human_UniRef50_AT$Organisms <- "Arabidopsis thaliana"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_AT_RBP <- df_Human_UniRef50_AT[df_Human_UniRef50_AT$Cluster_Members %in% table_AT_Dataset$Uniprot_ID,]
df_Human_UniRef50_AT_RBP$Times_Listed <- apply(df_Human_UniRef50_AT_RBP,1,function(x) {
  table_AT_Dataset[table_AT_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_AT_Non_RBP <- df_Human_UniRef50_AT[!df_Human_UniRef50_AT$Cluster_Members %in% table_AT_Dataset$Uniprot_ID,]
df_Human_UniRef50_AT_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Trypanosoma brucei
df_Human_UniRef50_TB <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_TB_alias$alias,]
df_Human_UniRef50_TB$Organisms <- "Trypanosoma brucei"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_TB_RBP <- df_Human_UniRef50_TB[df_Human_UniRef50_TB$Cluster_Members %in% table_TB_Dataset$Uniprot_ID,]
df_Human_UniRef50_TB_RBP$Times_Listed <- apply(df_Human_UniRef50_TB_RBP,1,function(x) {
  table_TB_Dataset[table_TB_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_TB_Non_RBP <- df_Human_UniRef50_TB[!df_Human_UniRef50_TB$Cluster_Members %in% table_TB_Dataset$Uniprot_ID,]
df_Human_UniRef50_TB_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Leishmania donovani
df_Human_UniRef50_LD <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_LD_alias$alias,]
df_Human_UniRef50_LD$Organisms <- "Leishmania donovani"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_LD_RBP <- df_Human_UniRef50_LD[df_Human_UniRef50_LD$Cluster_Members %in% table_LD_Dataset$Uniprot_ID,]
df_Human_UniRef50_LD_RBP$Times_Listed <- apply(df_Human_UniRef50_LD_RBP,1,function(x) {
  table_LD_Dataset[table_LD_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_LD_Non_RBP <- df_Human_UniRef50_LD[!df_Human_UniRef50_LD$Cluster_Members %in% table_LD_Dataset$Uniprot_ID,]
df_Human_UniRef50_LD_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Leishmania mexicana
df_Human_UniRef50_LM <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_LM_alias$alias,]
df_Human_UniRef50_LM$Organisms <- "Leishmania mexicana"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_LM_RBP <- df_Human_UniRef50_LM[df_Human_UniRef50_LM$Cluster_Members %in% table_LM_Dataset$Uniprot_ID,]
df_Human_UniRef50_LM_RBP$Times_Listed <- apply(df_Human_UniRef50_LM_RBP,1,function(x) {
  table_LM_Dataset[table_LM_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_LM_Non_RBP <- df_Human_UniRef50_LM[!df_Human_UniRef50_LM$Cluster_Members %in% table_LM_Dataset$Uniprot_ID,]
df_Human_UniRef50_LM_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Plasmodium falciparum
df_Human_UniRef50_PF <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_PF_alias$alias,]
df_Human_UniRef50_PF$Organisms <- "Plasmodium falciparum"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_PF_RBP <- df_Human_UniRef50_PF[df_Human_UniRef50_PF$Cluster_Members %in% table_PF_Dataset$Uniprot_ID,]
df_Human_UniRef50_PF_RBP$Times_Listed <- apply(df_Human_UniRef50_PF_RBP,1,function(x) {
  table_PF_Dataset[table_PF_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_PF_Non_RBP <- df_Human_UniRef50_PF[!df_Human_UniRef50_PF$Cluster_Members %in% table_PF_Dataset$Uniprot_ID,]
df_Human_UniRef50_PF_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Salmonella enterica subsp. enterica ser. Typhimurium
df_Human_UniRef50_SaE <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_SaE_alias$alias,]
df_Human_UniRef50_SaE$Organisms <- "Salmonella Typhimurium"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_SaE_RBP <- df_Human_UniRef50_SaE[df_Human_UniRef50_SaE$Cluster_Members %in% table_SaE_Dataset$Uniprot_ID,]
df_Human_UniRef50_SaE_RBP$Times_Listed <- apply(df_Human_UniRef50_SaE_RBP,1,function(x) {
  table_SaE_Dataset[table_SaE_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_SaE_Non_RBP <- df_Human_UniRef50_SaE[!df_Human_UniRef50_SaE$Cluster_Members %in% table_SaE_Dataset$Uniprot_ID,]
df_Human_UniRef50_SaE_Non_RBP$Times_Listed <- 0

# Check for proteins from other species - Escherichia coli
df_Human_UniRef50_EC <- df_Human_UniRef50[df_Human_UniRef50$Cluster_Members %in% df_EC_alias$alias,]
df_Human_UniRef50_EC$Organisms <- "Escherichia coli"
# Get the information whether the protein is RBP-listed or not
df_Human_UniRef50_EC_RBP <- df_Human_UniRef50_EC[df_Human_UniRef50_EC$Cluster_Members %in% table_EC_Dataset$Uniprot_ID,]
df_Human_UniRef50_EC_RBP$Times_Listed <- apply(df_Human_UniRef50_EC_RBP,1,function(x) {
  table_EC_Dataset[table_EC_Dataset$Uniprot_ID == x["Cluster_Members"],]$Times_Listed
})
df_Human_UniRef50_EC_Non_RBP <- df_Human_UniRef50_EC[!df_Human_UniRef50_EC$Cluster_Members %in% table_EC_Dataset$Uniprot_ID,]
df_Human_UniRef50_EC_Non_RBP$Times_Listed <- 0

# Table for all the proteins in the 13 species
df_Human_UniRef50_All <- rbind(df_Human_UniRef50_HS_RBP,
                               df_Human_UniRef50_HS_Non_RBP,
                               df_Human_UniRef50_MM_RBP,
                               df_Human_UniRef50_MM_Non_RBP,
                               df_Human_UniRef50_AT_RBP,
                               df_Human_UniRef50_AT_Non_RBP,
                               df_Human_UniRef50_SC_RBP,
                               df_Human_UniRef50_SC_Non_RBP,
                               df_Human_UniRef50_EC_RBP,
                               df_Human_UniRef50_EC_Non_RBP,
                               df_Human_UniRef50_SaE_RBP,
                               df_Human_UniRef50_SaE_Non_RBP,
                               df_Human_UniRef50_TB_RBP,
                               df_Human_UniRef50_TB_Non_RBP,
                               df_Human_UniRef50_LD_RBP,
                               df_Human_UniRef50_LD_Non_RBP,
                               df_Human_UniRef50_LM_RBP,
                               df_Human_UniRef50_LM_Non_RBP,
                               df_Human_UniRef50_PF_RBP,
                               df_Human_UniRef50_PF_Non_RBP,
                               df_Human_UniRef50_DM_RBP,
                               df_Human_UniRef50_DM_Non_RBP,
                               df_Human_UniRef50_CE_RBP,
                               df_Human_UniRef50_CE_Non_RBP,
                               df_Human_UniRef50_DaR_RBP,
                               df_Human_UniRef50_DaR_Non_RBP)
save(df_Human_UniRef50_All,file = "df_Human_UniRef50_All.rda")

# Table for all the proteins with RBP orthologs in the 13 species
df_Human_UniRef50_All_RBP <- rbind(df_Human_UniRef50_HS_RBP,
                                   df_Human_UniRef50_MM_RBP,
                                   df_Human_UniRef50_AT_RBP,
                                   df_Human_UniRef50_SC_RBP,
                                   df_Human_UniRef50_EC_RBP,
                                   df_Human_UniRef50_SaE_RBP,
                                   df_Human_UniRef50_TB_RBP,
                                   df_Human_UniRef50_LD_RBP,
                                   df_Human_UniRef50_LM_RBP,
                                   df_Human_UniRef50_PF_RBP,
                                   df_Human_UniRef50_DM_RBP,
                                   df_Human_UniRef50_CE_RBP,
                                   df_Human_UniRef50_DaR_RBP)
save(df_Human_UniRef50_All_RBP,file = "df_Human_UniRef50_All_RBP.rda")

# Some statistics
# total number of proteins with RBP orthologs
length(unique(df_Human_UniRef50_All[df_Human_UniRef50_All$Times_Listed >0,]$UniProt_ID))
# Total number of RBPs with RBP orthologs
length(unique(df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]$UniProt_ID))
# Total number of non-RBPs with RBP orthologs
length(unique(df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]$UniProt_ID))


#Number of proteins with orthologs in human
length(unique(df_Human_UniRef50_HS$UniProt_ID))

# Number of proteins with RBP orthologs in human
length(unique(df_Human_UniRef50_HS_RBP$UniProt_ID))

#Number of proteins with orthologs in mouse
length(unique(df_Human_UniRef50_MM$UniProt_ID))

# Number of proteins with RBP orthologs in mouse
length(unique(df_Human_UniRef50_MM_RBP$UniProt_ID))

#Number of proteins with orthologs in Zebrafish
length(unique(df_Human_UniRef50_DaR$UniProt_ID))

# Number of proteins with RBP orthologs in Zebrafish
length(unique(df_Human_UniRef50_DaR_RBP$UniProt_ID))

# other species
length(unique(df_Human_UniRef50_AT$UniProt_ID))
length(unique(df_Human_UniRef50_AT_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_SC$UniProt_ID))
length(unique(df_Human_UniRef50_SC_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_EC$UniProt_ID))
length(unique(df_Human_UniRef50_EC_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_SaE$UniProt_ID))
length(unique(df_Human_UniRef50_SaE_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_TB$UniProt_ID))
length(unique(df_Human_UniRef50_TB_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_LD$UniProt_ID))
length(unique(df_Human_UniRef50_LD_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_LM$UniProt_ID))
length(unique(df_Human_UniRef50_LM_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_PF$UniProt_ID))
length(unique(df_Human_UniRef50_PF_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_DM$UniProt_ID))
length(unique(df_Human_UniRef50_DM_RBP$UniProt_ID))
length(unique(df_Human_UniRef50_CE$UniProt_ID))
length(unique(df_Human_UniRef50_CE_RBP$UniProt_ID))

####################################################################
# Analysis of RBPs and Non-RBPs without the intraspecies information
# Restriction to the 13 species found in RBP2GO - file df_Human_UniRef50_All.rda

setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
load("df_Human_UniRef50_All.rda")
load("df_Human_UniRef50_All_RBP.rda")
HS_proteins <- c(table_HS_Dataset$Uniprot_ID,table_HS_Non_Listed_Proteins$Uniprot_ID)

# Total proteins with orthologs (no intra-species)
length(unique(df_Human_UniRef50_All[!df_Human_UniRef50_All$Cluster_Members %in% HS_proteins,]$UniProt_ID))
# Total proteins RBPs with orthologs
HS_RBP <- df_Human_UniRef50_All[df_Human_UniRef50_All$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]
length(unique(HS_RBP[!HS_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))
RBP_All <- length(unique(HS_RBP[!HS_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))
# Total proteins non-RBPs with orthologs
HS_Non_RBP <- df_Human_UniRef50_All[df_Human_UniRef50_All$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]
length(unique(HS_Non_RBP[!HS_Non_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))
Non_RBP_All <- length(unique(HS_Non_RBP[!HS_Non_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))

# Proteins with RBP orthologs
length(unique(df_Human_UniRef50_All_RBP[!df_Human_UniRef50_All_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))
# Non-RBPs with RBP orthologs
HS_Non_RBP_RBP <- df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]
length(unique(HS_Non_RBP_RBP[!HS_Non_RBP_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))
Non_RBP_RBP <- length(unique(HS_Non_RBP_RBP[!HS_Non_RBP_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))
# RBPs with RBP orthologs
HS_RBP_RBP <- df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]
length(unique(HS_RBP_RBP[!HS_RBP_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))
RBP_RBP <- length(unique(HS_RBP_RBP[!HS_RBP_RBP$Cluster_Members %in% HS_proteins,]$UniProt_ID))

# Enrichment / Fisher's extact test 
RBP_enrichment <- matrix(c(Non_RBP_All,RBP_All,Non_RBP_RBP,RBP_RBP),
                         nrow = 2,
                         dimnames = list(Proteome = c("Non-RBPs","RBPs"),
                                         Size = c("All","RBP-orthologs")
                         )
)
RBP_enrichment
# attributes(fisher.test(GO_enrichment))
fisher.test(RBP_enrichment)$estimate
fisher.test(RBP_enrichment)$p.value



####################################################################
# Analysis of RBPs and Non-RBPs WITH the intraspecies information
# Restriction to the 13 species found in RBP2GO - file df_Human_UniRef50_All.rda

load("df_Human_UniRef50_All.rda")
load("df_Human_UniRef50_All_RBP.rda")

# Total proteins with orthologs (no intra-species)
length(unique(df_Human_UniRef50_All$UniProt_ID))
# Total proteins non-RBPs with orthologs
HS_Non_RBP <- df_Human_UniRef50_All[df_Human_UniRef50_All$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]
length(unique(HS_Non_RBP$UniProt_ID))
Non_RBP_All <- length(unique(HS_Non_RBP$UniProt_ID))
# Total proteins RBPs with orthologs
HS_RBP <- df_Human_UniRef50_All[df_Human_UniRef50_All$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]
length(unique(HS_RBP$UniProt_ID))
RBP_All <- length(unique(HS_RBP$UniProt_ID))
# Proteins with RBP orthologs
length(unique(df_Human_UniRef50_All_RBP$UniProt_ID))
# Non-RBPs with RBP orthologs
HS_Non_RBP_RBP <- df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]
length(unique(HS_Non_RBP_RBP$UniProt_ID))
Non_RBP_RBP <- length(unique(HS_Non_RBP_RBP$UniProt_ID))
# RBPs with RBP orthologs
HS_RBP_RBP <- df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]
length(unique(HS_RBP_RBP$UniProt_ID))
RBP_RBP <- length(unique(HS_RBP_RBP$UniProt_ID))

# Enrichment / Fisher's extact test 
RBP_enrichment <- matrix(c(Non_RBP_All,RBP_All,Non_RBP_RBP,RBP_RBP),
                         nrow = 2,
                         dimnames = list(Proteome = c("Non-RBPs","RBPs"),
                                         Size = c("All","RBP-orthologs")
                         )
)
RBP_enrichment
# attributes(fisher.test(GO_enrichment))
fisher.test(RBP_enrichment)$estimate
fisher.test(RBP_enrichment)$p.value






####################################################################
####################################################################
# Analysis of RBPs and Non-RBPs othologs (other species) and paralogs (same species)
# Restriction to the 13 species found in RBP2GO - file df_Human_UniRef50_All.rda

setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
load("df_Human_UniRef50_All.rda")
load("df_Human_UniRef50_All_RBP.rda")
HS_proteins <- c(table_HS_Dataset$Uniprot_ID,table_HS_Non_Listed_Proteins$Uniprot_ID)

# Total proteins with orthologs (no intra-species)
length(unique(df_Human_UniRef50_All[!df_Human_UniRef50_All$Cluster_Members %in% HS_proteins,]$UniProt_ID))

# Total number of RBP proteins with homologs
HS_RBP <- df_Human_UniRef50_All[df_Human_UniRef50_All$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]
length(unique(HS_RBP$UniProt_ID))
# Total number of RBP proteins with orthologs
HS_RBP_ortho <- HS_RBP[!HS_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_RBP_ortho$UniProt_ID))
# Total number of RBP proteins with paralogs
HS_RBP_para <- HS_RBP[HS_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_RBP_para$UniProt_ID))
# Check
length(unique(c(HS_RBP_ortho$UniProt_ID,HS_RBP_para$UniProt_ID)))

# Total number of Non-RBP proteins with homologs
HS_Non_RBP <- df_Human_UniRef50_All[df_Human_UniRef50_All$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]
length(unique(HS_Non_RBP$UniProt_ID))
# Total number of Non-RBP proteins with orthologs
HS_Non_RBP_ortho <- HS_Non_RBP[!HS_Non_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_Non_RBP_ortho$UniProt_ID))
# Total number of Non-RBP proteins with paralogs
HS_Non_RBP_para <- HS_Non_RBP[HS_Non_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_Non_RBP_para$UniProt_ID))
# Check
length(unique(c(HS_Non_RBP_ortho$UniProt_ID,HS_Non_RBP_para$UniProt_ID)))

# Total number of RBP proteins with RBP homologs
HS_RBP_RBP <- df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Dataset$Uniprot_ID,]
length(unique(HS_RBP_RBP$UniProt_ID))
# Total number of RBP proteins with RBP orthologs
HS_RBP_ortho_RBP <- HS_RBP_RBP[!HS_RBP_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_RBP_ortho_RBP$UniProt_ID))
# Total number of RBP proteins with RBP paralogs
HS_RBP_para_RBP <- HS_RBP_RBP[HS_RBP_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_RBP_para_RBP$UniProt_ID))
# Check
length(unique(c(HS_RBP_ortho_RBP$UniProt_ID,HS_RBP_para_RBP$UniProt_ID)))

# Total number of Non-RBP proteins with RBP homologs
HS_Non_RBP_RBP <- df_Human_UniRef50_All_RBP[df_Human_UniRef50_All_RBP$UniProt_ID %in% table_HS_Non_Listed_Proteins$Uniprot_ID,]
length(unique(HS_Non_RBP_RBP$UniProt_ID))
# Total number of Non-RBP proteins with RBP orthologs
HS_Non_RBP_ortho_RBP <- HS_Non_RBP_RBP[!HS_Non_RBP_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_Non_RBP_ortho_RBP$UniProt_ID))
# Total number of Non-RBP proteins with RBP paralogs
HS_Non_RBP_para_RBP <- HS_Non_RBP_RBP[HS_Non_RBP_RBP$Cluster_Members %in% HS_proteins,]
length(unique(HS_Non_RBP_para_RBP$UniProt_ID))
# Check
length(unique(c(HS_Non_RBP_ortho_RBP$UniProt_ID,HS_Non_RBP_para_RBP$UniProt_ID)))

# Enrichment / Fisher's extact test orthologs
RBP_enrichment <- matrix(c(length(unique(HS_Non_RBP_ortho$UniProt_ID)),
                           length(unique(HS_RBP_ortho$UniProt_ID)),
                           length(unique(HS_Non_RBP_ortho_RBP$UniProt_ID)),
                           length(unique(HS_RBP_ortho_RBP$UniProt_ID))),
                         nrow = 2,
                         dimnames = list(Proteome = c("Non-RBPs","RBPs"),
                                         Size = c("All","RBP-orthologs")
                         )
)
RBP_enrichment
# attributes(fisher.test(GO_enrichment))
fisher.test(RBP_enrichment)$estimate
fisher.test(RBP_enrichment)$p.value

# Enrichment / Fisher's extact test paralogs
RBP_enrichment <- matrix(c(length(unique(HS_Non_RBP_para$UniProt_ID)),
                           length(unique(HS_RBP_para$UniProt_ID)),
                           length(unique(HS_Non_RBP_para_RBP$UniProt_ID)),
                           length(unique(HS_RBP_para_RBP$UniProt_ID))),
                         nrow = 2,
                         dimnames = list(Proteome = c("Non-RBPs","RBPs"),
                                         Size = c("All","RBP-paralogs")
                         )
)
RBP_enrichment
# attributes(fisher.test(GO_enrichment))
fisher.test(RBP_enrichment)$estimate
fisher.test(RBP_enrichment)$p.value


####################################################################
####################################################################
# Preparation of the table for the RBP2GO database

# Load the saved information
setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
load("df_Human_UniRef50_All.rda")
load("df_Human_UniRef50_All_RBP.rda")
df_HS_UniRef50_All <- df_Human_UniRef50_All
df_HS_UniRef50_All_RBP <- df_Human_UniRef50_All_RBP
save(df_HS_UniRef50_All,file = "df_HS_UniRef50_All.rda")
save(df_HS_UniRef50_All_RBP,file = "df_HS_UniRef50_All_RBP.rda")

setwd("~/Desktop/RBP2GO/RDA_data_RBP2GO_New")
load("table_HS_Dataset.rda")
load("table_HS_Non_Listed_Proteins.rda")
load("table_MM_Dataset.rda")
load("table_MM_Non_Listed_Proteins.rda")
load("table_DaR_Dataset.rda")
load("table_DaR_Non_Listed_Proteins.rda")
load("table_DM_Dataset.rda")
load("table_DM_Non_Listed_Proteins.rda")
load("table_CE_Dataset.rda")
load("table_CE_Non_Listed_Proteins.rda")
load("table_SC_Dataset.rda")
load("table_SC_Non_Listed_Proteins.rda")
load("table_AT_Dataset.rda")
load("table_AT_Non_Listed_Proteins.rda")
load("table_TB_Dataset.rda")
load("table_TB_Non_Listed_Proteins.rda")
load("table_LD_Dataset.rda")
load("table_LD_Non_Listed_Proteins.rda")
load("table_LM_Dataset.rda")
load("table_LM_Non_Listed_Proteins.rda")
load("table_PF_Dataset.rda")
load("table_PF_Non_Listed_Proteins.rda")
load("table_SaE_Dataset.rda")
load("table_SaE_Non_Listed_Proteins.rda")
load("table_EC_Dataset.rda")
load("table_EC_Non_Listed_Proteins.rda")


# Make a huge table with all the proteins from all the organisms
table_RBP2GO <- rbind(data.frame(Entry_Name = table_HS_Dataset$Entry_Name, 
                                 Uniprot_ID = table_HS_Dataset$Uniprot_ID,
                                 Protein_Name = table_HS_Dataset$Protein_Name,
                                 Mass_kDa = table_HS_Dataset$Mass_kDa,
                                 Length_AA = table_HS_Dataset$Length_AA,
                                 RBP2GO_Score = table_HS_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_HS_Dataset$Nb_Datasets,
                                 Listing_Count = table_HS_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_HS_Dataset$AVG10_Int_Listing_Count
                                 ),
                      data.frame(Entry_Name = table_HS_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_HS_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_HS_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_HS_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_HS_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_HS_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_HS_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_HS_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_HS_Non_Listed_Proteins$AVG10_Int_Listing_Count
                                 ),
                      data.frame(Entry_Name = table_MM_Dataset$Entry_Name, 
                                 Uniprot_ID = table_MM_Dataset$Uniprot_ID,
                                 Protein_Name = table_MM_Dataset$Protein_Name,
                                 Mass_kDa = table_MM_Dataset$Mass_kDa,
                                 Length_AA = table_MM_Dataset$Length_AA,
                                 RBP2GO_Score = table_MM_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_MM_Dataset$Nb_Datasets,
                                 Listing_Count = table_MM_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_MM_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_MM_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_MM_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_MM_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_MM_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_MM_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_MM_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_MM_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_MM_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_MM_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_DaR_Dataset$Entry_Name, 
                                 Uniprot_ID = table_DaR_Dataset$Uniprot_ID,
                                 Protein_Name = table_DaR_Dataset$Protein_Name,
                                 Mass_kDa = table_DaR_Dataset$Mass_kDa,
                                 Length_AA = table_DaR_Dataset$Length_AA,
                                 RBP2GO_Score = table_DaR_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_DaR_Dataset$Nb_Datasets,
                                 Listing_Count = table_DaR_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_DaR_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_DaR_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_DaR_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_DaR_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_DaR_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_DaR_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_DaR_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_DaR_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_DaR_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_DaR_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_DM_Dataset$Entry_Name, 
                                 Uniprot_ID = table_DM_Dataset$Uniprot_ID,
                                 Protein_Name = table_DM_Dataset$Protein_Name,
                                 Mass_kDa = table_DM_Dataset$Mass_kDa,
                                 Length_AA = table_DM_Dataset$Length_AA,
                                 Nb_Datasets = table_DM_Dataset$Nb_Datasets,
                                 RBP2GO_Score = table_DM_Dataset$RBP2GO_Score, 
                                 Listing_Count = table_DM_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_DM_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_DM_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_DM_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_DM_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_DM_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_DM_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_DM_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_DM_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_DM_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_DM_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_CE_Dataset$Entry_Name, 
                                 Uniprot_ID = table_CE_Dataset$Uniprot_ID,
                                 Protein_Name = table_CE_Dataset$Protein_Name,
                                 Mass_kDa = table_CE_Dataset$Mass_kDa,
                                 Length_AA = table_CE_Dataset$Length_AA,
                                 RBP2GO_Score = table_CE_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_CE_Dataset$Nb_Datasets,
                                 Listing_Count = table_CE_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_CE_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_CE_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_CE_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_CE_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_CE_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_CE_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_CE_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_CE_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_CE_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_CE_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_SC_Dataset$Entry_Name, 
                                 Uniprot_ID = table_SC_Dataset$Uniprot_ID,
                                 Protein_Name = table_SC_Dataset$Protein_Name,
                                 Mass_kDa = table_SC_Dataset$Mass_kDa,
                                 Length_AA = table_SC_Dataset$Length_AA,
                                 RBP2GO_Score = table_SC_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_SC_Dataset$Nb_Datasets,
                                 Listing_Count = table_SC_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_SC_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_SC_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_SC_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_SC_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_SC_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_SC_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_SC_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_SC_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_SC_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_SC_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_AT_Dataset$Entry_Name, 
                                 Uniprot_ID = table_AT_Dataset$Uniprot_ID,
                                 Protein_Name = table_AT_Dataset$Protein_Name,
                                 Mass_kDa = table_AT_Dataset$Mass_kDa,
                                 Length_AA = table_AT_Dataset$Length_AA,
                                 RBP2GO_Score = table_AT_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_AT_Dataset$Nb_Datasets,
                                 Listing_Count = table_AT_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_AT_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_AT_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_AT_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_AT_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_AT_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_AT_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_AT_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_AT_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_AT_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_AT_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_TB_Dataset$Entry_Name, 
                                 Uniprot_ID = table_TB_Dataset$Uniprot_ID,
                                 Protein_Name = table_TB_Dataset$Protein_Name,
                                 Mass_kDa = table_TB_Dataset$Mass_kDa,
                                 Length_AA = table_TB_Dataset$Length_AA,
                                 RBP2GO_Score = table_TB_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_TB_Dataset$Nb_Datasets,
                                 Listing_Count = table_TB_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_TB_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_TB_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_TB_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_TB_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_TB_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_TB_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_TB_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_TB_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_TB_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_TB_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_LD_Dataset$Entry_Name, 
                                 Uniprot_ID = table_LD_Dataset$Uniprot_ID,
                                 Protein_Name = table_LD_Dataset$Protein_Name,
                                 Mass_kDa = table_LD_Dataset$Mass_kDa,
                                 Length_AA = table_LD_Dataset$Length_AA,
                                 RBP2GO_Score = table_LD_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_LD_Dataset$Nb_Datasets,
                                 Listing_Count = table_LD_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_LD_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_LD_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_LD_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_LD_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_LD_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_LD_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_LD_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_LD_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_LD_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_LD_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_LM_Dataset$Entry_Name, 
                                 Uniprot_ID = table_LM_Dataset$Uniprot_ID,
                                 Protein_Name = table_LM_Dataset$Protein_Name,
                                 Mass_kDa = table_LM_Dataset$Mass_kDa,
                                 Length_AA = table_LM_Dataset$Length_AA,
                                 RBP2GO_Score = table_LM_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_LM_Dataset$Nb_Datasets,
                                 Listing_Count = table_LM_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_LM_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_LM_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_LM_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_LM_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_LM_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_LM_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_LM_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_LM_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_LM_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_LM_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_PF_Dataset$Entry_Name, 
                                 Uniprot_ID = table_PF_Dataset$Uniprot_ID,
                                 Protein_Name = table_PF_Dataset$Protein_Name,
                                 Mass_kDa = table_PF_Dataset$Mass_kDa,
                                 Length_AA = table_PF_Dataset$Length_AA,
                                 RBP2GO_Score = table_PF_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_PF_Dataset$Nb_Datasets,
                                 Listing_Count = table_PF_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_PF_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_PF_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_PF_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_PF_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_PF_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_PF_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_PF_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_PF_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_PF_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_PF_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_SaE_Dataset$Entry_Name, 
                                 Uniprot_ID = table_SaE_Dataset$Uniprot_ID,
                                 Protein_Name = table_SaE_Dataset$Protein_Name,
                                 Mass_kDa = table_SaE_Dataset$Mass_kDa,
                                 Length_AA = table_SaE_Dataset$Length_AA,
                                 RBP2GO_Score = table_SaE_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_SaE_Dataset$Nb_Datasets,
                                 Listing_Count = table_SaE_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_SaE_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_SaE_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_SaE_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_SaE_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_SaE_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_SaE_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_SaE_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_SaE_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_SaE_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_SaE_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_EC_Dataset$Entry_Name, 
                                 Uniprot_ID = table_EC_Dataset$Uniprot_ID,
                                 Protein_Name = table_EC_Dataset$Protein_Name,
                                 Mass_kDa = table_EC_Dataset$Mass_kDa,
                                 Length_AA = table_EC_Dataset$Length_AA,
                                 RBP2GO_Score = table_EC_Dataset$RBP2GO_Score,
                                 Nb_Datasets = table_EC_Dataset$Nb_Datasets,
                                 Listing_Count = table_EC_Dataset$Listing_Count,
                                 AVG10_Int_Listing_Count = table_EC_Dataset$AVG10_Int_Listing_Count
                      ),
                      data.frame(Entry_Name = table_EC_Non_Listed_Proteins$Entry_Name, 
                                 Uniprot_ID = table_EC_Non_Listed_Proteins$Uniprot_ID,
                                 Protein_Name = table_EC_Non_Listed_Proteins$Protein_Name,
                                 Mass_kDa = table_EC_Non_Listed_Proteins$Mass_kDa,
                                 Length_AA = table_EC_Non_Listed_Proteins$Length_AA,
                                 RBP2GO_Score = table_EC_Non_Listed_Proteins$RBP2GO_Score,
                                 Nb_Datasets = table_EC_Non_Listed_Proteins$Nb_Datasets,
                                 Listing_Count = table_EC_Non_Listed_Proteins$Listing_Count,
                                 AVG10_Int_Listing_Count = table_EC_Non_Listed_Proteins$AVG10_Int_Listing_Count
                      )
)
setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
save(table_RBP2GO,file = "table_RBP2GO.rda")


##################################################################################################
# Update the df_HS_UniRef50_All table
setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
load("table_RBP2GO.rda")
load("df_HS_UniRef50_All.rda")

# Only keep the UniProt IDs which are "TRUE" UniProt IDs
df_HS_UniRef50_All <- df_HS_UniRef50_All[df_HS_UniRef50_All$Cluster_Members %in% table_RBP2GO$Uniprot_ID,]

Entry_Name <- rep("",length(df_HS_UniRef50_All$UniProt_ID))
Protein_Name <- rep("",length(df_HS_UniRef50_All$UniProt_ID))
Mass_kDa <- rep(0,length(df_HS_UniRef50_All$UniProt_ID))
Length_AA <- rep(0,length(df_HS_UniRef50_All$UniProt_ID))
RBP2GO_Score <- rep("",length(df_HS_UniRef50_All$UniProt_ID))
Nb_Datasets <- rep("",length(df_HS_UniRef50_All$UniProt_ID))
Listing_Count <- rep("",length(df_HS_UniRef50_All$UniProt_ID))
AVG10_Int_Listing_Count <- rep("",length(df_HS_UniRef50_All$UniProt_ID))
# Add the missing information into the table
for (x in 1:dim(df_HS_UniRef50_All)[1]) {
  #for (x in 949:950) {
  Cluster_Member <- df_HS_UniRef50_All[x,]$Cluster_Members
  Entry_Name[x] <- as.character(table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$Entry_Name)
  Protein_Name[x] <- as.character(table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$Protein_Name)
  Mass_kDa[x] <- table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$Mass_kDa
  Length_AA[x] <- table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$Length_AA
  RBP2GO_Score[x] <- table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$RBP2GO_Score
  Nb_Datasets[x] <- table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$Nb_Datasets
  Listing_Count[x] <- as.character(table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$Listing_Count)
  AVG10_Int_Listing_Count[x] <- as.character(table_RBP2GO[table_RBP2GO$Uniprot_ID == Cluster_Member,]$AVG10_Int_Listing_Count)
}

df_HS_UniRef50_All$Entry_Name <- Entry_Name
df_HS_UniRef50_All$Protein_Name <- Protein_Name
df_HS_UniRef50_All$Mass_kDa <- Mass_kDa
df_HS_UniRef50_All$Length_AA <- Length_AA
df_HS_UniRef50_All$RBP2GO_Score <- RBP2GO_Score
df_HS_UniRef50_All$Nb_Datasets <- Nb_Datasets
df_HS_UniRef50_All$Listing_Count <- Listing_Count
df_HS_UniRef50_All$AVG10_Int_Listing_Count <- AVG10_Int_Listing_Count

save(df_HS_UniRef50_All,file = "df_HS_UniRef50_All.rda")

setwd("~/Desktop/RBP2GO/RBP2GO_Orthologs/")
load("df_HS_UniRef50_All.rda")
load("table_RBP2GO.rda")

