###########################################################################################################
#### Prepare the alias table including RBP and non-RBP ####

####################
# For all species  #
####################

species_list <- c("HS", "MM", "SC", "DM", "AT", "CE", "DaR", "EC", "LD", "LM", "PF", "SaE","TB")

for (species in species_list){
  print(paste0("Alias for ",species,""))
  table_alias <- NULL

  # Get and combine the datasets
  data_NLP <- readRDS(paste0("~/Desktop/230809_Database_RBP2GO_V2.0/data/table_",species,"_Non_Listed_Proteins.RDS"))
  data_RBP <- readRDS(paste0("~/Desktop/230809_Database_RBP2GO_V2.0/data/table_",species,"_Dataset.RDS"))
  table_alias <- rbind(data_RBP[,c("Entry_Name","Uniprot_ID","Protein_Name","Gene_Name","Alias_Names")],
                       data_NLP[,c("Entry_Name","Uniprot_ID","Protein_Name","Gene_Name","Alias_Names")])
  
  # Create the first line of the alias dataset
  df_alias <- NULL
  df_alias <- data.frame(
    Entry_Name = table_alias$Entry_Name[1],
    alias = c(table_alias$Entry_Name[1],
              table_alias$Uniprot_ID[1],
              table_alias$Protein_Name[1],
              table_alias$Gene_Name[1],
              unlist(strsplit(table_alias$Alias_Names[1], " "))
              ))
  
  # Complete the table for all entries
 for (i in 2:length(table_alias$Entry_Name)) {
#for (i in 2:100) {
   #print(i)
   df_alias_i <- NULL
   df_alias_i <- data.frame(Entry_Name = table_alias$Entry_Name[i],
                               alias = c(table_alias$Entry_Name[i],
                                         table_alias$Uniprot_ID[i],
                                         table_alias$Protein_Name[i],
                                         table_alias$Gene_Name[i],
                                         unlist(strsplit(table_alias$Alias_Names[i], " "))
                                         ))
   df_alias <- rbind(df_alias,df_alias_i)
   #print(dim(df_alias))
 }

  # remove duplicates
  df_alias <- unique(df_alias)
  # Save files
  saveRDS(df_alias, paste0("~/Desktop/Data_alias/df_",species,"_alias.RDS"))
}



####################
# For Homo sapiens #
####################

# RBP and non-RBP tables should already exist - see above
table_HS_Non_Listed_Proteins <- readRDS("~/Desktop/230809_Database_RBP2GO_V2.0/data/table_HS_Non_Listed_Proteins.RDS")
table_HS_Dataset <- readRDS("~/Desktop/230809_Database_RBP2GO_V2.0/data/table_HS_Dataset.RDS")

# Combine the two tables
table_HS_alias <- rbind(table_HS_Dataset[,c("Entry_Name","Uniprot_ID","Protein_Name","Gene_Name","Alias_Names")],table_HS_Non_Listed_Proteins[,c("Entry_Name","Uniprot_ID","Protein_Name","Gene_Name","Alias_Names")])

# Take the UniProt ID as reference ID
# Get the table with alias for each Uniprot ID - they are unique.
# First protein
Uniprot_ID_HS <- table_HS_alias$Uniprot_ID[1]
# Get the corresponding Entry Name
prot_EntryName <- table_HS_alias$Entry_Name[1]
# Get the corresponding Protein Name
prot_Prot_Name <- table_HS_alias$Protein_Name[1]
# Get the corresponding Gene Name
prot_Gene_Name <- table_HS_alias$Gene_Name[1]
# Get the corresponding Alias (possibly several)
prot_Alias_List <- table_HS_alias$Alias_Names[1]
prot_Alias_List <- unlist(strsplit(prot_Alias_List, " "))
# Make the dataframe
df_HS_alias <- data.frame(Entry_Name = prot_EntryName, alias = c(prot_EntryName,Uniprot_ID_HS,prot_Prot_Name,prot_Gene_Name,prot_Alias_List))

# Complete the table for all entries
#for (i in 2:3) {  # line to test the following lines
df_HS_alias_i <- NULL
for (i in 2:length(table_HS_alias$Entry_Name)) {
  Uniprot_ID_HS_i <- table_HS_alias$Uniprot_ID[i]
  prot_EntryName_i <- table_HS_alias$Entry_Name[i]
  prot_Prot_Name_i <- table_HS_alias$Protein_Name[i]
  prot_Gene_Name_i <- table_HS_alias$Gene_Name[i]
  prot_Alias_List_i <- unlist(strsplit(table_HS_alias$Alias_Names[i], " "))
  df_HS_alias_i <- data.frame(Entry_Name = prot_EntryName_i, alias = c(prot_EntryName_i,Uniprot_ID_HS_i,prot_Prot_Name_i,prot_Gene_Name_i,prot_Alias_List_i))
  df_HS_alias <- rbind(df_HS_alias,df_HS_alias_i)
}
# remove duplicates
df_HS_alias <- unique(df_HS_alias)


# Save files
write.table(df_HS_alias, file="df_HS_alias.csv", row.names = FALSE, col.names = TRUE, sep = ";")
saveRDS(df_HS_alias, file = "df_HS_alias.RDS")


########
for (species in species_list){
  table <- readRDS(paste0("~/Desktop/Data_alias/df_",species,"_alias.RDS"))
  print(c(species," ",length(unique(table$Entry_Name))))
  }



