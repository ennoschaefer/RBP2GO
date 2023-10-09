import pandas as pd
import pyreadr
import numpy as np
import time
import os
import subprocess

t0 = time.time()

# Enter Filepath -> More needed further down
os.chdir("")
# *****************************************************************************************
##########################################################################################
# *****************************************************************************************
"""
How to use this code:
1. Enter Filepaths of relevant data: 
    a) Old Database of both the RBP and the non-RBP [.RDS]
    b) New Datasets to be implemented [.RDS][.csv]
2. Use the "datamerger" to merge the new datasets onto the old database 
(beware, that the column to be merged on needs unique identifiers, eg. "Uniprot_ID" 
-> Gene Name may contain duplicates in the Database)
3. Use the "sort_columns" to put the new datasets into alignment with the old ones already 
present in the database. Oldest studies will appear as the first, the newest study will be the last column
4. Use "recalculate_scores" to properly calculate the average 10 interaction partner listing scores, 
the RBP2GO score and the RBP2GO composite score
5. Export the new database as .RDS and implement the new tables onto the server
6. To fix an issue with the pyreadr package, an additional R script is called (update_RDS_db_files.R)

The final database files can be found in the "R_export" folder in your current working-directory.
"""


# *****************************************************************************************
##########################################################################################
# *****************************************************************************************
print("Loading functions...")


# Functions to transform the data
def RDS_reader(filepath):
    """
    Reads a .RDS file and exports it as a pandas dataframe
    filepath: Complete filepath + filename of the file of interest
    """
    _importfile = pyreadr.read_r(filepath)
    _exportfile = _importfile[None]
    return _exportfile


def csv_reader(filepath, separator=";"):
    """
    Reads a .csv and exports it as a pandas dataframe
    """
    _import_file = pd.read_csv(filepath, sep=separator)
    return _import_file


def datamerger(dataset, db_file, db_file_nonRBPs):
    """
    Merges new dataframe with existing database-file
    dataset: pandas dataframe;  The new study that wants to be integrated into the database
    db_file: pandas dataframe; File of the database of the corresponding species
    db_file_nonRBPs: pandas dataframe; File of non RBPs
    column_to_merge: string; Column-name of the column in the db_file of the corresponding data in the dataset
    new_RBPs: pandas dataframe; Table of all new RBPs found throughout the studies
    """
    new_df_db = None
    new_df_db_nonRBPs = db_file_nonRBPs.copy()
    new_RBP_buffer = None
    unknown_proteins_df = None
    column_to_merge = dataset.columns[0]

    dataset = dataset.drop_duplicates(dataset.columns[0])

    new_df_db = pd.merge(db_file, dataset, on=column_to_merge, how="left")
    col_list = new_df_db.columns.tolist()
    col_list_ordered = col_list[: col_list.index("GO_BP")]
    col_list_ordered.append(col_list[-1])
    col_list_ordered.extend(col_list[col_list.index("GO_BP") : -1])
    new_df_db = new_df_db[col_list_ordered]

    # Get all RBPs that were not in the RBP dataframe -> Will be returned at the end of the function
    new_RBP_buffer = pd.merge(
        dataset, db_file, on=column_to_merge, how="outer", indicator=True
    )
    new_RBP_buffer = new_RBP_buffer[new_RBP_buffer["_merge"] == "left_only"]
    new_RBP_buffer = new_RBP_buffer.iloc[:, :2].reset_index(drop=True)

    # Add the new RBPs to the existing database using the Non_RBP_Table
    if len(new_RBP_buffer) > 0:
        buffer_df = pd.merge(
            db_file_nonRBPs,
            new_RBP_buffer,
            on=column_to_merge,
            how="outer",
            indicator=True,
        )
        unknown_proteins_df = buffer_df[buffer_df["_merge"] == "right_only"]
        if buffer_df["_merge"].str.count("both").sum() > 0:
            buffer_df_both = buffer_df[buffer_df["_merge"] == "both"]
            new_df_db = pd.concat([new_df_db, buffer_df_both.iloc[:, :-1]]).reset_index(
                drop=True
            )
            new_df_db_nonRBPs = buffer_df[buffer_df["_merge"] == "left_only"]
            new_df_db_nonRBPs = new_df_db_nonRBPs.iloc[:, :-2].reset_index(drop=True)

        if len(unknown_proteins_df) > 0:
            raise ValueError(
                "There are unknown proteins in " + unknown_proteins_df.columns[-2]
            )
            # print("Unknown proteins in " + unknown_proteins_df.columns[-2])       # activate when needed and deactivate the ValueError above.
    if len(db_file) + len(db_file_nonRBPs) != len(new_df_db) + len(new_df_db_nonRBPs):
        if len(db_file) + len(db_file_nonRBPs) > len(new_df_db) + len(
            new_df_db_nonRBPs
        ):
            raise ValueError("Entries have been deleted from the database!")
        if len(db_file) + len(db_file_nonRBPs) < len(new_df_db) + len(
            new_df_db_nonRBPs
        ):
            raise ValueError("There are more entries in the database than before!")

    return new_df_db, new_df_db_nonRBPs


def recalculate_scores(db_file, db_file_nonRBPs):
    """
    Recalculates the following scores from a RBP2GO database-file:
    Listing_Count, Nb_Datasets, AVG10_Int_Listing_Count, RBP2GO_Score, Composite_Score.
    If there are no String interaction partners given in the database-files,
    the AVG10_Int_Listing_Count as well as the RBP2GO_Score and the Composite_Score will not be calculated.

    Takes both the database-file of RBPs and the database-file of the non-RBPs in pandas.DataFrame format
    -> Faster and easier calculation with both dataframes at once.
    Returns both database-files in pandas.DataFrame format.
    """

    index_first_ds = db_file.columns.get_loc("Only without poly(A) enrichment") + 1
    index_last_ds = db_file.columns.get_loc("GO_BP")
    total_ds_number = index_last_ds - index_first_ds

    db_file["Nb_Datasets"] = total_ds_number
    db_file_nonRBPs["Nb_Datasets"] = total_ds_number

    # Replace the listing count
    db_file["Listing_Count"] = (
        db_file.iloc[:, index_first_ds:index_last_ds]
        .replace("X", 1)
        .replace("", np.nan)
        .replace(" ", np.nan)
        .sum(1)
    )

    # Replace the AVG10 interaction partners score (this takes a little while)
    # RBP-File
    if "String_PPI" in db_file.columns:
        string_ppi_df = db_file[["Uniprot_ID", "String_PPI"]]
        string_ppi_df["String_PPI"] = string_ppi_df["String_PPI"].str.split("; ")
        string_ppi_df = string_ppi_df.rename(
            columns={"Uniprot_ID": "Uniprot_ID_ori", "String_PPI": "Uniprot_ID"}
        ).explode("Uniprot_ID")
        string_ppi_df = string_ppi_df[~(string_ppi_df["Uniprot_ID"] == "")]

        string_ppi_merged_df = pd.merge(
            string_ppi_df,
            db_file[["Uniprot_ID", "Listing_Count"]],
            on="Uniprot_ID",
            how="left",
        )
        string_ppi_merged_df["occurence"] = (
            string_ppi_merged_df.groupby("Uniprot_ID_ori").cumcount() + 1
        )
        string_ppi_merged_df["Listing_Count"] = string_ppi_merged_df[
            "Listing_Count"
        ].replace(np.nan, 0)

        string_ppi_pivot_df = string_ppi_merged_df.pivot(
            values="Listing_Count", columns="occurence", index="Uniprot_ID_ori"
        )
        string_ppi_pivot_df = string_ppi_pivot_df.iloc[:, :10]
        string_ppi_pivot_df["NA_count"] = string_ppi_pivot_df.isna().sum(1)

        string_ppi_pivot_df["AVG10_Int_Listing_Count_new"] = round(
            string_ppi_pivot_df.iloc[:, :10].sum(1)
            / (10 - string_ppi_pivot_df["NA_count"]),
            1,
        )

        string_ppi_avg10_df = (
            string_ppi_pivot_df["AVG10_Int_Listing_Count_new"]
            .reset_index(drop=False)
            .rename(columns={"Uniprot_ID_ori": "Uniprot_ID"})
        )

        db_file["AVG10_Int_Listing_Count"] = pd.merge(
            db_file, string_ppi_avg10_df, on="Uniprot_ID", how="left"
        )["AVG10_Int_Listing_Count_new"].fillna(0)

        ## Non RBP File
        string_ppi_df = db_file_nonRBPs[["Uniprot_ID", "String_PPI"]]
        string_ppi_df["String_PPI"] = string_ppi_df["String_PPI"].str.split("; ")
        string_ppi_df = string_ppi_df.rename(
            columns={"Uniprot_ID": "Uniprot_ID_ori", "String_PPI": "Uniprot_ID"}
        ).explode("Uniprot_ID")
        string_ppi_df = string_ppi_df[~(string_ppi_df["Uniprot_ID"] == "")]

        string_ppi_merged_df = pd.merge(
            string_ppi_df,
            db_file[["Uniprot_ID", "Listing_Count"]],
            on="Uniprot_ID",
            how="left",
        )
        string_ppi_merged_df["occurence"] = (
            string_ppi_merged_df.groupby("Uniprot_ID_ori").cumcount() + 1
        )
        string_ppi_merged_df["Listing_Count"] = string_ppi_merged_df[
            "Listing_Count"
        ].replace(np.nan, 0)

        string_ppi_pivot_df = string_ppi_merged_df.pivot(
            values="Listing_Count", columns="occurence", index="Uniprot_ID_ori"
        )
        string_ppi_pivot_df = string_ppi_pivot_df.iloc[:, :10]
        string_ppi_pivot_df["NA_count"] = string_ppi_pivot_df.isna().sum(1)

        string_ppi_pivot_df["AVG10_Int_Listing_Count_new"] = round(
            string_ppi_pivot_df.iloc[:, :10].sum(1)
            / (10 - string_ppi_pivot_df["NA_count"]),
            1,
        )

        string_ppi_avg10_df = (
            string_ppi_pivot_df["AVG10_Int_Listing_Count_new"]
            .reset_index(drop=False)
            .rename(columns={"Uniprot_ID_ori": "Uniprot_ID"})
        )

        db_file_nonRBPs["AVG10_Int_Listing_Count"] = pd.merge(
            db_file_nonRBPs, string_ppi_avg10_df, on="Uniprot_ID", how="left"
        )["AVG10_Int_Listing_Count_new"].fillna(0)

        # Recalculate the RBP2GO_Score
        db_file["RBP2GO_Score"] = (
            (db_file["Listing_Count"] / db_file["Nb_Datasets"]) * 50
        ) + ((db_file["AVG10_Int_Listing_Count"] / db_file["Nb_Datasets"]) * 50)

        db_file_nonRBPs["RBP2GO_Score"] = (
            (db_file_nonRBPs["Listing_Count"] / db_file_nonRBPs["Nb_Datasets"]) * 50
        ) + (
            (
                db_file_nonRBPs["AVG10_Int_Listing_Count"]
                / db_file_nonRBPs["Nb_Datasets"]
            )
            * 50
        )

        # Recalculate the RBP2GO Composite Score
        db_file["Composite_Score"] = (
            ((db_file["Listing_Count"] / db_file["Nb_Datasets"]) * 50)
            + ((db_file["AVG10_Int_Listing_Count"] / db_file["Nb_Datasets"]) * 25)
            + db_file["Domain_score"]
        )
        db_file_nonRBPs["Composite_Score"] = (
            ((db_file_nonRBPs["Listing_Count"] / db_file_nonRBPs["Nb_Datasets"]) * 50)
            + (
                (
                    db_file_nonRBPs["AVG10_Int_Listing_Count"]
                    / db_file_nonRBPs["Nb_Datasets"]
                )
                * 25
            )
            + db_file_nonRBPs["Domain_score"]
        )

        # Sort Rows by the highest RBP2GO score
        db_file.sort_values(
            "Composite_Score", ascending=False, inplace=True, ignore_index=True
        )
        db_file_nonRBPs.sort_values(
            "Composite_Score", ascending=False, inplace=True, ignore_index=True
        )

        # Fill NaN again with white space
        db_file = db_file.fillna("")
        db_file_nonRBPs = db_file_nonRBPs.fillna("")
    return db_file, db_file_nonRBPs


def sort_columns(db_file):
    """
    Takes a RBP2GO database file and resorts the columns of the different studies by year. Most columns are named "Author_Cell-line_year".
    Studies without a year will be listed as first columns after the column "Only without poly(A) enrichment".
    Tales a database-file in pandas.DataFrame format.

    Returns database-file in pandas.DataFrame format with resorted columns.
    """
    index_first_ds = db_file.columns.get_loc("Only without poly(A) enrichment") + 1
    index_last_ds = db_file.columns.get_loc("GO_BP")

    old_column_sorting = db_file.columns.tolist()
    new_column_sorting = old_column_sorting[:(index_first_ds)]

    needs_sorting = sorted(old_column_sorting[index_first_ds:index_last_ds])
    is_sorted = sorted(
        needs_sorting,
        key=lambda x: (
            len(x.rsplit("_", 1)) == 1,
            int(x.rsplit("_", 1)[-1])
            if x.rsplit("_", 1)[-1].isdigit()
            else float("-inf"),
        ),
    )

    new_column_sorting.extend(is_sorted)
    new_column_sorting.extend(old_column_sorting[index_last_ds:])

    sorted_db_file = db_file[new_column_sorting]
    return sorted_db_file


# *****************************************************************************************
##########################################################################################
# *****************************************************************************************
print("Loading database and dataset files...")
# Relevant File-Paths
filepath_ds = ""
filepath_db = ""

ori_db_HS_path = "table_HS_Dataset.RDS"
ori_db_HS_nonRBP_path = "table_HS_Non_Listed_Proteins.RDS"

ori_db_SaE_path = "table_SaE_Dataset.RDS"
ori_db_SaE_nonRBP_path = "table_SaE_Non_Listed_Proteins.RDS"

ori_db_EC_path = "table_EC_Dataset.RDS"
ori_db_EC_nonRBP_path = "table_EC_Non_Listed_Proteins.RDS"

ori_db_MM_path = "table_MM_Dataset.RDS"
ori_db_MM_nonRBP_path = "table_mm_Non_Listed_Proteins.RDS"

ori_db_AT_path = "table_AT_Dataset.RDS"
ori_db_AT_nonRBP_path = "table_AT_Non_Listed_Proteins.RDS"

ori_db_LM_path = "table_LM_Dataset.RDS"
ori_db_LM_nonRBP_path = "table_LM_Non_Listed_Proteins.RDS"

ori_db_SC_path = "table_SC_Dataset.RDS"
ori_db_SC_nonRBP_path = "table_SC_Non_Listed_Proteins.RDS"

ori_db_PF_path = "table_PF_Dataset.RDS"  # doesn't get updated, but needs manual changes
ori_db_PF_nonRBP_path = "table_PF_Non_Listed_Proteins.RDS"

# Paths of datasets
gerovac_ds_name = "Original_Data//1.Gerovac_SaE.csv"
bach_pages_ds_name = "Original_Data//3.Bach-Pages_ST1_AT.csv"
qin_ds_name = "Original_Data//5.Qin_Human.csv"
kamel_ds_name = "Original_Data//6.Kamel_Human.csv"
dvir_ds_name = "Original_Data//8.Dvir_ST1_Human_ESC.csv"
kalesh_ds_name = "Original_Data//9.Kalesh_LM.csv"
vieira_ds_name = "Original_Data//10.Vieira-Vieira_Human.csv"
bressin_ds_HS_name = "Original_Data//11.Bressin_ST13_Human.csv"
# bressin_ds_SaE_name = "Original_Data//11.Bressin_ST14_Salmonella.csv"
bressin_ds_EC_name = "Original_Data//11.Bressin_ST15_Ecoli.csv"
gosh_ds_HS_name = "Original_Data//14.Gosh_ST1_Human.csv"
mestre_farras_HS_name = "Original_Data//15.Mestre-Farras_melanoma_Human.csv"
stenum_df_EC_name = "Original_Data//17.Stenum_Ecoli.csv"
perez_perri_ds_MM_brain_Name = "Original_Data//19.Perrez-Perri_Mouse_Brain.csv"
perez_perri_ds_MM_kidney_Name = "Original_Data//19.Perrez-Perri_Mouse_Kidney.csv"
perez_perri_ds_MM_liver_Name = "Original_Data//19.Perrez-Perri_Mouse_Liver.csv"
malhotra_ds_HS_name = "Original_Data//21.Malhotra_Human.csv"
asencio_df_SC_name = "Original_Data//22.Asencio_SC.csv"
schmidt_df_HS_name = "Original_Data//24.Schmidt_Human.csv"
rajagopal_df_HS_name = "Original_Data//25.Rajagopal_R-DeeP_HS.csv"
whitworth_df_MM_name = "Original_Data//26.Whitworth_MM_Liver.csv"
whitworth_df_HS_XG77_name = "Original_Data//26.Whitworth_Human_XG77.csv"
whitworth_df_HS_XG147_name = "Original_Data//26.Whitworth_Human_XG147.csv"
bao_df_HS_name = "Original_Data//27.Bao_Human.csv"


# Read original database-files to dataframes
ori_db_HS_df = RDS_reader(filepath_db + ori_db_HS_path)
ori_db_HS_nonRBP_df = RDS_reader(filepath_db + ori_db_HS_nonRBP_path)

ori_db_SaE_df = RDS_reader(filepath_db + ori_db_SaE_path)
ori_db_SaE_nonRBP_df = RDS_reader(filepath_db + ori_db_SaE_nonRBP_path)

ori_db_EC_df = RDS_reader(filepath_db + ori_db_EC_path)
ori_db_EC_nonRBP_df = RDS_reader(filepath_db + ori_db_EC_nonRBP_path)

ori_db_MM_df = RDS_reader(filepath_db + ori_db_MM_path)
ori_db_MM_nonRBP_df = RDS_reader(filepath_db + ori_db_MM_nonRBP_path)

ori_db_AT_df = RDS_reader(filepath_db + ori_db_AT_path)
ori_db_AT_nonRBP_df = RDS_reader(filepath_db + ori_db_AT_nonRBP_path)

ori_db_LM_df = RDS_reader(filepath_db + ori_db_LM_path)
ori_db_LM_nonRBP_df = RDS_reader(filepath_db + ori_db_LM_nonRBP_path)

ori_db_SC_df = RDS_reader(filepath_db + ori_db_SC_path)
ori_db_SC_nonRBP_df = RDS_reader(filepath_db + ori_db_SC_nonRBP_path)

ori_db_PF_df = RDS_reader(filepath_db + ori_db_PF_path)
ori_db_PF_nonRBP_df = RDS_reader(filepath_db + ori_db_PF_nonRBP_path)

# Read new datasets to dataframes
gerovac_SaE_df = csv_reader(filepath_ds + gerovac_ds_name)
qin_HS_df = csv_reader(filepath_ds + qin_ds_name)
kamel_HS_df = csv_reader(filepath_ds + kamel_ds_name)
bach_pages_AT_df = csv_reader(filepath_ds + bach_pages_ds_name)
dvir_HS_df = csv_reader(filepath_ds + dvir_ds_name)
kalesh_LM_df = csv_reader(filepath_ds + kalesh_ds_name)
vieira_HS_df = csv_reader(filepath_ds + vieira_ds_name)
bressin_HS_df = csv_reader(filepath_ds + bressin_ds_HS_name)
# bressin_SaE_df = csv_reader(filepath_ds + bressin_ds_SaE_name) # implement when dataset is ready again. > See docu
bressin_EC_df = csv_reader(filepath_ds + bressin_ds_EC_name)
gosh_HS_df = csv_reader(filepath_ds + gosh_ds_HS_name)
mestre_farras_HS_df = csv_reader(filepath_ds + mestre_farras_HS_name)
stenum_EC_df = csv_reader(filepath_ds + stenum_df_EC_name)
perez_perri_MM_brain_df = csv_reader(filepath_ds + perez_perri_ds_MM_brain_Name)
perez_perri_MM_kidney_df = csv_reader(filepath_ds + perez_perri_ds_MM_kidney_Name)
perez_perri_MM_liver_df = csv_reader(filepath_ds + perez_perri_ds_MM_liver_Name)
malhotra_HS_df = csv_reader(filepath_ds + malhotra_ds_HS_name)
asencio_SC_df = csv_reader(filepath_ds + asencio_df_SC_name)
schmidt_HS_df = csv_reader(filepath_ds + schmidt_df_HS_name)
rajagopal_HS_df = csv_reader(filepath_ds + rajagopal_df_HS_name)
whitworth_MM_liver_df = csv_reader(filepath_ds + whitworth_df_MM_name)
whitworth_HS_XG77_df = csv_reader(filepath_ds + whitworth_df_HS_XG77_name)
whitworth_HS_XG147_df = csv_reader(filepath_ds + whitworth_df_HS_XG147_name)
bao_HS_df = csv_reader(filepath_ds + bao_df_HS_name)

# Merge the datasets with the respective database sets
print("Updating database with dataset files...")
## Homo Sapiens
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    qin_HS_df, ori_db_HS_df, ori_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    kamel_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    dvir_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    vieira_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    bressin_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    gosh_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    mestre_farras_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    malhotra_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    schmidt_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    rajagopal_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    whitworth_HS_XG77_df, new_db_HS_df, new_db_HS_nonRBP_df
)
new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    whitworth_HS_XG147_df, new_db_HS_df, new_db_HS_nonRBP_df
)

# Need to change column names before adding this dataset to the database
# Will do this manually, because its a one-time thing.
new_db_HS_df = new_db_HS_df.rename(columns={"Bao_HeLa_2018": "Bao_HeLa_RICK_2018"})

new_db_HS_df, new_db_HS_nonRBP_df = datamerger(
    bao_HS_df, new_db_HS_df, new_db_HS_nonRBP_df
)


## Arabidopsis
new_db_AT_df, new_db_AT_nonRBP_df = datamerger(
    bach_pages_AT_df, ori_db_AT_df, ori_db_AT_nonRBP_df
)

## Salmonella
new_db_SaE_df, new_db_SaE_nonRBP_df = datamerger(
    gerovac_SaE_df, ori_db_SaE_df, ori_db_SaE_nonRBP_df
)
# new_db_SaE_df, new_db_SaE_nonRBP_df = datamerger(
#     bressin_SaE_df, new_db_SaE_df, new_db_SaE_nonRBP_df
# )

## Leishmania M.
new_db_LM_df, new_db_LM_nonRBP_df = datamerger(
    kalesh_LM_df, ori_db_LM_df, ori_db_LM_nonRBP_df
)

## E.coli
new_db_EC_df, new_db_EC_nonRBP_df = datamerger(
    bressin_EC_df, ori_db_EC_df, ori_db_EC_nonRBP_df
)
new_db_EC_df, new_db_EC_nonRBP_df = datamerger(
    stenum_EC_df, new_db_EC_df, new_db_EC_nonRBP_df
)

## S. cerevisiae
new_db_SC_df, new_db_SC_nonRBP_df = datamerger(
    asencio_SC_df, ori_db_SC_df, ori_db_SC_nonRBP_df
)

## Mus musculus
new_db_MM_df, new_db_MM_nonRBP_df = datamerger(
    perez_perri_MM_brain_df, ori_db_MM_df, ori_db_MM_nonRBP_df
)
new_db_MM_df, new_db_MM_nonRBP_df = datamerger(
    perez_perri_MM_kidney_df, new_db_MM_df, new_db_MM_nonRBP_df
)
new_db_MM_df, new_db_MM_nonRBP_df = datamerger(
    perez_perri_MM_liver_df, new_db_MM_df, new_db_MM_nonRBP_df
)
new_db_MM_df, new_db_MM_nonRBP_df = datamerger(
    whitworth_MM_liver_df, new_db_MM_df, new_db_MM_nonRBP_df
)

# Sorting the new database-table columns by year of study
print("Resorting the columns of the database files...")
new_db_HS_df_sorted = sort_columns(new_db_HS_df)
new_db_AT_df_sorted = sort_columns(new_db_AT_df)
new_db_SaE_df_sorted = sort_columns(new_db_SaE_df)
new_db_LM_df_sorted = sort_columns(new_db_LM_df)
new_db_EC_df_sorted = sort_columns(new_db_EC_df)
new_db_SC_df_sorted = sort_columns(new_db_SC_df)
new_db_MM_df_sorted = sort_columns(new_db_MM_df)


# Recalculating the scores in all database tables
print("Recalculating the scores of the database files...")
new_db_HS_df_RBP2GO, new_db_HS_nonRBP_df_RBP2GO = recalculate_scores(
    new_db_HS_df_sorted, new_db_HS_nonRBP_df
)
new_db_AT_df_RBP2GO, new_db_AT_nonRBP_df_RBP2GO = recalculate_scores(
    new_db_AT_df_sorted, new_db_AT_nonRBP_df
)
new_db_SaE_df_RBP2GO, new_db_SaE_nonRBP_df_RBP2GO = recalculate_scores(
    new_db_SaE_df_sorted, new_db_SaE_nonRBP_df
)
new_db_LM_df_RBP2GO, new_db_LM_nonRBP_df_RBP2GO = recalculate_scores(
    new_db_LM_df_sorted, new_db_LM_nonRBP_df
)
new_db_EC_df_RBP2GO, new_db_EC_nonRBP_df_RBP2GO = recalculate_scores(
    new_db_EC_df_sorted, new_db_EC_nonRBP_df
)
new_db_SC_df_RBP2GO, new_db_SC_nonRBP_df_RBP2GO = recalculate_scores(
    new_db_SC_df_sorted, new_db_SC_nonRBP_df
)
new_db_MM_df_RBP2GO, new_db_MM_nonRBP_df_RBP2GO = recalculate_scores(
    new_db_MM_df_sorted, new_db_MM_nonRBP_df
)

# Some manual changes on the datasets > Remove when done!
new_db_HS_df[new_db_HS_df.Entry_Name == "TITIN_MOUSE"]["Mass_kDa"] = 3816.030
new_db_HS_df[new_db_HS_df.Entry_Name == "SYNE1_HUMAN"]["Mass_kDa"] = 1011.086
new_db_MM_df[new_db_MM_df.Entry_Name == "TITIN_MOUSE"]["Mass_kDa"] = 3906.488
ori_db_PF_df[ori_db_PF_df.Entry_Name == "MLRR1_PLAF7"]["Mass_kDa"] = 1187.585

# Timeevaluation before the export of the data
t1 = time.time()
total = t1 - t0
print("This took " + str(round(total, 2)) + " seconds, before the export.")

# Export the data to .RDS
print("Starting the export of the database files...")
if os.path.isdir("Py_export"):
    print("Exportfolder ready.")
else:
    try:
        os.mkdir("Py_export")
    except OSError:
        path = os.getcwd()
        raise LookupError(
            """Exportfolder doesn't exist and can't be created! 
            Please create exportfolder ('Py_export') manually into the current directory. 
            The current working directory is %s"""
            % path
        )
    else:
        print("Successfully created exportfolder in directory.")


## RBP data
pyreadr.write_rdata("Py_export//table_HS_Dataset.RDS", new_db_HS_df_RBP2GO)
pyreadr.write_rdata("Py_export//table_AT_Dataset.RDS", new_db_AT_df_RBP2GO)
pyreadr.write_rdata("Py_export//table_SaE_Dataset.RDS", new_db_SaE_df_RBP2GO)
pyreadr.write_rdata("Py_export//table_LM_Dataset.RDS", new_db_LM_df_RBP2GO)
pyreadr.write_rdata("Py_export//table_EC_Dataset.RDS", new_db_EC_df_RBP2GO)
pyreadr.write_rdata("Py_export//table_SC_Dataset.RDS", new_db_SC_df_RBP2GO)
pyreadr.write_rdata("Py_export//table_MM_Dataset.RDS", new_db_MM_df_RBP2GO)
pyreadr.write_rdata("Py_export//table_PF_Dataset.RDS", ori_db_PF_df)

## Non-RBP data
pyreadr.write_rdata(
    "Py_export//table_HS_Non_Listed_Proteins.RDS", new_db_HS_nonRBP_df_RBP2GO
)
pyreadr.write_rdata(
    "Py_export//table_AT_Non_Listed_Proteins.RDS", new_db_AT_nonRBP_df_RBP2GO
)
pyreadr.write_rdata(
    "Py_export//table_SaE_Non_Listed_Proteins.RDS", new_db_SaE_nonRBP_df_RBP2GO
)
pyreadr.write_rdata(
    "Py_export//table_LM_Non_Listed_Proteins.RDS", new_db_LM_nonRBP_df_RBP2GO
)
pyreadr.write_rdata(
    "Py_export//table_EC_Non_Listed_Proteins.RDS", new_db_EC_nonRBP_df_RBP2GO
)
pyreadr.write_rdata(
    "Py_export//table_SC_Non_Listed_Proteins.RDS", new_db_SC_nonRBP_df_RBP2GO
)
pyreadr.write_rdata(
    "Py_export//table_MM_Non_Listed_Proteins.RDS", new_db_MM_nonRBP_df_RBP2GO
)
pyreadr.write_rdata("Py_export//table_PF_Non_Listed_Proteins.RDS", ori_db_PF_nonRBP_df)

# Because of an export issue of the package pyreadr, an R script is called, that fixes the export data.
print("Starting extra R file for bug prevention...")

if os.path.isdir("R_export"):
    print("R folder ready.")
else:
    try:
        os.mkdir("R_export")
    except OSError:
        path = os.getcwd()
        raise LookupError(
            """Exportfolder doesn't exist and can't be created! 
            Please create exportfolder ('R_export') manually into the current directory. 
            The current working directory is %s"""
            % path
        )
    else:
        print("Successfully created R folder in directory.")


subprocess.call(
    [
        "Rscript",
        "--vanilla",
        "update_RDS_db_files.R",
    ],
    shell=False,
)

# Lets take a look, how long all of this takes...
t1 = time.time()
total = t1 - t0
print("This took " + str(round(total, 2)) + " seconds in total.")
print("Done!")
