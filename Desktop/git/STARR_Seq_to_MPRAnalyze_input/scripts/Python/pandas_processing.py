# File: pandas_processing.py
# AUTHOR: CHRIS HSU <chrhsu@uw.edu>
# CREATED DATE: APRIL 12, 2021

import os
import time
import argparse
import pandas as pd

"""
This script is to process the files output from the pipeline.sh script.
-Input: 09_UMI_group directory
-Output: pipeline/12_pandas_processed/DNA_RNA_count.csv file containing DNA and RNA counts

I UNDERSTAND THAT MANY OF THESE FUNCTIONS ARE INEFFICIENT AND WORDY
IF THERE'S A BETTER WAY, PLEASE LET ME KNOW!
"""


# Parser Commands # ------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Inputs/Output for processing UMI files.')
requiredNamed = parser.add_argument_group('Required named arguments')
requiredNamed.add_argument('-d', '--directory', type=str, help='Directory to pipeline_input', required=True)
requiredNamed.add_argument('-o', '--outputdir', type=str, help='Output Directory)', required=True)
args = parser.parse_args()

input_dir = args.directory
output_dir = args.outputdir

# Default locations
# input_dir = "../../pipeline_output/09_UMI_group"
# output_dir = "../../pipeline_output/12_pandas_processed"


# START FUNCTIONS # ------------------------------------------------------------------------

# SORT INTO RNA OR DNA LIST
def split_to_DNA_RNA(lst):
    RNA_lst, DNA_lst = [], []
    for name in lst:

        # Check DNA
        if name[0] == 'd':
            # print("Found DNA")
            DNA_lst.append(name)

        # Check RNA
        elif name[0] == 'r':
            # print("Found RNA")
            RNA_lst.append(name)

        # Unknown type
        else:
            print(f"UNKNOWN TYPE FOUND: {name}")

    return DNA_lst, RNA_lst


def clean_pd_df(lst, name):
    counter = 1
    lst_clean = []
    for i in range(0, len(lst) - 1):
        df = lst[i]
        try:
            df["rsID_UMI"] = df["contig"] + ":" + df["final_umi"]
            df = df[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
            df = df.sort_values(by=["contig"])
            dfcontig = df["contig"].str.split("_", n=1, expand=True)
            df["rsID"] = dfcontig[0]
            df["variant"] = dfcontig[1]
            df_drop_dup = df.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
            df_drop_dup = df_drop_dup.rename({'final_umi_count': f'{name}rep{counter}'}, axis=1)
        except:
            df_drop_dup = None
        counter += 1
        lst_clean.append(df_drop_dup)
    return lst

# MERGE DNA RNA DFs
def merge_DNA_RNA(DNA_df, RNA_df):
    DNA_RNA_merged = []

    for i in range(0, len(DNA_df)):
        dna_name = DNA_df[i].columns[3]
        rna_name = RNA_df[i].columns[3]
        merged_df = pd.merge(left=DNA_df[i], right=RNA_df[i], how="left", left_on='rsID_UMI', right_on='rsID_UMI')
        merged_sort = merged_df.sort_values(["contig_x", "rsID_x"])
        DNA_RNA_merged.append(merged_sort[["rsID_UMI", f"{dna_name}", f"{rna_name}"]])

    return DNA_RNA_merged


def main():
    # Open the file
    filename_lst = [filename for filename in os.listdir(input_dir)]

    # Split list into DNA or RNA
    DNA_lst, RNA_lst = split_to_DNA_RNA(filename_lst)

    pd_DNA = [pd.read_csv(f"{input_dir}/{name}", delimiter='\t') for name in DNA_lst]
    pd_RNA = [pd.read_csv(f"{input_dir}/{name}", delimiter='\t') for name in RNA_lst]

    n = 10 - len(DNA_lst)
    pd_DNA += [None] * n
    pd_RNA += [None] * n

    print("OPENED ALL FILES IN PANDAS FOR pandas_processing.py in", round(time.time() - start_time, 2), " sec.")

    # Clean DNA dataframes
    DNA_sorted = clean_pd_df(pd_DNA, name="DNA")

    # Clean RNA dataframes
    RNA_sorted = clean_pd_df(pd_RNA, name="RNA")

    print("CLEANED ALL DNA/RNA DATAFRAMES FOR pandas_processing.py in", round(time.time() - start_time, 2), " sec.")

    # Merge the RNA and DNA files
    merged_DNA_RNA = merge_DNA_RNA(DNA_sorted, RNA_sorted)
    
    print("FINISHED MERGING ALL DATAFRAMES for pandas_processing.py in ", round(time.time() - start_time, 2), " sec.")

    # Merge to single file
    finaldf = pd.concat(DNA_RNA_merged, axis=1, join='inner')
    finaldf = finaldf.loc[:,~finaldf.columns.duplicated()].fillna(0)

    # Write to csv
    finaldf.to_csv(f"{output_dir}/DNA_RNA_count.csv", index=False)


###############
## EXECUTION ##
###############
if __name__ == '__main__':
    start_time = time.time()
    main()
    print(f'Finished pandas_processing.py in {round(time.time() - start_time, 2)} sec.')
