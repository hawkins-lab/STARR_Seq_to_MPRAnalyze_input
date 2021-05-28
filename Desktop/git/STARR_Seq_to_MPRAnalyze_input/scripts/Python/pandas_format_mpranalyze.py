# File: pandas_format_mpranalyze.py
# AUTHOR: CHRIS HSU <chrhsu@uw.edu>
# CREATED DATE: APRIL 13, 2021

import os
import time
import argparse
import pandas as pd
import numpy as np

"""
This script is to process the files output from the pipeline.sh script.
-Input: 12_UMI_group directory
-Output: .csv file containing DNA and RNA counts

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

# input_dir = "../../pipeline_output/12_pandas_processed"
# output_dir = "../../pipeline_output/13_final_mpranalyze"


# START FUNCTIONS # ------------------------------------------------------------------------

def filter_values_DNA_RNA(df):
    for col in df.columns:
        # Remove values in DNA less than 5
        if 'DNA' in col:
            df = df[(df[col] >= 5)]

        # Remove all 0 values in RNA
        elif 'RNA' in col:
            df = df[(df[col] >= 1)]
    return df


def add_rep_df(df, val, idx):
    df.insert(idx, 'rep', f'rep{val}')
    df.columns = ["rsID_UMI", "count", "batch"]
    return df


def manipulate_pd_dfs(df, nt_type):
    # of samples
    num_columns = int((len(df.columns) - 1) / 2 + 1)

    dna_lst = [f"DNArep{i}" for i in range(1,num_columns)]
    dna_lst.insert(0,"rsID_UMI")
    rna_lst = [f"RNArep{i}" for i in range(1,num_columns)]
    rna_lst.insert(0,"rsID_UMI")
    
    # Split dataframe into RNA or DNA
    if nt_type == 'DNA':
        df2 = df[dna_lst[0:num_columns]]
    elif nt_type == 'RNA':
        df2 = df[rna_lst[0:num_columns]]
    
    # Each replicate has own dataframe
    count_all_lst = []
    for idx in range(1, num_columns):
        rep = df2[["rsID_UMI", f"{nt_type}rep{idx}"]]
        rep = add_rep_df(rep, idx, 2)
        count_all_lst.append(rep)
       

    count_all_df = pd.concat(count_all_lst).reset_index(drop=True)
    count_all_df[["rsID", "barcode"]] = count_all_df.rsID_UMI.str.split(":", expand=True)
    count_all_df[["SNP", "allele"]] = count_all_df.rsID.str.split("_", expand=True)
    count_all_df["barcode_rep"] = count_all_df["rsID_UMI"] + "," + count_all_df["batch"]
    count_all_df.columns = ["rsID_UMI", "count", "batch", "rsID", "UMI", "SNP", "allele", "barcode"]

    final_count_df = count_all_df[["SNP", "barcode", "count"]]


    ## final dna counts
    fdc_rs = final_count_df[final_count_df["SNP"].str.match(r'(^rs*.*)') == True]

    ## fdc to pivot table to make matrix
    fdcrs_mat = fdc_rs.pivot_table(columns='barcode', index='SNP', values='count').reset_index()
    fdcrs_mat = fdcrs_mat.fillna(0)

    fdc_rs_1 = fdc_rs['barcode'].str.extract("(?P<rsID>.*?)_(?P<allele>.*?):(?P<UMI>.*?),(?P<batch>.*)", expand=True)
    fdc_rs_a = fdc_rs.join(fdc_rs_1)

    return fdc_rs_a, fdcrs_mat


def main():
    # Open file
    input_file = f"{input_dir}/DNA_RNA_count.csv"
    df = pd.read_csv(input_file, index_col=False)

    # Drop unwanted value rows
    df_filter = filter_values_DNA_RNA(df)

    final_dna_count, final_dna_matrix = manipulate_pd_dfs(df_filter, "DNA")
    final_rna_count, final_rna_matrix = manipulate_pd_dfs(df_filter, 'RNA')

    # Get Column Annotations
    colAnnotation_rs = final_dna_count[["barcode", "batch", "UMI", "allele"]]
    colAnnotation_rs.columns = ["observation", "batch", "barcode", "allele"]
    colAnnotation_rs["observation"] = colAnnotation_rs['observation'].str.replace(',', '.')

    # Write annotations to csv
    colAnnotation_rs.to_csv(f"{output_dir}/colAnnot_rs_dnagrt5_rnagrt1.csv", index=False, header=True)

    # Write DNA to csv
    final_dna_matrix.to_csv(f"{output_dir}/DNACounts_rs_dnagrt5_rnagrt1.csv", index=False, header=True)
    final_dna_count.to_csv(f"{output_dir}/dnacount_dnagrt5.csv", index=False, header=True)

    # Write RNA to csv
    final_rna_matrix.to_csv(f"{output_dir}/RNACounts_rs_dnagrt5_rnagrt1.csv", index=False, header=True)
    final_rna_count.to_csv(f"{output_dir}/rnacount_rnagrt1.csv", index=False, header=True)


###############
## EXECUTION ##
###############

if __name__ == '__main__':
    start_time = time.time()
    main()
    print(f'Finished pandas_format_mpranalyze.py in {round(time.time() - start_time, 2)} sec.')




