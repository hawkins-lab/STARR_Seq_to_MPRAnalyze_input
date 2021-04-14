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


def clean_pd_df(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, name):
    # X1
    try:
        x1["rsID_UMI"] = x1["contig"] + ":" + x1["final_umi"]
        x1 = x1[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x1 = x1.sort_values(by=["contig"])
        x1contig = x1["contig"].str.split("_", n=1, expand=True)
        x1["rsID"] = x1contig[0]
        x1["variant"] = x1contig[1]
        x1_drop_dup = x1.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x1_drop_dup = x1_drop_dup.rename({'final_umi_count': f'{name}rep1'}, axis=1)
    except:
        x1_drop_dup = None

    # X2
    try:
        x2["rsID_UMI"] = x2["contig"] + ":" + x2["final_umi"]
        x2 = x2[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x2 = x2.sort_values(by=["contig"])
        x2contig = x2["contig"].str.split("_", n=1, expand=True)
        x2["rsID"] = x2contig[0]
        x2["variant"] = x2contig[1]
        x2_drop_dup = x2.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x2_drop_dup = x2_drop_dup.rename({'final_umi_count': f'{name}rep2'}, axis=1)
    except:
        x2_drop_dup = None

    # X3
    try:
        x3["rsID_UMI"] = x3["contig"] + ":" + x3["final_umi"]
        x3 = x3[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x3 = x3.sort_values(by=["contig"])
        x3contig = x3["contig"].str.split("_", n=1, expand=True)
        x3["rsID"] = x3contig[0]
        x3["variant"] = x3contig[1]
        x3_drop_dup = x3.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x3_drop_dup = x3_drop_dup.rename({'final_umi_count': f'{name}rep3'}, axis=1)
    except:
        x3_drop_dup = None

    # X4
    try:
        x4["rsID_UMI"] = x4["contig"] + ":" + x4["final_umi"]
        x4 = x4[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x4 = x4.sort_values(by=["contig"])
        x4contig = x4["contig"].str.split("_", n=1, expand=True)
        x4["rsID"] = x4contig[0]
        x4["variant"] = x4contig[1]
        x4_drop_dup = x4.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x4_drop_dup = x4_drop_dup.rename({'final_umi_count': f'{name}rep4'}, axis=1)
    except:
        x4_drop_dup = None

    # X5
    try:
        x5["rsID_UMI"] = x5["contig"] + ":" + x5["final_umi"]
        x5 = x5[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x5 = x5.sort_values(by=["contig"])
        x5contig = x5["contig"].str.split("_", n=1, expand=True)
        x5["rsID"] = x5contig[0]
        x5["variant"] = x5contig[1]
        x5_drop_dup = x5.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x5_drop_dup = x5_drop_dup.rename({'final_umi_count': f'{name}rep5'}, axis=1)
    except:
        x5_drop_dup = None

    # X6
    try:
        x6["rsID_UMI"] = x6["contig"] + ":" + x6["final_umi"]
        x6 = x6[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x6 = x6.sort_values(by=["contig"])
        x6contig = x6["contig"].str.split("_", n=1, expand=True)
        x6["rsID"] = x6contig[0]
        x6["variant"] = x6contig[1]
        x6_drop_dup = x6.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x6_drop_dup = x6_drop_dup.rename({'final_umi_count': f'{name}rep6'}, axis=1)
    except:
        x6_drop_dup = None

    # X7
    try:
        x7["rsID_UMI"] = x7["contig"] + ":" + x7["final_umi"]
        x7 = x7[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x7 = x7.sort_values(by=["contig"])
        x7contig = x7["contig"].str.split("_", n=1, expand=True)
        x7["rsID"] = x7contig[0]
        x7["variant"] = x7contig[1]
        x7_drop_dup = x7.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x7_drop_dup = x7_drop_dup.rename({'final_umi_count': f'{name}rep7'}, axis=1)
    except:
        x7_drop_dup = None

    # X8
    try:
        x8["rsID_UMI"] = x8["contig"] + ":" + x8["final_umi"]
        x8 = x8[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x8 = x8.sort_values(by=["contig"])
        x8contig = x8["contig"].str.split("_", n=1, expand=True)
        x8["rsID"] = x8contig[0]
        x8["variant"] = x8contig[1]
        x8_drop_dup = x8.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x8_drop_dup = x8_drop_dup.rename({'final_umi_count': f'{name}rep8'}, axis=1)
    except:
        x8_drop_dup = None

    # X9
    try:
        x9["rsID_UMI"] = x9["contig"] + ":" + x9["final_umi"]
        x9 = x9[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x9 = x9.sort_values(by=["contig"])
        x9contig = x9["contig"].str.split("_", n=1, expand=True)
        x9["rsID"] = x9contig[0]
        x9["variant"] = x9contig[1]
        x9_drop_dup = x9.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x9_drop_dup = x9_drop_dup.rename({'final_umi_count': f'{name}rep9'}, axis=1)
    except:
        x9_drop_dup = None

    # X10
    try:
        x10["rsID_UMI"] = x10["contig"] + ":" + x9["final_umi"]
        x10 = x10[["contig", "rsID_UMI", "final_umi", "final_umi_count"]]
        x10 = x10.sort_values(by=["contig"])
        x10contig = x10["contig"].str.split("_", n=1, expand=True)
        x10["rsID"] = x10contig[0]
        x10["variant"] = x10contig[1]
        x10_drop_dup = x10.drop_duplicates(subset="rsID_UMI", keep="first")  # drop duplicates
        x10_drop_dup = x10_drop_dup.rename({'final_umi_count': f'{name}rep10'}, axis=1)
    except:
        x10_drop_dup = None

    return x1_drop_dup, x2_drop_dup, x3_drop_dup, x4_drop_dup, x5_drop_dup, \
           x6_drop_dup, x7_drop_dup, x8_drop_dup, x9_drop_dup, x10_drop_dup


def merged_all_df(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10,
                  r1, r2, r3, r4, r5, r6, r7, r8, r9, r10):
    # DF1
    try:
        merged_df1 = merged_pd_df(d1, r1, "DNArep1", "RNArep1")
    except:
        print(f"DNA_1 AND RNA_1 FILES DOES NOT EXIST!")
        merged_df1 = None

    # DF2
    try:
        merged_df2 = merged_pd_df(d2, r2, "DNArep2", "RNArep2")
    except:
        print(f"DNA_2 AND RNA_2 FILES DOES NOT EXIST!")
        merged_df2 = None

    # DF3
    try:
        merged_df3 = merged_pd_df(d3, r3, "DNArep3", "RNArep3")
    except:
        print(f"DNA_3 AND RNA_3 FILES DOES NOT EXIST!")
        merged_df3 = None

    # DF4
    try:
        merged_df4 = merged_pd_df(d4, r4, "DNArep4", "RNArep4")
    except:
        print(f"DNA_4 AND RNA_4 FILES DOES NOT EXIST!")
        merged_df4 = None

    # DF5
    try:
        merged_df5 = merged_pd_df(d5, r5, "DNArep5", "RNArep5")
    except:
        print(f"DNA_5 AND RNA_5 FILES DOES NOT EXIST!")
        merged_df5 = None

    # DF6
    try:
        merged_df6 = merged_pd_df(d6, r6, "DNArep6", "RNArep6")
    except:
        print(f"DNA_6 AND RNA_6 FILES DOES NOT EXIST!")
        merged_df6 = None

    # DF7
    try:
        merged_df7 = merged_pd_df(d7, r7, "DNArep7", "RNArep7")
    except:
        print(f"DNA_7 AND RNA_7 FILES DOES NOT EXIST!")
        merged_df7 = None

    # DF8
    try:
        merged_df8 = merged_pd_df(d8, r8, "DNArep8", "RNArep8")
    except:
        print(f"DNA_8 AND RNA_8 FILES DOES NOT EXIST!")
        merged_df8 = None

    # DF9
    try:
        merged_df9 = merged_pd_df(d9, r9, "DNArep9", "RNArep9")
    except:
        print(f"DNA_9 AND RNA_9 FILES DOES NOT EXIST!")
        merged_df9 = None

    # DF10
    try:
        merged_df10 = merged_pd_df(d10, r10, "DNArep10", "RNArep10")
    except:
        print(f"DNA_10 AND RNA_10 FILES DOES NOT EXIST!")
        merged_df10 = None

    return merged_df1, merged_df2, merged_df3, merged_df4, merged_df5, \
           merged_df6, merged_df7, merged_df8, merged_df9, merged_df10


def merged_pd_df(DNA_df, RNA_df, dna_name, rna_name):
    merged_df = pd.merge(left=DNA_df, right=RNA_df, how="left", left_on='rsID_UMI', right_on='rsID_UMI')

    merged_sort = merged_df.sort_values(["contig_x", "rsID_x"])
    return merged_sort[["rsID_UMI", f"{dna_name}", f"{rna_name}"]]


# MERGE THE DNA+RNA MERGED FILES
def merge_merged_pd_df(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10):

    overall_merged = m1
    try:
        # M1 and 2
        if m2 is not None:
            overall_merged = pd.merge(left=m1, right=m2, how="outer", left_on='rsID_UMI', right_on='rsID_UMI')
        else:
            print("NO 2nd FILE TO MERGE.")
            return m1.fillna(0)

        # M3
        if m3 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m3, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 3rd FILE TO MERGE.")
            return overall_merged.fillna(0)

        # M4
        if m4 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m4, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 4th FILE TO MERGE.")
            return overall_merged.fillna(0)

        # M5
        if m5 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m5, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 5th FILE TO MERGE.")
            return overall_merged.fillna(0)

        # M6
        if m6 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m6, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 6th FILE TO MERGE.")
            return overall_merged.fillna(0)

        # M7
        if m7 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m7, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 7th FILE TO MERGE.")
            return overall_merged.fillna(0)

        # M8
        if m8 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m8, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 8th FILE TO MERGE.")
            return overall_merged.fillna(0)

        # M9
        if m9 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m9, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 9th FILE TO MERGE.")
            return overall_merged.fillna(0)

        # M10
        if m10 is not None:
            overall_merged = pd.merge(left=overall_merged, right=m10, how="outer", left_on='rsID_UMI',
                                      right_on='rsID_UMI')
        else:
            print("NO 10th FILE TO MERGE.")
            return overall_merged.fillna(0)

        return overall_merged.fillna(0)

    except:
        print("AN ERROR HAS OCCURRED WHILE MERGING THE MERGED FILES!")

    pass


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

    # Each RNA/DNA file as pandas dataframe
    d1_sort, d2_sort, d3_sort, d4_sort, d5_sort, d6_sort, d7_sort, d8_sort, d9_sort, d10_sort = pd_DNA
    r1_sort, r2_sort, r3_sort, r4_sort, r5_sort, r6_sort, r7_sort, r8_sort, r9_sort, r10_sort = pd_RNA

    print("OPENED ALL FILES IN PANDAS FOR pandas_processing.py in", round(time.time() - start_time, 2), " sec.")

    # Clean DNA dataframes
    d1_sorted, d2_sorted, d3_sorted, d4_sorted, d5_sorted, d6_sorted, d7_sorted, d8_sorted, d9_sorted, d10_sorted = \
        clean_pd_df(d1_sort, d2_sort, d3_sort, d4_sort, d5_sort,
                    d6_sort, d7_sort, d8_sort, d9_sort, d10_sort, name="DNA")

    # Clean RNA dataframes
    r1_sorted, r2_sorted, r3_sorted, r4_sorted, r5_sorted, r6_sorted, r7_sorted, r8_sorted, r9_sorted, r10_sorted = \
        clean_pd_df(r1_sort, r2_sort, r3_sort, r4_sort, r5_sort,
                    r6_sort, r7_sort, r8_sort, r9_sort, r10_sort, name="RNA")

    print("CLEANED ALL DNA/RNA DATAFRAMES FOR pandas_processing.py in", round(time.time() - start_time, 2), " sec.")

    # Merge the RNA and DNA files
    merged1, merged2, merged3, merged4, merged5, merged6, merged7, merged8, merged9, merged10 = \
        merged_all_df(d1_sorted, d2_sorted, d3_sorted, d4_sorted, d5_sorted,
                      d6_sorted, d7_sorted, d8_sorted, d9_sorted, d10_sorted,
                      r1_sorted, r2_sorted, r3_sorted, r4_sorted, r5_sorted,
                      r6_sorted, r7_sorted, r8_sorted, r9_sorted, r10_sorted)

    print("FINISHED MERGING ALL DATAFRAMES for pandas_processing.py in ", round(time.time() - start_time, 2), " sec.")

    # Merge the merged DNA/RNA files
    final_merged_df = merge_merged_pd_df(merged1, merged2, merged3, merged4, merged5, merged6, merged7, merged8, merged9, merged10)

    # Write to csv
    final_merged_df.to_csv(f"{output_dir}/DNA_RNA_count.csv", index=False)


###############
## EXECUTION ##
###############
if __name__ == '__main__':
    start_time = time.time()
    main()
    print(f'Finished pandas_processing.py in {round(time.time() - start_time, 2)} sec.')
