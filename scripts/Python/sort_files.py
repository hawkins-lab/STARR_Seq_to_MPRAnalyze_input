# File: sort_files.py
# AUTHOR: CHRIS HSU <chrhsu@uw.edu>
# CREATED DATE: APRIL 01, 2021

import os
import time
import shutil
import argparse

"""
This script is to sort all of the files into either DNA or RNA. Then, submit cluster job for pipeline.sh
-Input: data/fastq dump
-Output_DNA: pipeline_output/01_sorted_files/DNA
-Output_RNA: pipeline_output/01_sorted_files/RNA
"""

# Parser Commands # ------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Inputs for sorting fastq.gz input files.')
requiredNamed = parser.add_argument_group('Required named arguments')
requiredNamed.add_argument('-d', '--directory', type=str, help='Directory to pipeline_input', required=True)
requiredNamed.add_argument('-s', '--scripts', type=str, help='Directory to scripts', required=True)
requiredNamed.add_argument('-o', '--outputdir', type=str, help='Output Directory)', required=True)
args = parser.parse_args()

input_dir = args.directory
output_dir = args.outputdir
scripts_dir = args.scripts


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


# COPIES FILES AND MOVES THEM TO A NEW DIRECTORY
def copy_files_new_dir(lst, nt_type):
    # Check if directories exist, if not, make them
    if not os.path.exists(f'{output_dir}/01_sorted_files'):
        os.mkdir(f'{output_dir}/01_sorted_files')

    # Check if even number of files
    if len(lst) % 2 != 0:
        print("ERROR: Missing R1 or R2 files!")
        exit()

    # Make new directories
    for i in range(0, int(len(lst) / 2)):
        if not os.path.exists(f'{output_dir}/01_sorted_files/{nt_type}_{i + 1}'):
            os.mkdir(f'{output_dir}/01_sorted_files/{nt_type}_{i + 1}')

    # Copying Files to new directory
    counter = 1
    for i in range(0, len(lst) - 1):

        print("Counter: ", i)

        # Check current and next file names
        curr_name = lst[i].split('_')[0]
        next_name = lst[i + 1].split('_')[0]

        # Name for the string
        name_2 = str(lst[i].split('_')[0] + '_' + lst[i].split('_')[1])

        # Move file to new directory
        if curr_name == next_name:
            shutil.copyfile(f'{input_dir}/{lst[i]}', f'{output_dir}/01_sorted_files/{nt_type}_{counter}/{lst[i]}')
            shutil.copyfile(f'{input_dir}/{lst[i + 1]}',
                            f'{output_dir}/01_sorted_files/{nt_type}_{counter}/{lst[i + 1]}')

            # qsub another script
            print(f"Submitted {name_2}.")
            os.system(f"qsub {scripts_dir}/Shell/pipeline.sh {output_dir} {nt_type} {counter} {name_2}")

            counter += 1
            continue


def main():
    # Parse through the datafile names
    filename_lst = [filename for filename in os.listdir(input_dir)]

    # Split list into DNA or RNA
    DNA_lst, RNA_lst = split_to_DNA_RNA(filename_lst)

    # Need to sort here or downstream won't work on Linux
    DNA_lst_sorted = sorted(DNA_lst)
    RNA_lst_sorted = sorted(RNA_lst)

    # Copy to new directories
    copy_files_new_dir(DNA_lst_sorted, "DNA")
    copy_files_new_dir(RNA_lst_sorted, "RNA")


###############
## EXECUTION ##
###############

if __name__ == '__main__':
    start_time = time.time()
    main()
    print(f'Finished sort_files.py in {round(time.time() - start_time, 2)} sec.')

