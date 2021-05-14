# STARR Seq to MPRAnalyze input


## Overview

STARR_Seq_to_MPRAnalyze_input is a pipeline that processes STARR Seq data (*.fastq.gz*) to an input file for MPRAnalyze (*.csv*). This repository contains the workflow and scripts for processing this data, and many of these scripts were adopted and curated from **Arpit Misha** (post-doc in the Hawkins lab). If there are any questions/bugs/errors, please contact me at <chrhsu@uw.edu> or <arpitm@uw.edu>.

## Installation 
1. Confirm that [conda](https://docs.conda.io/en/latest/miniconda.html) is installed.
2. Clone this repository into the location you want to run the pipeline.
3. Create and activate the provided [environment](https://github.com/hawkins-lab/STARR_Seq_to_MPRAnalyze_input/tree/main/envs
): 

```
git clone https://github.com/hawkins-lab/STARR_Seq_to_MPRAnalyze_input.git \ 
&& cd STARR_Seq_to_MPRAnalyze_input/env/ \ 
&& conda env create -f STARRSeq2MPRAnalyze_env \ 
&& conda env create -f umitools_env
```

## Running the pipeline
1. Navigate to the STARR_Seq_to_MPRAnalyze_input directory.
2. Add any Starr Seq files into the **pipeline_input** directory. These should be .fasta.gz files. 
3. Submit a job to the computation cluster with the **run_pipeline.sh** script.
4. Wait ~8 hours for the data to be processed
5. Check pipeline_output/13_final_mpranalyze_input/ directory for the MPRAnalyze input .csv files.

Just a few notes for reference:
- We run all heavy computational processes on University of Washington's Genome Sciences cluster computing.


