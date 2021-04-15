# STARR Seq to MPRAnalyze input


## Overview

STARR_Seq_to_MPRAnalyze_input is a pipeline that processes STARR Seq data (*.fastq.gz*) to an input file for MPRAnalyze (*.csv*). This repository contains the workflow and scripts for processing. Many of these scripts were adopted and curated from Arpit Misha (post-doc in the Hawkins lab). 

## Installation 
1. Confirm that [conda](https://docs.conda.io/en/latest/miniconda.html) is installed.
2. Clone this repository into the location you want to run the pipeline.
3. Create and activate the provided [environment](https://github.com/hawkins-lab/STARR_Seq_to_MPRAnalyze_input/tree/main/envs
): 

```
git clone https://github.com/hawkins-lab/STARR_Seq_to_MPRAnalyze_input.git \ 
&& cd STARR_Seq_to_MPRAnalyze_input/ \ 
&& conda env create -f STARRSeq2MPRAnalyze_env \ 
&& conda env create -f umitools_env
```

## Running the pipeline

We run all heavy computational processes on 


