#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=10.0G
#$ -R y
#$ -l h_rt=24:0:0

# Go to this working directory, change as need
cd /net/hawkins/vol1/home/{$USER}/STARR_Seq_to_MPRAnalyze

# Make all sub-directories
bash scripts/Shell/make_dir.sh pipeline_output

# Sort files with Python script and start first part of data processing pipeline
python3 scripts/Python/sort_files.py -d pipeline_input -o pipeline_output -s scripts

# Sleep
sleep 8h

# Continue with jupyter notebook portion of the data processing
python3 scripts/Python/pandas_processing.py -d pipeline_output/09_UMI_group -o pipeline_output/12_pandas_processed
python3 scripts/Python/pandas_format_mpranalyze.py -d pipeline_output/12_pandas_processed -o pipeline_output/13_final_mpranalyze_input


##### Other R scripts you want to put here for MPRAnalyze
# conda activate r_env
#R CMD BATCH --no-save --no-restore scripts/R/mpra_barcodeallelic.R

