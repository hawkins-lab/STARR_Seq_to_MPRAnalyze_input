#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=5.0G
#$ -R y
#$ -l h_rt=1:0:0

#Retrieve my CPU conda environment from my bashrc
source /net/hawkins/vol1/home/chrhsu/.bashrc
cd /net/hawkins/vol1/home/chrhsu/
source activate STARRSeq2MPRAnalyze_snakemake


# configure file paths
SNAKE_FILE='../workflow/Snakefile'


# Execution
snakemake --snakefile $SNAKE_FILE

