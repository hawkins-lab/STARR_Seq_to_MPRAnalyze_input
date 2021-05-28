#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=2.0G
#$ -R y
#$ -l h_rt=1:0:0

output_dir=$1

mkdir -p ${output_dir}/02_trimmed_files
mkdir -p ${output_dir}/03_pandaseq_files
mkdir -p ${output_dir}/04_cutadapt_files/
mkdir -p ${output_dir}/05_umi_extract/
mkdir -p ${output_dir}/06_bwa_align/
mkdir -p ${output_dir}/07_sam2bam/
mkdir -p ${output_dir}/08_sam_sort/
mkdir -p ${output_dir}/09_UMI_group/
mkdir -p ${output_dir}/10_sam2txt/
mkdir -p ${output_dir}/11_text2tsv/
mkdir -p ${output_dir}/11_text2tsv/log
mkdir -p ${output_dir}/12_pandas_processed
mkdir -p ${output_dir}/13_pandas_mpranalyze
