#$ -S "/bin/bash"
#$ -cwd
#$ -o "job.stdout"
#$ -e "job.stderr"
#$ -l m_mem_free=5.0G
#$ -R y
#$ -P sage
#$ -l h_rt=32:0:0

output_dir=$1
nt_type=$2
counter=$3
filename=$4

# Go to this working directory, change as need
cd /net/hawkins/vol1/home/chrhsu/proj/STARR_Seq_to_MPRAnalyze
conda activate snakemake_STARRSeq2MPRAnalyze

mkdir ${output_dir}/02_trimmed_files/${nt_type}_${counter}
mkdir ${output_dir}/03_pandaseq_files/${nt_type}_${counter}
mkdir ${output_dir}/04_cutadapt_files/${nt_type}_${counter}
mkdir ${output_dir}/05_umi_extract/${nt_type}_${counter}


# Takes long time (30 minutes)
trim_galore --no_report_file --cores 4 --gzip -o ${output_dir}/02_trimmed_files/${nt_type}_${counter} --paired ${output_dir}/01_sorted_files/${nt_type}_${counter}/${filename}_R1.fastq.gz ${output_dir}/01_sorted_files/${nt_type}_${counter}/${filename}_R2.fastq.gz

# Takes long time (25 minutes)
pandaseq -F -f ${output_dir}/02_trimmed_files/${nt_type}_${counter}/${filename}_R1_val_1.fq.gz -r ${output_dir}/02_trimmed_files/${nt_type}_${counter}/${filename}_R2_val_2.fq.gz -d rbfkms -u ${output_dir}/03_pandaseq_files/${nt_type}_${counter}/unmerged_pandaseq${file}.fa 2> ${output_dir}/03_pandaseq_files/${nt_type}_${counter}/${filename}pandastat.txt 1> ${output_dir}/03_pandaseq_files/${nt_type}_${counter}/merged_pandaseq_${filename}.fastq

# cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq
# Fast
cutadapt  -j 0 -u 50 -u -40 -o ${output_dir}/04_cutadapt_files/${nt_type}_${counter}/trimmed_merged_pandaseq_${filename}.fastq ${output_dir}/03_pandaseq_files/${nt_type}_${counter}/merged_pandaseq_${filename}.fastq

# UMITOOLS Conda env
conda deactivate
conda activate umitools_latest

umi_tools extract --extract-method=regex --stdin=${output_dir}/04_cutadapt_files/${nt_type}_${counter}/trimmed_merged_pandaseq_${filename}.fastq --bc-pattern='.+(?P<umi_1>AATGAT(.{5}))' --log=${output_dir}/05_umi_extract/${nt_type}_${counter}/${filename}_UMI.log --stdout ${output_dir}/05_umi_extract/${nt_type}_${counter}/${filename}_UMI.fastq.gz


# snakemake_STARRSeq2MPRAnalyze
conda deactivate
conda activate snakemake_STARRSeq2MPRAnalyze
bwa mem -t 4 bwa_index_aid/aid_bwaindex ${output_dir}/05_umi_extract/${nt_type}_${counter}/${filename}_UMI.fastq.gz > ${output_dir}/06_bwa_align/${filename}_UMI_bwa.err > ${output_dir}/06_bwa_align/${filename}_UMI.sam

##################
#### SAMTOOLS ####
##################
samtools view -bS  ${output_dir}/06_bwa_align/${filename}_UMI.sam >  ${output_dir}/07_sam2bam/${filename}_UMI.bam
samtools sort ${output_dir}/07_sam2bam/${filename}_UMI.bam -o ${output_dir}/08_sam_sort/${filename}_UMI_sorted.bam
samtools index ${output_dir}/08_sam_sort/${filename}_UMI_sorted.bam
samtools view ${output_dir}/06_bwa_align/${filename}_UMI.sam | cut -f 1,3 > ${output_dir}/10_sam2txt/${filename}_UMIsam.txt

###################
#### UMI_TOOLS ####
###################
conda deactivate
conda activate umitools_latest

umi_tools group -I ${output_dir}/08_sam_sort/${filename}_UMI_sorted.bam  --group-out=${output_dir}/09_UMI_group/${filename}_UMI_sorted.tsv
umi_tools count_tab -I ${output_dir}/10_sam2txt/${filename}_UMIsam.txt -S ${output_dir}/11_text2tsv/${filename}_UMIsam.tsv -L ${output_dir}/11_text2tsv/log/${filename}_UMIsamcount.lo


