#$ -S /bin/bash
#$ -pe serial 1
#$ -l h_rt=100:00:00 -l mfree=5G
# "/net/hawkins/vol1/arpit/STARR_seq_fastafiles/STARRseqdata/Th1_Failed_PairedendRUN_18June_2020/catfiles"

### trimgalore for trimming paired end reads

trim_galore --no_report_file --gzip --paired ${file}_R1.fastq ${file}_R2.fastq

### pandaseq to make single read from trimmed paired end read

pandaseq -F -f ${file}_R1_val_1.fq.gz  -r ${file}_R2_val_2.fq.gz  -d rbfkms -u unmerged_pandaseq${file}.fa 2> ${file}pandastat.txt 1>merged_pandaseq_${file}.fastq


### removes the 5 prime and 3 prime extra region (do MSA before this to confirm how much to chop from right or left

cutadapt -u 50 -u -40 -o trimmed_merged_pandaseq_d1Th1_S4.fastq merged_pandaseq_d1Th1_S4.fastq

###works with cutadapt -u 50 -u -40 but adds AATGAT to barcode

conda activate umitools_latest

umi_tools extract  --extract-method=regex --stdin=trimmed_merged_pandaseq_d1Th1_S4.fastq --bc-pattern='.+(?P<umi_1>AATGAT(.{5}))' --log=d1Th1S4UMI.log --stdout d1Th1_S4_UMI.fastq.gz

conda deactivate
### next steps are allignment using BWA and SAM BAM generation
bwa mem -t 4 bwa_index_aid/aid_bwaindex  r1Th1_S1_UMI.fastq.gz > r1Th1_S1_UMI_bwa.err > r1Th1_S1_UMI.sam

samtools view -bS  r1Th1_S1_UMI.sam > r1Th1_S1_UMI.bam

samtools sort r1Th1_S1_UMI.bam -o r1Th1_S1_UMI_sorted.bam

samtools index r1Th1_S1_UMI_sorted.bam

### generates UMIsam TSV for generating MPRAnalyze input in pandas

conda activate umitools_latest

umi_tools group -I r1Th1_S1_UMI_sorted.bam  --group-out=r1Th1_S1_UMI_sorted.tsv  --output-bam --log =r1Th1_S1_UMI_sorted.log

samtools view r1Th1_S1_UMI.sam | cut -f 1,3 > r1Th1_S1_UMIsam.txt

umi_tools count_tab -I r1Th1_S1_UMIsam.txt -S r1Th1_S1_UMIsam.tsv -L r1Th1_S1_UMIsamcount.lo