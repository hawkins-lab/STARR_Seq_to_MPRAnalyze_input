#$ -S /bin/bash
#$ -pe serial 1
#$ -l mfree=8G -l h_rt=100:00:00



bwa index -p aid_bwaindex aid_snipa.fasta
