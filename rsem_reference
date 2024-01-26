#!/bin/bash 
#SBATCH --job-name=rsem
#SBATCH --time=2-00:00:00
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH -c 2
#SBATCH --account=user_id
#SBATCH --mem=40G

module load rsem/1.2.30

rsem-prepare-reference  --gtf /reference/RefGenomes/H_sapiens/hg38/illumina_igenomes_homo_sapiens_ucsc_hg38.gtf /reference/RefGenomes/GATK_Resource_Bundle/hg38/Homo_sapiens_assembly38.fasta hg38.rsem

