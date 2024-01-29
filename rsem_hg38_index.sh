#!/bin/bash 
#SBATCH --job-name=rsem
#SBATCH --time=2-00:00:00
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH -c 2
#SBATCH --account=user_id # enter HPC user account information
#SBATCH --mem=40G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user_email # enter HPC user email

module load rsem/1.2.30

rsem-prepare-reference  --gtf /<path_to_genome_GTF_file>/H_sapiens/hg38/illumina_igenomes_homo_sapiens_ucsc_hg38.gtf /<path_to_genome_FastaFiles>/hg38/Homo_sapiens_assembly38.fasta hg38.rsem

