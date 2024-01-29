#!/bin/bash 
#SBATCH --job-name=STAR
#SBATCH --time=2-00:00:00
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH -c 2
#SBATCH --account=user_id # enter user account information
#SBATCH --mem=40G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user_email # enter user email

module load STAR

STAR --runThreadN 20
--runMode genomeGenerate 
--genomeDir /<path_to_genome_directory>/RNA-seq/NCBI_refseq/ 
--genomeFastaFiles /<path_to_genome_FastaFiles>/hg38/Homo_sapiens_assembly38.fasta 
--sjdbGTFfile /<path_to_genome_GTF_file>/H_sapiens/hg38/illumina_igenomes_homo_sapiens_ucsc_hg38.gtf 
--sjdbOverhang 100 
--outFileNamePrefix star_index /

