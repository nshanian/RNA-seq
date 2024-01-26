#!/bin/bash -l
# Set the name of the job
#SBATCH --job-name=kallisto2bam
#
# Set the maximum memory allowed
#SBATCH --mem-per-cpu=12G
#
# Set the maximum run time
#SBATCH -t 48:00:00
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user_email
#
# The number of threads we will require
#
# Set output and error log files
#SBATCH -o /<path to working directory>/log/out.txt
#SBATCH -e /<path to working directory>log/error.txt
#
# set the account for billing (user_id)
#SBATCH --account=user_id
#
#
#SBATCH --export=ALL

############## kallisto is a program for quantifying abundances of transcripts from bulk and single-cell RNA-seq data, or more generally of target sequences using high-throughput sequencing reads. 
############## It is based on the novel idea of pseudoalignment for rapidly determining the compatibility of reads with targets, without the need for alignment.
###### For documentation see https://pachterlab.github.io/kallisto/about
###### To submit this script, type this on the command line in shell from the directory where this script is:    CHANGE variables here before submitting
#
#  for sampleID in cnt1 cnt2 cnt3 rep1 rep2 rep3
#  do
#      sbatch -J kallisto_${sampleID} --export=sample=${sampleID},workingDir=/<path to working directory>,fastqDir=/<path to fastq files>/,barcodes=/<path to barcodes.txt file> kallisto.sh
#  done
#######

module purge
module load anaconda
source activate kallisto_0.45.0


# the human kallisto index file (e.g. Homo sapiens hg38)
index=/<path to directory>/Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx

# output directory
outputdir=${workingDir}/${sample}

cd ${outputdir}

# collect fastq files
sampleBarcode=`awk 'BEGIN { FS = "\t" } ; ($2 == "'${sample}'") {print $1}' ${barcodes}`

fastqOneFilePairA=`ls ${fastqDir}/*${sampleBarcode}*_1_pf.fastq.gz | head -1`
fastqTwoFilePairA=`ls ${fastqDir}/*${sampleBarcode}*_2_pf.fastq.gz | head -1`
#fastqOneFilePairB=`ls ${fastqDir}/*${sampleBarcode}*_1_pf.fastq.gz | tail -1`
#fastqTwoFilePairB=`ls ${fastqDir}/*${sampleBarcode}*_2_pf.fastq.gz | tail -1`




# run kallisto for normal output
    kallisto quant -i ${index} -o ${outputdir} -b 100 ${fastqOneFilePairA} ${fastqTwoFilePairA}

# run kallisto again to get bam files
    kallisto quant -i ${index} -o ${outputdir} -b 30 --genomebam --gtf /<path_to_gtf_file>/Homo_sapiens.GRCh38.94.gtf.gz  --chromosomes /<path_to_chr_sizes_file>/hg38_chr_sizes.txt  ${fastqOneFilePairA} ${fastqTwoFilePairA} #| samtools view -Sb -> ${sample}.genome.bam

# sort and index bam files (kallisto does sorting and indexing automatically)
#    samtools sort -m 12G ${sample}.genome.bam ${sample}.genome.sorted
#    samtools index ${sample}.genome.sorted.bam

module purge
module add samtools/1.2
    samtools flagstat pseudoalignments.bam > ${sample}.pseudoalignments.bam_flagstats


