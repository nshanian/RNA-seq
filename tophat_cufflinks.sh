#!/bin/sh
#$ -o /<path_to_directory>/log/out.txt         
#$ -e /<path_to_directory>/log/error.txt
# maximum memory
#$ -l h_vmem=6G
# e-mail user when finished or aborted
#$ -m ae
#$ -M user_email
# check for errors
#$ -w e
# time allowed
#$ -l h_rt=24:00:00
#shared memory, number of threads
#$ -pe shm 6


############## TopHat is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie
############## It then analyzes the mapping results to identify splice junctions between exons.
############## Cufflinks assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-seq samples, expressed as Fragments Per Kilobase of transcript per Million (FPKM) mapped reads.
######For ducumentation see https://ccb.jhu.edu/software/tophat/index.shtml
######For ducumentation see https://cole-trapnell-lab.github.io/cufflinks/
######To submit this script, type this on the shell in the directory where this shell script is:    CHANGE variables here before submitting
#
#  for sampleID in cnt1 cnt2 cnt3 rep1 rep2 rep3
#  do
#      qsub -N cufflinks_${sampleID} -v sample=${sampleID},workingDir=/<path to working directory>,fastqDir=/<path_to_fastq_files>,barcodes=/<path_to_barcodes>/barcodes.txt tophat_plus_cuffFPKMs.sh
#  done
##############


module add picard-tools/1.92
module add samtools/0.1.19
module add bedtools/2.18.0
module add bowtie/2.2.1
module add tophat/2.0.11
module add fastqc/0.10.1
module add r/3.1.0
module add igvtools/2.3.3
module add cufflinks/2.2.1


# the human transcriptome file
geneGTF=/<path_to_reference_genome>/gencode.v19.annotation.gtf


cd ${workingDir}/${sample}

if [ ! -d RNAseq ]
then
    mkdir RNAseq
    mkdir RNAseq/tophatAlignment
    mkdir RNAseq/cufflinks
fi
cd RNAseq


sampleBarcode=`awk 'BEGIN { FS = "\t" } ; ($2 == "'${sample}'") {print $1}' ${barcodes}`

fastqOneFile=`ls ${fastqDir}/*${sampleBarcode}*_1_pf.fastq.gz`
fastqTwoFile=`ls ${fastqDir}/*${sampleBarcode}*_2_pf.fastq.gz`



### TopHat alignment of reads to the genome and discovery of transcript splice sites.
if [ ! -f tophatAlignment/accepted_hits.bam ]
then
    tophat -p 6 -G ${geneGTF} -o tophatAlignment /<path_to_bowtie_reference_file> ${fastqOneFile} ${fastqTwoFile}

    cd tophatAlignment
    samtools sort -m 12G accepted_hits.bam ${sample}.accepted_hits.sorted
    samtools index ${sample}.accepted_hits.sorted.bam
    samtools flagstat ${sample}.accepted_hits.sorted.bam > ${sample}.accepted_hits.sorted.bam_flagstats

fi


cd ${workingDir}/${sample}/RNAseq

# Cufflinks uses this map against the genome to assemble the reads into transcripts.

cufflinks -p 6 -o cufflinks -G ${geneGTF} tophatAlignment/${sample}.accepted_hits.sorted.bam


