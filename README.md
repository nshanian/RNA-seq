This workflow will perform STAR mapping / RSEM quantification

It is written to run as a batch script <sbatch> using SLURM in bash/shell on RNA-seq data in <.fastq> format. 

In addition to STAR and RSEM, this workflow will require bedGraphToBigWig, ucsc_tools, samtools, bedtools and r programs.
Reference genome files for STAR and RSEM are also required. 

https://github.com/alexdobin/STAR  ,   https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

http://deweylab.biostat.wisc.edu/rsem/

https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads

https://www.htslib.org/download/

https://bedtools.readthedocs.io/en/latest/

https://www.r-project.org/

usage: from an empty working directory, run ./STAR_RSEM.sh (read1) (read2 or "") (STARgenomeDir) (RSEMrefDir) (dataType) (nThreadsSTAR) (nThreadsRSEM) input: gzipped fastq file read1 [read2 for paired-end] /

STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script / 

read1=$1 #gzipped fastq file for read1 /

read2=$2 #gzipped fastq file for read1, use "" if single-end /

STARgenomeDir=/< > # Full path to STAR index /

RSEMrefDir=/< >  # Full path to RSEM index /

dataType=$3 # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE /

nThreadsSTAR=$4 # number of threads for STAR /

nThreadsRSEM=$5 # number of threads for RSEM /

output: all in the working directory

{filename}.sortedByCoord.bam                   alignments, standard sorted BAM, agreed upon formatting

{filename}.toTranscriptome.bam                 transcriptome alignments, used for expression quantification

Log.final.out                                  mapping statistics to be used for QC, text, STAR formatting

{filename}.genes.results                       RSEM gene quantifications, tab separated text, RSEM formatting

{filename}.isoforms.results                    RSEM transcript quantifications, tab separated text, RSEM formatting

{filename}.pdf                                 RSEM diagnostic plots

Signal.{Unique,UniqueMultiple}.strand{+,-}.bw  4 bigWig files for stranded data

Signal.{Unique,UniqueMultiple}.unstranded.bw   2 bigWig files for unstranded data


