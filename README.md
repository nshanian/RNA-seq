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

