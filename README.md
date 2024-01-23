This workflow will perform STAR mapping / RSEM quantification on paired-end RNA-seq data in <.fastq> format.

It is written to run as a batch script (sbatch .sh) using SLURM in bash/shell. 

In addition to STAR and RSEM, this workflow will require trim_galore, bedGraphToBigWig, ucsc_tools, samtools, bedtools and r programs.

Reference genome files for STAR and RSEM are also required. 

This worfklow uses STAR v2.5.3a and RSEM v1.2.30 and GRCh38/hg38 assembly of the human genome: 

https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/. 

Please see below for infromation on each module used in the workflow:

https://github.com/alexdobin/STAR  ,   https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

http://deweylab.biostat.wisc.edu/rsem/

https://github.com/FelixKrueger/TrimGalore

https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads

https://www.htslib.org/download/

https://bedtools.readthedocs.io/en/latest/

https://www.r-project.org/
