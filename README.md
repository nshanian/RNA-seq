This repository contain tools for performing alignment and quantitation of transcript abundances from RNA-seq data in <.fastq> format.

They are written to run as a batch scripts (sbatch .sh) using SLURM in bash/shell. 

In addition to STAR and RSEM, this workflow will require trim_galore, bedGraphToBigWig, ucsc_tools, samtools, bedtools and r programs.

Reference genome files for STAR and RSEM are also required. 

This workflow uses STAR v2.5.3a and RSEM v1.2.30 and GRCh38/hg38 assembly of the human genome: 

https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/. 

https://github.com/ENCODE-DCC/rna-seq-pipeline

Please see below for information on each module used in the workflow:

https://github.com/alexdobin/STAR  ,   https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

http://deweylab.biostat.wisc.edu/rsem/

https://pachterlab.github.io/kallisto/about

https://github.com/FelixKrueger/TrimGalore

https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads

https://www.htslib.org/download/

https://bedtools.readthedocs.io/en/latest/

https://www.r-project.org/
