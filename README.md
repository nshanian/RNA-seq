## RNA-seq data analysis 

## Part 1: Alignment and Quantitation

This repository contains workflows and tools for performing alignment and quantitation of transcript abundances from RNA-seq data.

The workflows are written to run as batch shell scripts using `SLURM` on a high performance computing (HPC) cluster. 

In addition to `STAR` and `RSEM`, the modules in the workflows will require: 

`trim_galore`, `bedGraphToBigWig`, `ucsc_tools`, `samtools`, `bedtools` and `r` programs.

Reference genome files `<.gtf>` and `<.fasta>` for `STAR` and `RSEM` are also required for generating the reference indexes. 

Appropriate reference genome files can be downloaded from GENCODE's website:

https://www.gencodegenes.org/human/release_38.html

The `STAR_RSEM.sh` workflow uses STAR v2.5.3a and RSEM v1.2.30 and GRCh38/hg38 assembly of the human genome:
  
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/.

https://github.com/ENCODE-DCC/rna-seq-pipeline

## Documentation and References:

Please see below for information on each module used in the provided workflows.
  
https://github.com/alexdobin/STAR , https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

http://deweylab.biostat.wisc.edu/rsem/
  
https://pachterlab.github.io/kallisto/about

https://github.com/FelixKrueger/TrimGalore

https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads

https://www.htslib.org/download/
  
https://bedtools.readthedocs.io/en/latest/
  
https://broadinstitute.github.io/picard/
  
https://bowtie-bio.sourceforge.net/index.shtml

https://ccb.jhu.edu/software/tophat/index.shtml

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  
https://cole-trapnell-lab.github.io/cufflinks/
  
https://github.com/igvteam/igv

https://www.r-project.org/
  
