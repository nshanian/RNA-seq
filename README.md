## RNA-seq data analysis 

## Part 1: Alignment and Quantitation

This repository contains workflows and tools for performing alignment and quantitation of transcript abundances from RNA-seq data. 

The workflows can be run on paired-end (PE) or single-end (SE) next-generation sequcning (NGS) of bulk, as well as single-cell RNA-seq data in `.fastq` format.

The workflows are written to run as batch shell scripts using `SLURM` on a high performance computing (HPC) cluster. 

The main `STAR_RSEM.sh` workflow uses STAR v2.5.3a and `rsem` v1.2.30 and GRCh38/hg38 assembly of the human genome:
  
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/.

Reference genome files in `<.gtf>` and `<.fasta>` formats for `STAR` and `rsem` are required for generating the reference indexes. 

Appropriate reference genome files can be downloaded from GENCODE's website:

https://www.gencodegenes.org/human/release_38.html

Examples for generating `STAR` and `rsem` indexes from references genome files are provided for the hg38 assembly of the human genome.

For documentation on an additional `STAR_RSEM.ENCODE.sh` workflow from the ENCODE consortium see:

https://github.com/ENCODE-DCC/rna-seq-pipeline

In addition to `STAR` and `rsem`, the modules in the workflows will require: 

`trim_galore`, `bedGraphToBigWig`, `ucsc_tools`, `picard-tools`, `samtools`, `bedtools`, `bowtie`, `tophat`, `fastqc`, `igvtools`, `cufflinks` and `r` programs.

Please see below for information on each module used in the provided workflows.

## Documentation and References:

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
