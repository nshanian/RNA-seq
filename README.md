The files in this repository contain workflows and tools for performing alignment and quantitation of transcript abundances from RNA-seq data.

They are written to run as batch scripts (sbatch .sh) using SLURM in bash/shell. 

In addition to STAR and RSEM, STAR_RSEM.sh workflow will require trim_galore, bedGraphToBigWig, ucsc_tools, samtools, bedtools and r programs.

Reference genome files for STAR and RSEM are also required. 

The STAR_RSEM.sh workflow uses STAR v2.5.3a and RSEM v1.2.30 and GRCh38/hg38 assembly of the human genome: 

https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/. 

https://github.com/ENCODE-DCC/rna-seq-pipeline

Please see below for information on each module used in the provided workflows:

https://github.com/alexdobin/STAR  ,   https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

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
