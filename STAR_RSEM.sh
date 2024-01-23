#!/bin/bash 
#SBATCH --job-name=STAR
#SBATCH --time=2-00:00:00
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH -c 2
# set the account for billing (if running on a cluster)
#SBATCH --account=user id 
#SBATCH --mem=60G


# STAR mapping / RSEM quantification pipeline
# usage: from an empty working directory, run
# ./STAR_RSEM.sh (read1) (read2 or "") (STARgenomeDir) (RSEMrefDir) (dataType) (nThreadsSTAR) (nThreadsRSEM)
# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read1, use "" if single-end
STARgenomeDir=/< > # Full path to STAR index
RSEMrefDir=/< > # Full path to RSEM index
dataType=$3 # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsSTAR=$4 # number of threads for STAR
nThreadsRSEM=$5 # number of threads for RSEM

# output: all in the working directory
# ${filename}.sortedByCoord.bam                 # alignments, standard sorted BAM, agreed upon formatting
# ${filename}.toTranscriptome.bam               # transcriptome alignments, used for expression quantification
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# ${filename}.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# ${filename}.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# ${filename}.pdf                                     # RSEM diagnostic plots
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

##### executables
STAR=STAR                             
RSEM=rsem-calculate-expression        
bedGraphToBigWig=bedGraphToBigWig                          

module load trim_galore
##### Filenames and output directory
filename=$(basename "$1")
filename=`echo $filename | sed 's/.gz$//g' | sed 's/.fastq$//g'  | sed 's/.fq$//g' | sed 's/_R1//g' | sed 's/_1$//g' | sed 's/_r1$//g'` #"${filename%.*}"
indir=$(dirname "$1")
outdir=$indir/$filename\_star/
mkdir -p $outdir/logs

#### Reads trimming
case "$dataType" in
unstr_PE|str_PE)
  echo "Running QC filtering on PE reads ..."
  trim_galore --fastqc --paired --length 15 -o $outdir --dont_gzip $1 $2 2> $outdir/qualityfilter.log
  read1=`find $outdir/ -name "*_val_1.fq"`
  read2=`find $outdir/ -name "*_val_2.fq"`
  echo "Reads1 = $read1"
  echo "Reads2 = $read2"
  fastqc $read1 $read2
  ;;
unstr_SE|str_SE)
  echo "Running QC filtering on SE reads ..."
  trim_galore --fastqc --length 15 -o $outdir --dont_gzip $1 2> $outdir/qualityfilter.log
  read1=`find $outdir/ -name "*.fq"`
  read2=""
  echo "Reads1 = $read1"
  fastqc $read1
  ;;
esac

##### Unload/purge trim_galore and load the tools for the next steps
module purge
module load STAR/2.5.3a rsem/1.2.30 ucsc_tools samtools bedtools r

##### STAR
echo "Start mapping..."

# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn $read1 $read2   --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes All    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999  \
 --outFilterMismatchNoverReadLmax 0.07   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --outFileNamePrefix $outdir/"

#STARparCommon="--genomeDir $STARgenomeDir  --readFilesIn $read1 $read2 --outFileNamePrefix $outdir/ --outStd BAM_Unsorted   --outReadsUnmapped Fastx   --outSAMattributes All      --outSAMunmapped Within   --outSAMattrRGline ID:foo      --outFilterType BySJout   --outFilterMultimapNmax 10   --outFilterMultimapScoreRange 1   --outFilterScoreMin 10   --alignEndsType EndToEnd"
STARparCircRNA="--chimSegmentMin 12 --chimJunctionOverhangMin 12"

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN $nThreadsSTAR --limitBAMsortRAM 50000000000"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"


# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 
case "$dataType" in
str_SE|str_PE)
      #OPTION: stranded data
      STARparStrand=""
      STARparWig="--outWigStrand Stranded"
      ;;
      #OPTION: unstranded data
unstr_SE|unstr_PE)
      STARparStrand="--outSAMstrandField intronMotif"
      STARparWig="--outWigStrand Unstranded"
      ;;
esac

###### STAR command
echo $STAR $STARparCommon $STARparCircRNA $STARparRun $STARparBAM $STARparStrand
$STAR $STARparCommon $STARparCircRNA $STARparRun $STARparBAM $STARparStrand

###### bedGraph generation, now decoupled from STAR alignment step
cd $outdir

mkdir Signal

echo $STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
$STAR --runMode inputAlignmentsFromBAM   --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr

# move the signal files from the subdirectory
mv Signal/Signal*bg .




###### bigWig conversion commands
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

case "$dataType" in
str_SE|str_PE)
      # stranded data
      str[1]=-; str[2]=+;
      for istr in 1 2
      do
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str$istr.out.bg > sig.tmp
          bedSort sig.tmp sig.tmp
          $bedGraphToBigWig sig.tmp  chrNL.txt ${filename}.$imult.strand${str[istr]}.bw
      done
      done
      ;;
unstr_SE|unstr_PE)
#      # unstranded data
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str1.out.bg > sig.tmp
          bedSort sig.tmp sig.tmp
          $bedGraphToBigWig sig.tmp chrNL.txt  ${filename}.$imult.unstranded.bw
      done
      ;;
esac




######### RSEM

#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
trBAMsortRAM=55G

mv Aligned.toTranscriptome.out.bam Tr.bam 

case "$dataType" in
str_SE|unstr_SE)
      # single-end data
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | sort -S $trBAMsortRAM -T ./ ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.bam
      ;;
str_PE|unstr_PE)
      # paired-end data, merge mates into one line before sorting, and un-merge after sorting
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S $trBAMsortRAM -T ./ | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > Aligned.toTranscriptome.out.bam
      ;;
esac

'rm' Tr.bam


# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM --ci-memory 30000 "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac


###### RSEM command
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant >& Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType Aligned.toTranscriptome.out.bam $RSEMrefDir Quant >& Log.rsem

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file Quant.pdf, which contains multiple plots
echo rsem-plot-model Quant Quant.pdf
rsem-plot-model Quant Quant.pdf


###### Rename files and index BAM
mv Aligned.toTranscriptome.out.bam ${filename}.toTranscriptome.bam
mv Aligned.sortedByCoord.out.bam ${filename}.sortedByCoord.bam
mv Quant.genes.results ${filename}.genes.results
mv Quant.isoforms.results ${filename}.isoforms.results
mv Quant.pdf ${filename}.quant.pdf
mv Chimeric.out.junction ${filename}.chimeric.out.junction
mv Chimeric.out.sam ${filename}.chimeric.out.sam
samtools index ${filename}.toTranscriptome.bam
samtools index ${filename}.sortedByCoord.bam
rm *.fq

###### RSEM counts matrix generation: Combinding results from each sample/control into one data matrix for downstream analysis.

rsem-generate-data-matrix Cnt1.genes.results  Cnt2.genes.results  Cnt3.genes.results  Cond1.genes.results  Cond1.genes.results  Cond3.genes.results > rsem_prop_counts.matrix

