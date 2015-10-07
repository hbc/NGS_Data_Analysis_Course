#! /bin/bash    


fq=$1

# location to genome reference FASTA file
     genome=/groups/hbctraining/unix_oct2015_other/reference_STAR/
     gtf=data/reference_data/chr1-hg19_genes.gtf

# set up our software environment...
    module load seq/samtools
    module load seq/htseq

# make all of our output directories
     mkdir results/STAR
     mkdir results/counts


echo "Processing file $fq ..."

# grab base of filename for future naming
    base=$(basename $fq .qualtrim25.minlen35.fq)
    echo "basename is $base"


# set up output filenames and locations
    align_out=results/STAR/${base}_
    align_in=results/STAR/${base}_Aligned.sortedByCoord.out.bam
    counts=results/counts/${base}.counts

# Run STAR
STAR --runThreadN 6 --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS

# Create BAM index
samtools index $align_in

# Count mapped reads
htseq-count --stranded reverse --format bam $align_in $gtf  >  $counts


