#! /bin/bash    

fq=$1

# location to genome reference FASTA file
     genome=/groups/hbctraining/unix_oct2015_other/reference_STAR/
     gtf=~/unix_oct2015/rnaseq_project/data/reference_data/chr1-hg19_genes.gtf

# grab base of filename for future naming
    base=$(basename $fq .qualtrim25.minlen35.fq)

# set up output filenames and locations
    align_out=~/unix_oct2015/rnaseq_project/results/STAR/${base}_Aligned.bam
    align_in=~/unix_oct2015/rnaseq_project/results/STAR/${base}_Aligned.sortedByCoord.out.bam
    counts=~/unix_oct2015/rnaseq_project/results/counts/${base}.counts

echo "#! /bin/bash"
echo "#BSUB -q priority       # Partition to submit to (comma separated)"
echo "#BSUB -n 6                  # Number of cores"
echo "#BSUB -W 1:30               # Runtime in D-HH:MM (or use minutes)"
echo "#BSUB -R "rusage[mem=4000]"    # Memory in MB"
echo "#BSUB -J rnaseq_mov10         # Job name"
echo "#BSUB -o %J.out       # File to which standard out will be written"
echo "#BSUB -e %J.err       # File to which standard err will be written"

echo ""
echo "# set up our software environment..."
echo "module load seq/samtools"
echo "module load seq/htseq"

echo ""
echo "# Run STAR"
echo "STAR --runThreadN 6 --genomeDir $genome --readFilesIn $fq --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS"

echo ""
echo "# Create BAM index"
echo "samtools index $align_in"

echo ""
echo "# Count mapped reads"
echo "htseq-count --stranded reverse --format bam $align_in $gtf  >  $counts"


