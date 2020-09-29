#!/usr/bin/bash
#SBATCH -p short --mem 128gb -N 1 -n 32 --out logs/STAR.%a.log

module load STAR
module load subread
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=fastq
OUTDIR=results/STAR
IDX=genome/STAR
SAMPLEFILE=samples.tsv
GENOME=$(realpath genome/Homalodisca_vitripennis.A6A7A9_masurca_v1.masked_RModeler.fasta)
GFF=$(realpath ../annotation/annotate/Homalodisca_vitripennis.A6A7A9_masurca_v1/annotate_results/Homalodisca_vitripennis_A6A7A9_masurca_v1.gff3)
GTF=genome/Homalodisca_vitripennis_A6A7A9_masurca_v1.gtf
if [ ! -f $GTF ]; then
	grep -P "\texon\t" Homalodisca_vitripennis_A6A7A9_masurca_v1.gff3 | perl -p -e 's/ID=[^;]+;Parent=([^;]+);/gene_id $1/' > Homalodisca_vitripennis_A6A7A9_masurca_v1.gtf
fi
if [ ! -d $IDX ]; then
	STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $IDX --genomeFastaFiles $GENOME --sjdbGTFfile $GTF --genomeChrBinNbits 16
fi
mkdir -p $OUTDIR

tail -n +2 $SAMPLEFILE |  sed -n ${N}p | while read SAMPLE NAME REP LOCATION STATUS
do
 OUTNAME=$NAME.${REP}
 STAR --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --twopassMode Basic --runThreadN $CPU --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.08 --alignIntronMin 20 --alignIntronMax 500000 --alignMatesGapMax 500000 --genomeDir $IDX --outFileNamePrefix $OUTDIR/${OUTNAME}. --readFilesIn $INDIR/${SAMPLE}_R1.fastq.gz $INDIR/${SAMPLE}_R2.fastq.gz --readFilesCommand zcat
 
 featureCounts -a $GTF -G $GENOME -T 16 -o $OUTDIR/${OUTNAME}.featureCounts.tsv -g gene_id -J -F GTF $OUTDIR/${OUTNAME}.Aligned.sortedByCoord.out.bam
done
