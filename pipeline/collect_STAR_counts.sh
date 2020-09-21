#!/usr/bin/bash
#SBATCH -p short --mem 24gb -N 1 -n 16 --out logs/STARfeatureCounts.log

module load subread

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

OUTDIR=results/STAR
SAMPLEFILE=samples.tsv
GENOME=$(realpath genome/Homalodisca_vitripennis.A6A7A9_masurca_v1.masked_RModeler.fasta)
GTF=$(realpath genome/Homalodisca_vitripennis_A6A7A9_masurca_v1.gtf)
GENOME=genome/Homalodisca_vitripennis.A6A7A9_masurca_v1.masked_RModeler.fasta
GTF=genome/Homalodisca_vitripennis_A6A7A9_masurca_v1.gtf

if [ ! -s $OUTDIR/STAR_featureCounts.tsv ]; then
  featureCounts -a $GTF -o $OUTDIR/STAR_featureCounts.tsv -G $GENOME -J -g gene_id -F GTF $(find $OUTDIR -size +0 -name "*.Aligned.sortedByCoord.out.bam")
fi
