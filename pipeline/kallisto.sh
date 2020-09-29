#!/usr/bin/bash
#SBATCH -p short --mem 48gb -N 1 -n 16  --out logs/kallisto.log

module load kallisto/0.46.2

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=fastq
OUTDIR=results/kallisto
IDX=Hvit_masurca_v1.idx
TX=genome/Homalodisca_vitripennis_A6A7A9_masurca_v1.mrna-transcripts.fa
SAMPLEFILE=samples.tsv
mkdir -p $OUTDIR
if [ ! -f $IDX ]; then
    kallisto index -i $IDX $TX
fi

tail -n +2 $SAMPLEFILE |  while read SAMPLE NAME REP LOCATION STATUS
do
 OUTNAME=$NAME.${REP}
 if [ ! -f $OUTDIR/$OUTNAME/abundance.h5 ]; then
     kallisto quant -i $IDX -o $OUTDIR/$OUTNAME -t $CPU --bias $INDIR/${SAMPLE}_R1.fastq.gz $INDIR/${SAMPLE}_R2.fastq.gz
 fi
done
