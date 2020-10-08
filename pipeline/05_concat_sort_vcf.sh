#!/usr/bin/bash
#SBATCH -p intel --mem 64gb -N 1 -n 24 --out logs/concat_vcf.log --time 24:00:00

module load bcftools

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
if [ -f config.txt ]; then
	source config.txt
else
	echo "need a config.txt"
fi

IN=$SLICEVCF/$PREFIX
for TYPE in SNP INDEL
do
OUT=$FINALVCF/$PREFIX.$TYPE.combined_selected.vcf.gz
bcftools concat -Oz -o $OUT --threads $CPU $IN.*.$TYPE.selected.vcf.gz
done
