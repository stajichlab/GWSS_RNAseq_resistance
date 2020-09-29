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

TYPE=SNP
OUT=$FINALVCF/$PREFIX.combined_selected.vcf.gz
IN=$SLICEVCF/$PREFIX
bcftools concat -Oz -o $OUT --threads $CPU $SLICEVCF/$PREFIX.*.$TYPE.selected.vcf.gz
