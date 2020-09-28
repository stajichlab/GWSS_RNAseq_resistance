#!/usr/bin/bash
#SBATCH --mem 64G --nodes 1 --ntasks 8 -J GATK.GVCFGeno --out logs/GVCFGenoGATK4.log --time 7-0:00:00

MEM=64g
module unload java
module load java/8
module load picard
module load gatk/4
module load bcftools
module load parallel
TEMP=/scratch
FINALVCF=vcf
mkdir -p $FINALVCF
if [ -f config.txt ]; then
	source config.txt
fi
OUT=$FINALVCF/$PREFIX.all.vcf

if [ ! -f $REFGENOME ]; then
    module load samtools/1.9
    samtools faidx $REFGENOME
fi
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
    CPU=2
fi
if [[ $(ls $GVCFFOLDER | grep -c -P "\.g.vcf$") -gt "0" ]]; then
    parallel -j $CPU bgzip {} ::: $GVCFFOLDER/*.g.vcf
    parallel -j $CPU tabix -f {} ::: $GVCFFOLDER/*.g.vcf.gz
fi

FILES=$(ls $GVCFFOLDER/*.g.vcf.gz | sort | perl -p -e 's/(\S+)\n/-V $1 /')
INTERVALS=$(cut -f1 $REFGENOME.fai  | perl -p -e 's/(\S+)\n/--intervals $1 /g')
mkdir /scratch/$USER
DB=/scratch/$USER/${GVCFFOLDER}_db
if [ ! -f $OUT.gz ]; then
    if [ ! -f $OUT ]; then
	rm -rf $DB
        if [ ! -d $DB ]; then
		gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --consolidate --merge-input-intervals --genomicsdb-workspace-path $DB $FILES $INTERVALS --reader-threads $CPU --tmp-dir $TEMP
	#gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --genomicsdb-workspace-path $DB $FILES $INTERVALS  --reader-threads $CPU
	time gatk GenotypeGVCFs --reference $REFGENOME --output $OUT -V gendb://$DB
    fi
    if [ -f $OUT ]; then
    	bgzip $OUT
    	tabix $OUT.gz
    fi
fi
