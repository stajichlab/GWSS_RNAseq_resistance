#!/usr/bin/bash
#SBATCH --mem 24G --nodes 1 --ntasks 2 -J slice.GVCFGeno --out logs/GVCFGenoGATK4.slice_%a.log --time 48:00:00

MEM=24g
module unload java
module load picard
module load gatk/4
module load java/13
module load bcftools
module load parallel
TEMP=/scratch/$USER
FINALVCF=vcf_slices
INTERVAL=50
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
mkdir -p $FINALVCF
if [ -f config.txt ]; then
	source config.txt
fi
OUT=$FINALVCF/$PREFIX.$N.all.vcf
FILTERSNP=$FINALVCF/$PREFIX.$N.SNP.filter.vcf
FILTERINDEL=$FINALVCF/$PREFIX.$N.INDEL.filter.vcf
SELECTSNP=$FINALVCF/$PREFIX.$N.SNP.selected.vcf
SELECTINDEL=$FINALVCF/$PREFIX.$N.INDEL.selected.vcf

if [ ! -f $REFGENOME ]; then
    module load samtools/1.9
    samtools faidx $REFGENOME
fi
NSTART=$(perl -e "printf('%d',1 + $INTERVAL * ($N - 1))")
NEND=$(perl -e "printf('%d',$INTERVAL * $N)")
MAX=$(wc -l $REFGENOME.fai | awk '{print $1}')
if [ "$NSTART" -gt "$MAX" ]; then
	echo "NSTART ($NSTART) > $MAX"
	exit
fi
if [ "$NEND" -gt "$MAX" ]; then
	NEND=$MAX
fi
echo "$NSTART -> $NEND"

CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
    CPU=2
fi
if [[ $(ls $GVCFFOLDER | grep -c -P "\.g.vcf$") -gt "0" ]]; then
    parallel -j $CPU bgzip {} ::: $GVCFFOLDER/*.g.vcf
    parallel -j $CPU tabix -f {} ::: $GVCFFOLDER/*.g.vcf.gz
fi

FILES=$(ls $GVCFFOLDER/*.g.vcf.gz | sort | perl -p -e 's/(\S+)\n/-V $1 /')
INTERVALS=$(cut -f1 $REFGENOME.fai  | sed -n "${NSTART},${NEND}p" | perl -p -e 's/(\S+)\n/--intervals $1 /g')
mkdir -p $TEMP $FINAL
DB=/$TEMP/${GVCFFOLDER}_slice_$N
if [ ! -f $OUT.gz ]; then
    if [ ! -f $OUT ]; then
	rm -rf $DB
	gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --consolidate --merge-input-intervals --genomicsdb-workspace-path $DB $FILES $INTERVALS --reader-threads $CPU --tmp-dir $TEMP
	#gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --genomicsdb-workspace-path $DB $FILES $INTERVALS  --reader-threads $CPU
	time gatk GenotypeGVCFs --reference $REFGENOME --output $OUT -V gendb://$DB
	rm -rf $DB
    fi
    if [ -f $OUT ]; then
    	bgzip $OUT
    	tabix $OUT.gz
    fi
fi


TYPE=SNP
echo "VCF = $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz"
if [[ ! -f $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz ]]; then
    gatk SelectVariants \
	-R $REFGENOME \
	--variant $OUT.gz \
	-O $FINALVCF/$PREFIX.$N.$TYPE.vcf \
	--restrict-alleles-to BIALLELIC \
	--select-type-to-include $TYPE --create-output-variant-index false
    bgzip $FINALVCF/$PREFIX.$N.$TYPE.vcf
    tabix $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz
fi

if [[ ! -f $FILTERSNP.gz || $FINALVCF/$$PREFIX.$N.$TYPE.vcf.gz -nt $FILTERSNP.gz ]]; then
    gatk VariantFiltration --output $FILTERSNP \
	--variant $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz -R $REFGENOME \
	--cluster-window-size 10  \
	--filter-expression "QD < 2.0" --filter-name QualByDepth \
	--filter-expression "MQ < 40.0" --filter-name MapQual \
	--filter-expression "QUAL < 100" --filter-name QScore \
	--filter-expression "SOR > 4.0" --filter-name StrandOddsRatio \
	--filter-expression "FS > 60.0" --filter-name FisherStrandBias \
	--missing-values-evaluate-as-failing --create-output-variant-index false

#	--filter-expression "MQRankSum < -12.5" --filter-name MapQualityRankSum \
#	--filter-expression "ReadPosRankSum < -8.0" --filter-name ReadPosRank \

    bgzip $FILTERSNP
    tabix $FILTERSNP.gz
fi

if [[ ! -f $SELECTSNP.gz || $FILTERSNP.gz -nt $SELECTSNP.gz ]]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $FILTERSNP.gz \
	--output $SELECTSNP \
	--exclude-filtered --create-output-variant-index false
    bgzip $SELECTSNP
    tabix $SELECTSNP.gz
fi

TYPE=INDEL
if [ ! -f $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz ]; then
    gatk SelectVariants \
        -R $REFGENOME \
        --variant $OUT.gz \
        -O $FINALVCF/$PREFIX.$N.$TYPE.vcf  --select-type-to-include MIXED --select-type-to-include MNP \
        --select-type-to-include $TYPE --create-output-variant-index false
    bgzip $FINALVCF/$PREFIX.$N.$TYPE.vcf
    tabix $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz
fi

if [[ ! -f $FILTERINDEL.gz || $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz -nt $FILTERINDEL.gz ]]; then
    gatk VariantFiltration --output $FILTERINDEL \
	--variant $FINALVCF/$PREFIX.$N.$TYPE.vcf.gz -R $REFGENOME \
	--cluster-window-size 10  -filter "QD < 2.0" --filter-name QualByDepth \
	-filter "SOR > 10.0" --filter-name StrandOddsRatio \
	-filter "FS > 200.0" --filter-name FisherStrandBias \
	-filter "InbreedingCoeff < -0.8" --filter-name InbreedCoef \
	--create-output-variant-index false

#	-filter "ReadPosRankSum < -20.0" --filter-name ReadPosRank \
#	-filter "MQRankSum < -12.5" --filter-name MapQualityRankSum \

    bgzip $FILTERINDEL
    tabix $FILTERINDEL.gz
fi

if [[ ! -f $SELECTINDEL.gz || $FILTERINDEL.gz -nt $SELETINDEL.gz ]]; then
    gatk SelectVariants -R $REFGENOME \
	--variant $FILTERINDEL.gz \
	--output $SELECTINDEL \
	--exclude-filtered --create-output-variant-index false
    bgzip $SELECTINDEL
    tabix $SELECTINDEL.gz
fi
