#!/usr/bin/bash 
#SBATCH --mem=64G --nodes 1 --ntasks 2 --out logs/snpEff.log

module unload perl
module load bcftools/1.10
module load snpEff
module load tabix

SNPEFFOUT=snpEff
SNPEFFGENOME=Homalodisca_vitripennis.A6A7A9_masurca_v1
snpEffConfig=snpEff.config
GFFGENOME=Homalodisca_vitripennis_A6A7A9_masurca_v1.gff3
MEM=64g

# this module defines SNPEFFJAR and SNPEFFDIR
if [ -f config.txt ]; then
	source config.txt
fi
GFFGENOMEFILE=$GENOMEFOLDER/$GFFGENOME
FASTAGENOMEFILE=$GENOMEFOLDER/$GENOMEFASTA

if [ -z $SNPEFFJAR ]; then
 echo "need to defined \$SNPEFFJAR in module or config.txt"
 exit
fi
if [ -z $SNPEFFDIR ]; then
 echo "need to defined \$SNPEFFDIR in module or config.txt"
 exit
fi
# could make this a confi

if [ -z $FINALVCF ]; then
	echo "need a FINALVCF variable in config.txt"
	exit
fi

mkdir -p $SNPEFFOUT
if [ ! -e $SNPEFFOUT/$snpEffConfig ]; then
	rsync -a $SNPEFFDIR/snpEff.config $SNPEFFOUT/$snpEffConfig
	echo "# Hvit " >> $SNPEFFOUT/$snpEffConfig
  	echo "$SNPEFFGENOME.genome : Homalodisca_vitripennis.A6A7A9 masurca 1 " >> $SNPEFFOUT/$snpEffConfig
	chroms=$(grep '##sequence-region' $GFFGENOMEFILE | awk '{print $2}' | perl -p -e 's/\n/, /' | perl -p -e 's/,\s+$/\n/')
	echo -e "\t$SNPEFFGENOME.chromosomes: $chroms" >> $SNPEFFOUT/$snpEffConfig
	mkdir -p $SNPEFFOUT/data/$SNPEFFGENOME
	gzip -c $GFFGENOMEFILE > $SNPEFFOUT/data/$SNPEFFGENOME/genes.gff.gz
	rsync -aL $REFGENOME $SNPEFFOUT/data/$SNPEFFGENOME/sequences.fa

	java -Xmx$MEM -jar $SNPEFFJAR build -datadir `pwd`/$SNPEFFOUT/data -c $SNPEFFOUT/$snpEffConfig -gff3 -v $SNPEFFGENOME
fi
pushd $SNPEFFOUT
COMBVCF="../$FINALVCF/$PREFIX.SNP.combined_selected.vcf.gz ../$FINALVCF/$PREFIX.INDEL.combined_selected.vcf.gz"
for n in $COMBVCF
do
 echo $n
 st=$(echo $n | perl -p -e 's/\.gz//')
 if [ ! -f $n ]; then
	 bgzip $st
	 tabix $n
 fi
done
INVCF=$PREFIX.combined_all.vcf
OUTVCF=$PREFIX.snpEff.vcf
OUTTAB=$PREFIX.snpEff.tab
if [[ ! -f $INVCF || ]; then
	bcftools concat -a -d both -o $INVCF -O v $COMBVCF
fi
java -Xmx$MEM -jar $SNPEFFJAR eff -dataDir `pwd`/data -v $SNPEFFGENOME $INVCF > $OUTVCF

bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\t%INFO/ANN\n' $OUTVCF > $OUTTAB

#../scripts/map_snpEff2domains.py --vcf $OUTVCF --domains ../genome/FungiDB-39_AfumigatusAf293_InterproDomains.txt --output A_fumigiatus_Af293.Popgen8.snpEf.domain_variants.tsv

../scripts/snpEff_2_tab.py $OUTVCF > $PREFIX.matrix.tab
