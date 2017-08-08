#!/bin/bash
#set -x

###############################################################
#
#  Curate TFBS gene sets from human ENCODE data
#  by Mark Ziemann 2016 mark.ziemann@gmail.com
#
###############################################################

# Description: Will download Encode TFBS peak sets and will
# annotate TSS and enhancers bound by transcription factors
# and create a gene matrix that can be used in GSEA and other
# pathway analysis tools

###############################################################
echo "Testing dependancies"
###############################################################

PARALLEL=`which parallel | wc -l`
if [ $PARALLEL -eq "0" ] ; then
  echo 'Error: gnu parallel was not found.
It can be downloaded from the following URL:
https://www.gnu.org/software/parallel/
Also available from the Ubuntu software centre:
sudo apt-get install parallel'
  exit
fi

LIFTOVER=`which liftOver | wc -l`
if [ $LIFTOVER -eq "0" ] ; then
  echo 'Error: liftOver was not found.
It can be downloaded from the following URL:
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver'
  exit
fi

BEDTOOLS=`which bedtools | wc -l`
if [ $BEDTOOLS -eq "0" ] ; then
  echo 'Error: bedtools was not found.
It can be downloaded from the following URL:
http://bedtools.readthedocs.io/en/latest/
Older version is available from the Ubuntu software centre:
sudo apt-get install bedtools'
  exit
fi

###############################################################
# Lets declare a bunch of vars for later
###############################################################
GTFURL=ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
GTF=`basename $GTFURL`
TSSBED=`echo $GTF | sed 's/.gtf.gz/.tss.bed/'`
#SIZE_RANGE='40000,30000,25000,20000,17500,15000,12500,10000,9000,8000,7000,6000,5000,4000,3000,2000,1000,500,250,100'
SIZE_RANGE='5000,2500,1000,500'
CHAINURL=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
CHAIN=`basename $CHAINURL`
ENHANCERURL=https://genecards.weizmann.ac.il/geneloc/gh.gff
#ENHANCERURL=ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/dhs_gene_connectivity/genomewideCorrs_above0.7_promoterPlusMinus500kb_withGeneNames_32celltypeCategories.bed8.gz
ENHANCERZIP=`basename $ENHANCERURL`
GNAMES=$(basename $GTF .gtf.gz).gnames.txt
ENH_BED19=enhancers.hg19.bed
ENH_BED38=enhancers.hg38.bed
METADATA_URL="https://www.encodeproject.org/metadata/type=Experiment&assay_title=ChIP-seq&target.investigated_as=transcription+factor&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+narrowPeak/metadata.tsv"
METADATA=metadata.tsv
METADATA_SUMMARY=metadata_summary.tsv
URLLIST=files.txt

#NUMRANGE='100,250,500,1000,2000,3000,5000,100000'
NUMRANGE='250,500,1000,2000'
UPPER_LIMIT='5000'
TSSGMT=human_tss_TFBS.gmt
ENHGMT=human_enhancer_TFBS.gmt
NRCPU=`nproc`

###############################################################
echo "Dependancies OK. Downloading required annotation files"
###############################################################
wget -N $GTFURL
wget -N $CHAINURL
wget -N $ENHANCERURL

###############################################################
echo "extracting gene names from GTF"
###############################################################
zcat $GTF | grep -w gene | cut -f9 | cut -d '"' -f2,6 \
| tr '"' '\t' | sort -uk 2b,2 > $GNAMES

###############################################################
echo "extracting enhancer coordinates"
###############################################################

myfunc(){
LINE="$*"
COORD=$(echo $LINE | cut -d ' ' -f1)
echo $LINE | cut -d ' ' -f2- | tr ' ' '\n' | sed "s/^/${COORD}\t/;s/:/\t/;s/-/\t/"
}
export -f myfunc

sed 1d gh.gff  | tr '=' '\t' \
| awk -F"\t" '{print $1,$4,$5,$NF}' \
| sed 's/ /:/;s/ /-/;s/,//g' \
| parallel -X myfunc \
| sed 's/chr//' | bedtools sort > $ENH_BED38

grep -v ENSG $ENH_BED38  | sort -k 4b,4 \
| join -1 4 -2 2 - $GNAMES | tr ' ' '\t' | cut -f2- > $ENH_BED38.tmp

grep ENSG $ENH_BED38 >> $ENH_BED38.tmp

bedtools sort -i $ENH_BED38.tmp > $ENH_BED38 && rm $ENH_BED38.tmp

###############################################################
echo "extracting TSS positions from ensembl GTF"
###############################################################
#Get coordinates of plus strand genes
zgrep -w 'exon_number "1"' $GTF | awk '$7=="+"' \
| cut -f1,4,5,7,9 | cut -d '"' -f-2,12 \
| sed 's!gene_id "!!' | tr '"' '_'  \
| awk '{OFS="\t"} {print $1,$2,$2+1,$5,$4}' \
| awk '{OFS="\t"} { if ($2<1) print $1,"1",$3,$5,$4 ; else print $0 }' > $TSSBED

#Get coordinates of minus strand genes
zgrep -w 'exon_number "1"' $GTF | awk '$7=="-"' \
| cut -f1,4,5,7,9 | cut -d '"' -f-2,12 \
| sed 's!gene_id "!!' | tr '"' '_'  \
| awk '{OFS="\t"} {print $1,$3-1,$3,$5,$4}' \
| awk '{OFS="\t"} { if ($2<1) print $1,"1",$3,$5,$4 ; else print $0 }' >> $TSSBED

###############################################################
echo "download metadata and peak data"
###############################################################

wget -O $METADATA "https://www.encodeproject.org/metadata/type=Experiment&assay_title=ChIP-seq&target.investigated_as=transcription+factor&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+narrowPeak/metadata.tsv"

grep released $METADATA | grep GRCh38 | egrep '(optimal|psudorep)' \
| cut -f1,7,8,13,17,25,42 | tr ' ' '_' | sed 's/\t/ /2;s/\t/ /2;s/\t/ /2' \
| sed 's/  / untreated /' | sort -k6r | awk -F"\t" '!arr[$2,$3,$4]++ ' \
| tr ' ' '\t' | sed 's/$/\tGRCh38/'  > $METADATA_SUMMARY

cut -f7 $METADATA_SUMMARY > $URLLIST

if [ ! -d "peakfiles" ] ; then
  mkdir peakfiles
else
  mv peakfiles/*bed.gz .
fi

for URL in $(cat $URLLIST) ; do
  BASENAME=$(basename $URL)
  FILENAME=$BASENAME
  if [ ! -r $BASENAME ] ; then
    wget -N --quiet -O $FILENAME $URL
  fi

  REF_MD5SUM=$(grep -w $URL $METADATA | cut -f40)
  MY_MD5SUM=$(md5sum $FILENAME | awk '{print $1}')

  if [ $REF_MD5SUM != $MY_MD5SUM ] ; then
    rm $FILENAME
    wget -N --quiet -O $FILENAME $URL
    MY_MD5SUM=$(md5sum $FILENAME | awk '{print $1}')

  else
      mv $FILENAME peakfiles
  fi
done

if [ ! -d "trash" ] ; then
  mkdir trash
fi

mv *bed.gz trash

###############################################################
echo "convert hg19 files to GRCh38"
###############################################################
comm -23 <(grep hg19 $METADATA | cut -f7,17 | sort -u) <(grep GRCh38 $METADATA | cut -f7,17 | sort -u) \
| tr '\t' '@' \
| while read line ; do
  CELL=$(echo $line | cut -d '@' -f1)
  TF=$(echo $line | cut -d '@' -f2)
  grep hg19 $METADATA | grep "$CELL" | grep "$TF"
done | grep released | egrep '(pseudo|optimal)' \
| cut -f1,7,8,13,17,25,42 | tr ' ' '_' | sed 's/\t/ /2;s/\t/ /2;s/\t/ /2' \
| sed 's/  / untreated /' | sort -k6r | awk -F"\t" '!arr[$2,$3,$4]++ ' \
| tr ' ' '\t' | sed 's/$/\thg19/'  >> $METADATA_SUMMARY

grep hg19 $METADATA_SUMMARY | cut -f7 > $URLLIST

if [ ! -d "hg19_peakfiles" ] ; then
  mkdir hg19_peakfiles
else
  mv hg19_peakfiles/*bed.gz .
fi

for URL in $(cat $URLLIST) ; do
  BASENAME=$(basename $URL)
  FILENAME=$BASENAME
  if [ ! -r $BASENAME ] ; then
    wget -N --quiet -O $FILENAME $URL
  fi

  REF_MD5SUM=$(grep -w $URL $METADATA | cut -f40)
  MY_MD5SUM=$(md5sum $FILENAME | awk '{print $1}')

  if [ $REF_MD5SUM != $MY_MD5SUM ] ; then
    rm $FILENAME
    wget -N --quiet -O $FILENAME $URL
    MY_MD5SUM=$(md5sum $FILENAME | awk '{print $1}')

  else
      mv $FILENAME hg19_peakfiles
  fi
done

run_liftover(){
set -x
HG19PK=$1
CHAIN=$2
BASE=$(basename $HG19PK)
liftOver <(zcat $HG19PK | cut -f-3,7) $CHAIN /dev/stdout unmapped \
| sed 's/\t/\t.\t.\t.\t/3' | gzip > peakfiles/$BASE
}
export -f run_liftover
parallel -X run_liftover ::: hg19_peakfiles/*gz ::: $CHAIN

###############################################################
echo "Intersect TFBS with TSS and distal elements"
###############################################################
if [ ! -d "GMT" ] ; then
  mkdir GMT
else
  rm -rf GMT
  mkdir GMT
fi

echo map2
map2(){
set -x
PKZ_PATH=$1
PKZ=`basename $PKZ_PATH`
PK=`basename $PKZ_PATH .bed.gz`
METADATA_SUMMARY=$2
TSSBED=$3
ENH_BED38=$4
SIZE_RANGE=$5
MD5SUM=$(md5sum $PKZ| awk '{print $1}')

if [ $(grep -c $PK $METADATA_SUMMARY ) -eq "1" ] ; then

TF=$(grep -w $PK $METADATA_SUMMARY | cut -f5)
CELL=$(grep -w $PK $METADATA_SUMMARY | cut -f2,3 | sed 's/\t/_/' )
TRT=$(grep -w $PK $METADATA_SUMMARY | cut -f4 )
BASE=${TF}_${CELL}_${TRT}_EncodeDataset:${PK}

for UPSTREAM in `echo $SIZE_RANGE | tr ',' ' '` ; do
  #for DOWNSTREAM in `echo $SIZE_RANGE | tr ',' ' '` ; do
  DOWNSTREAM=$UPSTREAM
    GMT=GMT/$PK.tss.${UPSTREAM}bpUpstream.${DOWNSTREAM}bpDownstream.gmt

    bedtools intersect -wb \
    -a <(zcat $PKZ_PATH | sed 's/chr//' | sort -k7gr) \
    -b <(awk -v U=$UPSTREAM -v D=$DOWNSTREAM '{OFS="\t"} { if ($5=="+") print $1,$2-U,$3+D,$4 ; else print $1,$2-D,$3+U,$4}' $TSSBED \
    | awk '{OFS="\t"} { if ($2<1) print $1,"1",$3,$4 ; else print $0 }') \
    | awk '{print $NF}' | cut -d '_' -f1 \
    | awk '!arr[$1]++'  \
    | tr '\n' '\t' | sed 's!$!\n!' \
    | sed "s!^!${BASE}\tENCODE_Dataset:${PK}_TFBS_at_TSS\t!"  > $GMT

#  done

  GMT=GMT/$PK.enh.${UPSTREAM}bpFlanking.gmt

  bedtools intersect -wb \
  -a <(zcat $PKZ_PATH | sed 's/chr//' | sort -k7gr) \
  -b <(awk -v U=$UPSTREAM '{OFS="\t"} {print $1,$2-U,$3+U,$4}' $ENH_BED38 \
  | awk '{OFS="\t"} { if ($2<1) print $1,"1",$3,$4 ; else print $0 }' ) \
  | awk '{print $NF}' | awk '!arr[$1]++' \
  | tr '\n' '\t' | sed 's!$!\n!' \
  | sed "s!^!${BASE}\tENCODE_Dataset:${PK}_TFBS_at_distalelements\t!" > $GMT
done

fi
}
export -f map2
parallel -j $NRCPU map2 \
::: $(ls peakfiles/*bed.gz) \
::: $METADATA_SUMMARY \
::: $TSSBED \
::: $ENH_BED38 \
::: $SIZE_RANGE

###############################################################
echo "Trim down gene sets to the desired size"
###############################################################
for SET in `ls GMT/ | grep '.enh.' | awk -F"." '{print $(NF-1)"."$NF}' | sort -u ` ; do
  GMT=`echo $ENHGMT | sed "s#.gmt#.${SET}#" `
  echo $SET
  cat GMT/*.$SET | tr '/' '_' | cut -f-$((UPPER_LIMIT+2)) > $GMT
done

for SET in `ls GMT/ | grep tss | awk -F"." '{print $(NF-2)"."$(NF-1)"."$NF}' | sort -u ` ; do
  GMT=`echo $TSSGMT | sed "s#.gmt#.${SET}#" `
  echo $SET
  cat GMT/*.$SET | tr '/' '_' | cut -f-$((UPPER_LIMIT+2)) > $GMT
done


###############################################################
echo "Trim down gene sets to the desired size"
###############################################################
for GMT in `ls *gmt | grep -v genes.gmt` ; do
  for NUM in `echo $NUMRANGE | tr ',' ' '` ; do
    GMT2=`echo $GMT | sed "s#.gmt#.${NUM}genes.gmt#"`
    NUM=$((NUM+2))
    cut -f-$NUM $GMT > $GMT2
  done
done


###############################################################
echo "Translate Ensembl IDs into gene names"
###############################################################
GNAMES=$(basename $GTF .gtf.gz).gnames.txt
GENE_SYMBOLS_DIR=gene_symbols

zgrep -w gene $GTF | cut -d '"' -f2,6 \
| tr '"' '\t' | sort -k 1b,1 > $GNAMES

myfunc(){
GNAMES=*.gnames.txt
LINE="$*"
PFX=$(echo $LINE | cut -d ' ' -f-2)
GENES=$(echo $LINE | cut -d ' ' -f3- | tr ' ' '\n' | sort -k 1b,1 \
| join -1 1 -2 1 - $GNAMES  | cut -d ' ' -f2 | tr '\n' ' ' )
echo $PFX $GENES | tr ' ' '\t'
}
export -f myfunc

if [ ! -d $GENE_SYMBOLS_DIR ] ; then
  mkdir $GENE_SYMBOLS_DIR
else
  rm $GENE_SYMBOLS_DIR/*
fi

for GMT in *gmt ; do
  >$GENE_SYMBOLS_DIR/$GMT
  cat $GMT | parallel myfunc >> $GENE_SYMBOLS_DIR/$GMT
done



