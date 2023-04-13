#!/bin/bash
#set  -eu  -o pipefail

source bashrc_rdna_variant_calling

date
SPECIES = "mm39"
GTFtag= "Mus_musculus.GRCm39.104.gtf"

ORFDIR="/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/"
ORFTAG="orfs.$SPECIES.$GTFtag.orfs"
ORFPREF="$ORFDIR/$ORFTAG"





################################# copy forward
echo -e "copy forward...\n"

rsync  -tvh  $RSEMBAM   $TMPWF.rsem.bam
rsync  -tvh  $ORFPREF.bed.gz  $TMPWF.orfs.bed.gz



echo -e "index..."

samtools  index  $TMPWF.rsem.bam



echo -e "gunzip...\n"

gunzip  $TMPWF.orfs.bed.gz






echo -e "sort orfs...\n"

sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -o $TMPWF.orfs.bed  $TMPWF.orfs.bed



echo -e "merge orfs by class...\n"

awk  '{ $1 = $1 "====" $4 ; print $1 , $2 , $3  }'  FS="\t"  OFS="\t"  $TMPWF.orfs.bed  >  $TMPWF.orfs.classchr.bed

sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -o $TMPWF.orfs.classchr.bed  $TMPWF.orfs.classchr.bed

bedtools  merge  -i $TMPWF.orfs.classchr.bed  >  $TMPWF.orfs.classchr.merge.bed

sed  -i  's/====/\t/'  $TMPWF.orfs.classchr.merge.bed

awk  '{ print $1 , $3 , $4 , $2 }'  FS="\t"  OFS="\t"  $TMPWF.orfs.classchr.merge.bed  >  $TMPWF.orfs.classchr.merge.bed.tmp	\
&&  mv   $TMPWF.orfs.classchr.merge.bed.tmp  $TMPWF.orfs.classchr.merge.bed

sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -k4,4  -o $TMPWF.orfs.classchr.merge.bed  $TMPWF.orfs.classchr.merge.bed




echo -e "mapped only...\n"

samtools  view  -hb  -o $TMPWF.rsem.map.bam  $TMPWF.rsem.bam



echo -e "bam to bed..\n"

# bamtobed with -tag ZW  does not work for some unacceptable reason.

samtools  view  $TMPWF.rsem.map.bam		|	\
grep  -oP  'ZW:f:[^\t]+'			|	\
cut -d ':' -f 3 					\
>  $TMPWF.zw

bedtools  bamtobed  -i $TMPWF.rsem.map.bam	|	\
cut -f 1-3		|	\
paste  -  $TMPWF.zw		\
>  $TMPWF.rsem.bam.bed




echo -e "bambed to point only...\n"

awk  -v f=$OFFSET  '{ $3 = $2 + 1 + f ; print $0 }'  FS="\t"  OFS="\t"  $TMPWF.rsem.bam.bed  >  $TMPWF.rsem.bam.bed.pt



echo -e "sort bambed...\n"

sort  -t $'\t'  -k1,1  -k2,2n  -o $TMPWF.rsem.bam.bed.pt   $TMPWF.rsem.bam.bed.pt




echo -e "bed map...\n"

echo -e "\tmap to each orf...\n"

bedtools  map  -c 4  -o sum  -null 0  -a $TMPWF.orfs.bed  -b $TMPWF.rsem.bam.bed.pt  >  $TMPOF.bed

echo -e "\tcollapse by classification..\n"

bedtools  map  -c 4  -o sum  -null 0  -a $TMPWF.orfs.classchr.merge.bed  -b $TMPWF.rsem.bam.bed.pt  >  $TMPOF.class.bed












echo -e "gzip...\n"

gzip  $TMPOF.*







echo -e "copy back...\n"

rsync  -tvh  $TMPOF.*   $OUTDIR













echo -e "\n\nDONE!\n\n\n\n"

date






perl -pe 's/\r\n|\n|\r/\n/g' orfs.mm39.Mus_musculus.GRCm39.104.gtf.orfs.bed >TMPWF.orfs.bed
cp TMPWF.orfs.bed orfs.mm39.Mus_musculus.GRCm39.104.gtf.orfs.bed
awk  '{ $1 = $1 "====" $4 ; print $1 , $2 , $3  }'  FS="\t"  OFS="\t"  TMPWF.orfs.bed  >  TMPWF.orfs.classchr.bed
sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -o TMPWF.orfs.classchr.bed  TMPWF.orfs.classchr.bed
merge  -i TMPWF.orfs.classchr.bed  >  TMPWF.orfs.classchr.merge.bed
sed  -i  's/====/\t/'  TMPWF.orfs.classchr.merge.bed
awk  '{ print $1 , $3 , $4 , $2 }'  FS="\t"  OFS="\t"  TMPWF.orfs.classchr.merge.bed  >  TMPWF.orfs.classchr.merge.bed.tmp	\ 
   &&  mv   TMPWF.orfs.classchr.merge.bed.tmp  TMPWF.orfs.classchr.merge.bed
sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -k4,4  -o TMPWF.orfs.classchr.merge.bed  TMPWF.orfs.classchr.merge.bed
conda activate cutadaptenv
cd /research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM
cp /research_jude/rgs01_jude/groups/blancgrp/projects/rRNA_variation/common/EMT_analysis_Nitish/Riboseq/RD0175_S46_L007/STAR.transcript.bam TMPWF.rsem.bam
#samtools  index  TMPWF.rsem.bam
samtools  sort  -@ 5  -n -T TMPWF.sam-sort.tmp  -o TMPOF.transcript.bam  TMPWF.rsem.bam
samtools  view  -hb  -o  TMPWF.rsem.map.bam  TMPOF.transcript.bam 
samtools  view  TMPWF.rsem.map.bam | grep  -oP  'ZW:f:[^\t]+' | cut -d ':' -f 3 > TMPWF.zw
bedtools  bamtobed  -i TMPWF.rsem.map.bam	| cut -f 1-3 | paste  -  TMPWF.zw	>  TMPWF.rsem.bam.bed
awk  -v f=$OFFSET  '{ $3 = $2 + 1 + f ; print $0 }'  FS="\t"  OFS="\t"  TMPWF.rsem.bam.bed  >  TMPWF.rsem.bam.bed.pt
sort  -t $'\t'  -k1,1  -k2,2n  -o TMPWF.rsem.bam.bed.pt   TMPWF.rsem.bam.bed.pt
bedtools  map  -c 4  -o sum  -null 0  -a /research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/mm39/TMPWF.orfs.bed  -b TMPWF.rsem.bam.bed.pt  >  TMPOF.bed
bedtools  map  -c 4  -o sum  -null 0  -a /research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/mm39/TMPWF.orfs.classchr.merge.bed  -b TMPWF.rsem.bam.bed.pt  >  TMPOF.class.bed
