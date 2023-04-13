#!/bin/bash
#BSUB -n 4
#BSUB -q priority
#BSUB -R "span[ptile=2]"
#BSUB -M 25000
#BSUB -o UTRs.out
#BSUB -e UTRs.err
#BSUB -J UTRs


module load parallel/20220722
module load samtools
module load bedtools
module load R/4.2.2

time

RUNID=""

while  test  $# -gt 0;  do
        case "$1"  in
                -runid)
                        shift
                        if test $# -gt 0; then
                                RUNID="$1"
                        else
                                echo -e "\n\nERROR! runid not specified.\n\n"
                                exit 1
                        fi
                        shift
                        ;;
        esac
done

SPECIES="m39"
GTFtag="Mus_musculus.GRCm39.104.gtf"
ORFDIR="/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/$SPECIES"
ORFTAG="utrs.$SPECIES.$GTFtag"
ORFPREF="$ORFDIR/$ORFTAG"

TMPWF="$RUNID"
TMPOF="$RUNID"


if [ ! -z "$RUNID" ]; then
        echo -e "runid is given.  inferring all data...\n"

        $RPROF/src/retrieve_subset_from_tabular_database.sh  -r $RUNID  -c project,individual,runid,species,layout,strategy  -o $TMPWF.info.sub

        tail  -n +2  $TMPWF.info.sub  >  $TMPWF.info.sub.tmp        \
        &&  mv  $TMPWF.info.sub.tmp   $TMPWF.info.sub

        PROJECT=$(cat  $TMPWF.info.sub  |  cut -f 1)
        INDIV=$(cat  $TMPWF.info.sub  |  cut -f 2)
        SPECIES=$(cat  $TMPWF.info.sub  |  cut -f 4)
        LAYOUT=$(cat  $TMPWF.info.sub  |  cut -f 5)
        LAYOUT=$( echo "$LAYOUT"  |  tr '[:lower:]' '[:upper:]' )
        STRATEGY=$(cat  $TMPWF.info.sub  |  cut -f 6)

        RSEMDIR="$RPROF/out/$PROJECT/$RUNID"
        RSEMTAG="STAR.transcript"
        RSEMBAM="$RSEMDIR/$RSEMTAG"



        echo -e "\n\n\n"
        echo -e "PROJECT=[$PROJECT]"
        echo -e "INDIV=[$INDIV]"
        echo -e "SPECIES=[$SPECIES]"
        echo -e "LAYOUT=[$LAYOUT]"
        echo -e "STRATEGY=[$STRATEGY]"
        echo -e "RSEMDIR=[$RSEMDIR]"
        echo -e "RSEMTAG=[$RSEMTAG]"
        echo -e "RSEMBAM=[$RSEMBAM]"
        echo -e "\n\n\n"

fi # inferif

echo -e "PROJECT=[$ORFPREF.utrs.bed.gz]\n"


echo -e "Local name of BAM TMPWF=[$TMPWF.bam]"
cp $RSEMBAM.bam $TMPWF.bam


echo -e "\t samtool index ..\n"
samtools index $TMPWF.bam


#echo -e "\t samtool sort ..\n"
#samtools sort -@ 4 -n -o $TMPWF.bam $TMPWF.bam


echo -e "\t samtool view ..\n"
samtools view -@ 4 -hb -o $TMPWF.rsem.map.bam $TMPWF.bam



echo -e "\t samtool view and bedtool bamtobed..\n"
samtools view -@ 4 $TMPWF.rsem.map.bam | grep -oP 'ZW:f:[^\t]+' | cut -d ':' -f 3 > $TMPWF.zw
bedtools bamtobed -i $TMPWF.rsem.map.bam| cut -f 1-3 | paste - $TMPWF.zw > $TMPWF.rsem.bam.bed
awk -v f=$OFFSET  '{ $3 = $2 + 1 + f ; print $0 }'  FS="\t"  OFS="\t"  $TMPWF.rsem.bam.bed  >  $TMPWF.rsem.bam.bed.pt


echo -e "sort bambed...\n"
sort -t $'\t' -k1,1 -k2,2n --parallel=4 -o $TMPWF.rsem.bam.bed.pt $TMPWF.rsem.bam.bed.pt


echo -e "bed map...\n"
echo -e "\tmap to each orf...\n"


bedtools map -c 4 -o sum -null 0 -a $ORFTAG.utrs.bed -b $TMPWF.rsem.bam.bed.pt > $TMPWF.bed

time
