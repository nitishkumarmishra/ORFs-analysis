#!/bin/bash
#BSUB -n 8
#BSUB -q priority
#BSUB -R rusage[mem=8GB]
#BSUB -o ORFs.out
#BSUB -e ORFs.err
#BSUB -J ORFs


#set  -eu  -o pipefail
# This is the shell script for rpf counts for each ORFs. Modified Matt code and change it accordingly.

source /research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/bashrc_rdna_variant_calling

echo -e "\t Activate cutadapt conda environment..\n"
source /hpcf/apps/conda3/install/5.1.0/bin/activate cutadaptenv

echo -e "\t Start time for the job\n"
time

RUNID=""
OUTPREF=""
NUMTHREADS="5"

echo -e "\t copy BAM file..\n"

while  test  $# -gt 0;  do
	case "$1"  in
		-h|--help|-help)
			echo -e "\n\naccumulate_rpf_counts_per_putative_ORFs.sh\n"
			echo -e "This function will sum RPF read count per each putative ORF in the transcriptome, using reads aligned to the transcriptome as output by RSEM."
			echo -e "\n"

			echo -e "\nparameters:"
			echo -e "\t-runid r\trunid"
			echo -e "\t-offset f\toffset from 5' to determine read placement [default = $OFFSET]."
			echo -e "\n"

			echo -e "\t-o file\toutput file prefix.  writes to <out>.bam"

			echo -e "\t-p\tnumber of threads for multithreading [default: 1]."

			echo -e "\n\n"

			exit  0
			;;
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
		-o)
			# output prefix
			shift
			if test $# -gt 0; then
				OUTPREF="$1"
			else
				echo -e "\n\nERROR! output prefix not specified.\n\n"
				exit 1
			fi
			shift
			;;
		-p)
			# threads
			shift
			if test $# -gt 0; then
				NUMTHREADS="$1"
			else
				echo -e "\n\nERROR! number of threads not specified.\n\n"
				exit 1
			fi
			shift
			;;
		*)
			echo -e "\n\nERROR!  nothing requested!\n"
			exit 1
			;;
	esac
done



SPECIES="mm39"
GTFtag="Mus_musculus.GRCm39.104.gtf"
ORFDIR="/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/transcriptome-gtf/$SPECIES"
ORFTAG="orfs.$SPECIES.$GTFtag"
ORFPREF="$ORFDIR/$ORFTAG"

TMPWF="$RUNID"

if [ ! -z "$RUNID" ]; then
	echo -e "runid is given.  inferring all data...\n"

	$RPROF/src/retrieve_subset_from_tabular_database.sh  -r $RUNID  -c project,individual,runid,species,layout,strategy  -o $TMPWF.info.sub

	tail  -n +2  $TMPWF.info.sub  >  $TMPWF.info.sub.tmp	    \
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

echo -e "PROJECT=[$ORFPREF.orfs.bed.gz]\n"
# Run analysis
#gunzip  $ORFPREF.orfs.bed.gz


echo -e "sort orfs...\n"
# Sort is not required, as I already sorted BED file in R (Windows)
perl -pe 's/\r\n|\n|\r/\n/g' $ORFPREF.orfs.bed >$ORFTAG.orfs.bed
sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -o $ORFTAG.orfs.bed  $ORFTAG.orfs.bed



echo -e "merge orfs by class...\n"

awk  '{ $1 = $1 "====" $4 ; print $1 , $2 , $3  }'  FS="\t"  OFS="\t"  $ORFTAG.orfs.bed  >  $ORFTAG.orfs.classchr.bed


sort -t $'\t' -k1,1 -k2,2n --parallel=$NUMTHREADS -o $ORFTAG.orfs.classchr.bed  $ORFTAG.orfs.classchr.bed

#sort  -t $'\t'  -k1,1  -k2,2n  -k3,3n  -o $TMPWF.orfs.classchr.bed  $TMPWF.orfs.classchr.bed

bedtools  merge  -i $ORFTAG.orfs.classchr.bed  >  $ORFTAG.orfs.classchr.merge.bed

sed  -i  's/====/\t/'  $ORFTAG.orfs.classchr.merge.bed

awk  '{ print $1 , $3 , $4 , $2 }'  FS="\t"  OFS="\t"  $ORFTAG.orfs.classchr.merge.bed  >  $ORFTAG.orfs.classchr.merge.bed.tmp	\
&&  mv   $ORFTAG.orfs.classchr.merge.bed.tmp  $ORFTAG.orfs.classchr.merge.bed

sort --parallel=$NUMTHREADS -t $'\t' -k1,1 -k2,2n -k3,3n -k4,4 -o $ORFTAG.orfs.classchr.merge.bed $ORFTAG.orfs.classchr.merge.bed


############################################
echo -e "NUMTHREADS=[$NUMTHREADS]\n"
echo -e "Remote path of RSEMBAM=[$RSEMBAM.bam] \n"

#bsub -P hpcf_interactive -J hpcf_interactive -n 4 -q interactive -R "rusage[mem=5500]" -Is "bash"
############################################


echo -e "Local name of BAM TMPWF=[$TMPWF.bam]"
cp $RSEMBAM.bam $TMPWF.bam


echo -e "\t samtool index ..\n"
samtools index -@ 4 $TMPWF.bam


echo -e "\t samtool sort ..\n"
samtools sort -@ 4 -n -o $TMPWF.bam $TMPWF.bam


echo -e "\t samtool view ..\n"
samtools view -@ 4 -hb -o $TMPWF.rsem.map.bam $TMPOF.bam



echo -e "\t samtool view and bedtool bamtobed..\n"
samtools view -@ 4 $TMPWF.rsem.map.bam | grep -oP 'ZW:f:[^\t]+' | cut -d ':' -f 3 > $TMPWF.zw
bedtools bamtobed -i $TMPWF.rsem.map.bam| cut -f 1-3 | paste - $TMPWF.zw > $TMPWF.rsem.bam.bed
awk -v f=$OFFSET  '{ $3 = $2 + 1 + f ; print $0 }'  FS="\t"  OFS="\t"  $TMPWF.rsem.bam.bed  >  $TMPWF.rsem.bam.bed.pt


echo -e "sort bambed...\n"
sort -t $'\t' -k1,1 -k2,2n --parallel=4 -o $TMPWF.rsem.bam.bed.pt $TMPWF.rsem.bam.bed.pt


echo -e "bed map...\n"
echo -e "\tmap to each orf...\n"

bedtools map -c 4 -o sum -null 0 -a $ORFTAG.orfs.bed -b $TMPWF.rsem.bam.bed.pt > $TMPOF.bed


echo -e "\tcollapse by classification..\n"
bedtools map -c 4 -o sum -null 0 -a $ORFTAG.orfs.classchr.merge.bed -b $TMPWF.rsem.bam.bed.pt > $TMPOF.class.bed


time
