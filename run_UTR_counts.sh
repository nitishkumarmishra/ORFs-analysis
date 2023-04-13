#!/bin/bash

for sample in $(ls -d out/161021/* |cut -d "/" -f 3)
do
        echo "****************** Sample $sample submitted ******************"
        bsub -q priority -n 4 -P UTRs -R "rusage[mem=20000]" -J UTR_count -oo %J.out -eo %J.err "sh UTR_counts.sh -runid $sample"
        echo -e "#################### finish for sample $sample #####################""\n"

done

