#!/bin/bash
module load parallel
module load samtools
bsub -q priority -n 14 -P Samtools -R "rusage[mem=25000]" -J STAR_BAM -oo Bam_2_SAM.out -eo Bam_2_SAM.err "parallel --plus 'samtools view -h {} -o {...}.sam' ::: *.bam"
