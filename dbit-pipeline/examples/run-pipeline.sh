#!/bin/bash

set -euo pipefail

pipeline="/nfs/cellgeni/tickets/tic-2598/actions/main.nf"

##PLEASE NOTE:
#dbit will produce 4 fastqs. 2 will contian long sequences of bases and the others will have incredibly short (8bp)
#this pipeline requires the 2 fastqs containing long sequences of bases 
#using: https://github.com/cellgeni/nf-irods-to-fastq to download fastqs will mean you need the R1 and R2 fastq
#i do not believe fastq order matters but when I submitted it was: "SAMPLE-ID\tR1-fastq\tR2-fastq"

samplefile="/path/to/samplefile" #i.e. ./pipeline.samples

nextflow run $pipeline \
  --samplefile $samplefile \
  -resume --ansi-log false
