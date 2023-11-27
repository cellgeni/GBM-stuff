#!/bin/bash

set -euo pipefail

script="./iterable-qc.R" #i.e. /lustre/scratch127/cellgen/cellgeni/simon/gbm_qc/iterable-qc.R
manifest="/path/to/manifest.tsv" #i.e. /nfs/cellgeni/tickets/tic-2507/actions/AT19-Manifest.tsv
indir="/path/to/donor/directory" #i.e. /nfs/team283/GBM_LEAP/phase_2_multiome_data/AT19
outdir="/path/to/output/directory" #i.e. /lustre/scratch127/cellgen/cellgeni/simon/gbm_qc/output

CPU=16
MEMORY=40000
QUEUE="normal"
GROUP="cellgeni"
IMAGE="/nfs/cellgeni/singularity/images/gbm_qc.sif"

bsub \
  -n "${CPU}" \
  -M "${MEMORY}" \
  -R"span [hosts=1] select[mem>${MEMORY}] rusage[mem=${MEMORY}]" \
  -o GBM-QC.%J.out \
  -e GBM-QC.%J.err \
  -q "${QUEUE}" \
  singularity exec -B /lustre,/nfs,/software "${IMAGE}" Rscript "${script}" --samples="${manifest}" --input="${indir}" --output="${outdir}"
