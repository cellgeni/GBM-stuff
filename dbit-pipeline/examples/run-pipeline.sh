#!/bin/bash

set -euo pipefail

nextflow run /nfs/cellgeni/tickets/tic-2598/actions/main.nf \
  --samplefile "/nfs/cellgeni/tickets/tic-2598/actions/pipeline.samples" \
  -resume --ansi-log false
