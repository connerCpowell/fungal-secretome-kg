#!/bin/bash
set -e

FASTA=data/raw/y_secreted.fasta
OUTDIR=data/processed/

mkdir -p $OUTDIR

signalp \
  --fastafile $FASTA \
  --organism eukarya \
  --output_dir $OUTDIR

