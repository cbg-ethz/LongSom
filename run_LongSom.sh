#!/usr/bin/env bash
mkdir -p logs

bsub \
  -N \
  -R 'rusage[mem=2000]' \
  -W 24:00 \
  -oo logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out \
  -eo logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).err \
snakemake \
  -s snake/LongSom.snake.py \
  --cores 16 \
  --use-conda \
  --use-singularity \
  --singularity-args "-B `pwd` -B ${DATA_FOLDER}:/data -B /tmp:/tmp -B ${REF_FOLDER}:/ref:ro" \
  --configfile config/config.yaml \
  --profile ~/.config/snakemake/lsf/ \
  -pr \
  --latency-wait 30 \
  --rerun-incomplete \
  "$@"
