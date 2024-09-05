#!/usr/bin/env bash
mkdir -p logs
sbatch \
  --mem-per-cpu=2000 \
  --time=20:00:00 \
  -o logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out \
  -e logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).err \
snakemake \
  -s snakefiles/QC.smk \
  --configfile config/config_OvCa_LR.yml \
  --profile profile_simple/ \
  --use-conda \
  --use-singularity \
  --singularity-args "-B `pwd` -B ${DATA_FOLDER}:/data -B /tmp:/tmp -B ${REF_FOLDER}:/ref:ro" \
  -pr \
  --latency-wait 30 \
  --show-failed-logs \
  "$@" #--rerun-triggers mtime \
