#!/usr/bin/env bash
sbatch \
  --mem-per-cpu=2000 \
  --time=20:00:00 \
  -o snake.out -e snake.err \
snakemake \
  -s SNVCalling.smk \
  --configfile config.yml \
  --profile profile_simple/ \
  --use-conda \
  --use-singularity \
  --singularity-args "-B `pwd` -B ${DATA_FOLDER}:/data -B /tmp:/tmp -B ${REF_FOLDER}:/ref:ro" \
  -pr \
  --latency-wait 30 \
  --show-failed-logs \
  --rerun-triggers mtime \
  "$@"
