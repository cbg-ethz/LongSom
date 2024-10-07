#!/usr/bin/env bash
OUTPUT_DIR=/path/to/output_dir
REF_DIR=path/to/LongSom/ref/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir

snakemake \
  -s workflow/LongSom.smk \
  --configfile config/config.yaml \
  --use-conda \
  --use-singularity \
  --singularity-args "-B ${OUTPUT_DIR}:/output -B /tmp:/tmp -B ${REF_DIR}:/ref" \
  --show-failed-logs \
  --rerun-incomplete \
  -p \


