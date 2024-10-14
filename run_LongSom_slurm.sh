#!/usr/bin/env bash
OUTPUT_DIR=/path/to/output_dir
REF_DIR=REF_DIR=path/to/LongSom/ref/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir

mkdir -p logs
mkdir -p ${OUTPUT_DIR}
sbatch \
  --time=20:00:00 \
  -o logs/snakelog.out \
  -e logs/snakelog.err \
snakemake \
  -s workflow/Snakefile \
  --configfile config/config.yaml \
  --profile profile/ \
  --use-conda \
  --use-singularity \
  --singularity-args "-B ${OUTPUT_DIR}:/output -B /tmp:/tmp -B ${REF_DIR}:/ref" \
  --show-failed-logs \
  --rerun-incomplete \
  -p \
 "$@"

