#!/usr/bin/env Rscript

library(infercnv)
library(Seurat)

args = commandArgs(trailingOnly=TRUE)
raw_counts_matrix = args[1]
annotations_file = args[2]
gene_order_file = args[3]
out_dir = args[4]

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=raw_counts_matrix,
                                    annotations_file=annotations_file,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names=c("NonCancer"),
                                    min_max_counts_per_cell=c(1e3,1e7))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir=out_dir,  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
                             )
