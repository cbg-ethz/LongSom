library(infercnv)
library(Seurat)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/mnt/featurecount/B486_Tum.counts.formated.txt",
                                    annotations_file="/mnt/ctypes/B486_Tum.txt",
                                    delim="\t",
                                    gene_order_file="/mnt/gene_ordering_file.txt",
                                    ref_group_names=c("T.NK.cells", "Fibroblasts", "Endothelial.cells","Mesothelial.cells","Myeloid.cells","B.cells"),
                                    min_max_counts_per_cell=c(1e3,1e7))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="uncorrected_dir_0.1",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T
                             )
