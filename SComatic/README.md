# SComatic
SComatic is a tool that provides functionalities to detect somatic single-nucleotide mutations in high-throughput single-cell genomics and transcriptomics data sets, such as single-cell RNA-seq and single-cell ATAC-seq.

If you use SComatic (see **License** at the bottom of this page), please cite our publication [Muyas et al. 2023](https://www.nature.com/articles/s41587-023-01863-z).

For further details on SComatic, its assumptions, limitations and applications, please see [Muyas et al. 2023](https://www.nature.com/articles/s41587-023-01863-z).

![Algorithm](/docs/Algorithm.jpeg)

*Created with BioRender.com*

## Installation and requirements
SComatic requires Python version >=3.7.0, R version >=3.6.0, [samtools](https://github.com/samtools/samtools), [bedtools](https://bedtools.readthedocs.io/en/latest/) and datamash (>=v1.1.0, only for building your own panel of normals). 
- We strongly recommend to build your own [conda](https://docs.conda.io/en/latest/) environment as follows:
```
conda create -n SComatic -c bioconda python=3.7 r-base=3.6.1 samtools datamash bedtools
conda activate SComatic
```
- Additional dependencies can be installed by running the following commands:
```python
pip install -r requirements.txt
```
```R
Rscript r_requirements_install.R
```
If your R version is >= 3.6.0 but < 4, you might have issues installing the `VGAM` package. If this is the case, try this:
```R
Rscript r_requirements_install.v3_6.R
```

- Unpack PoN (Panel of Normals) files:
```
gunzip PoNs/PoN.scRNAseq.hg38.tsv.gz
gunzip PoNs/PoN.scATACseq.hg38.tsv.gz 
```
- Unpack RNA editing file:
```
gunzip RNAediting/AllEditingSites.hg38.txt.gz
```
If you use the RNA editing database above, please cite the following articles: [Tan et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/29022589), [Kiran et al. 2012](https://www.ncbi.nlm.nih.gov/pubmed/23074185) and [Picardi et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/27587585).

# Detection of somatic mutations in single-cell data sets using SComatic
We show below how to run SComatic for the detection of somatic mutations in scRNA-seq data. SComatic requires two data types as input:

- Aligned sequencing reads in BAM format for all cell types analysed. **The input BAM file must contain the cell type barcode information in the cell barcode tag “CB” (as reported by popular tools, such as Cell Ranger, 10x Genomics), and ideally, the "nM" tag (number of mismatches) and the "NH" tag (number of hits)**.
- A file listing precomputed cell type annotations where each row reports the cell barcode and cell type for each cell analysed. **Importanly, do not include doublets in this file**.

SComatic consists of the following 4 steps, each of which is run using a different Python script as indicated below.

## Step 1: Splitting alignment file into cell-type-specific bams
The first step consists of splitting the BAM file containing aligned sequencing reads for all cell types in a sample  into cell-type-specific BAM files using precomputed cell type annotations. **The BAM file must contain the cell type barcode information in the cell barcode tag “CB” (as reported by popular tools, such as Cell Ranger, 10x Genomics)**. In addition, we strongly suggest generating the input bam files to contain the **"nM" tag** (number of mismatches) and the **"NH" tag** (number of hits), which will be used to remove low-quality reads for downstream analysis (*--max_nM* and *--max_NH* parameters. However, we included the possibility of running SComatic without these last two tags (nM and NH), at the user's own risk. 

Step 1 is executed using the script SplitBam/SplitBamCellTypes.py, which has the following parameters:

- List of parameters:
```
python scripts/SplitBam/SplitBamCellTypes.py --help
usage: SplitBamCellTypes.py [-h] --bam BAM --meta META [--id ID]
                            [--max_nM MAX_NM] [--max_NH MAX_NH]
                            [--min_MQ MIN_MQ] [--n_trim N_TRIM]
                            [--outdir OUTDIR]

Split alignment file into cell type specific BAMs

optional arguments:
  -h, --help       show this help message and exit
  --bam BAM        BAM file to be analysed (Sorted by coordinate)
  --meta META      Metadata file mapping cell barcodes to cell type
                   information
  --id ID          Sample ID
  --max_nM MAX_NM  Maximum number of mismatches permitted to consider reads
                   for analysis. By default, this filter is switched off,
                   although we recommed using --max_nM 5. If applied, this
                   filter requires having the nM tag in the bam file.
                   [Default: Switched off]
  --max_NH MAX_NH  Maximum number of alignment hits permitted to consider
                   reads for analysis. By default, this filter is switched
                   off, although we recommend using --max_NH 1. This filter
                   requires having the NH tag in the bam file. [Default:
                   Switched off]
  --min_MQ MIN_MQ  Minimum mapping quality required to consider reads for
                   analysis. Set this value to 0 to switch this filter off.
                   --min_MQ 255 is recommended for RNA data, and --min_MQ 30
                   for DNA data. [Default: 255]
  --n_trim N_TRIM  Number of bases trimmed by setting the base quality to 0 at
                   the beginning and end of each read [Default: 0]
  --outdir OUTDIR  Out directory
```

The precomputed cell type annotation file provided with the --meta parameter must contain at least the following two columns (Index for cell barcode ID and Cell_type for the precomputed cell type annotation) and must be a tab-separated file. Cell type annotations containing whitespaces or any of the following special characters (~ . ` ! @ # $ % ^ & * ( ) { | } / \ : ; " ' < > ? , = +) are not supported. Dashes and underscores are supported. Whitespace characters in the filenames are not supported. **Importanly, do not include doublets in this file or cell types with unkown cell type annotations.**

```
Index Cell_type
AAACCTGCATGCTAGT  Epithelial
AAACCTGGTAGCCTAT  Epithelial
AAACCTGGTTGTCGCG  Epithelial
AAACCTGTCATGTGGT  Epithelial
AAACCTGTCCTTGGTC  Epithelial
AAACCTGTCGGATGTT  T_cell
AAACCTGTCGTACGGC  T_cell
AAACCTGTCTTGCAAG  T_cell
AAACGGGAGACGCACA  T_cell
```

In addition to the cell-type specific bam files, this script creates a txt (\*.report.txt) file showing the number of reads filter, and importantly, why there were filtered out. 

- **Example:** check [here](/docs/SComaticExample.md) to see how to run this step with an example sample.  

## Step 2: Collecting base count information
Base count information for each cell type and for every position in the genome is recorded in a base count matrix indexed by cell types and genomic coordinates. 

The command line to run this step is: 

- List of parameters:
```
python scripts/BaseCellCounter/BaseCellCounter.py --help
usage: BaseCellCounter.py [-h] --bam BAM --ref REF --chrom CHROM
                                   [--out_folder OUT_FOLDER] [--id ID]
                                   [--nprocs NPROCS] [--bin BIN] [--bed BED]
                                   [--bed_out BED_OUT] [--min_ac MIN_AC]
                                   [--min_af MIN_AF] [--min_dp MIN_DP]
                                   [--min_cc MIN_CC] [--min_bq MIN_BQ]
                                   [--min_mq MIN_MQ] [--tmp_dir TMP_DIR]

Script to obtain a list of base and cell counts in scRNA bam file

optional arguments:
  -h, --help            Show this help message and exit
  --bam BAM             Tumor bam file to be analysed
  --ref REF             Reference genome. *fai must be available in the same
                        folder as reference
  --chrom CHROM         Chromosome to be analysed. Use --chrom all to analyse
                        all chromosomes
  --out_folder OUT_FOLDER
                        Output folder
  --id ID               Prefix of out file. If provided, please use the following
                        format: *.[cell_type] . Example: sample1.t_cell. If
                        not provided, it is taken from the bam file
  --nprocs NPROCS       Number of processes [Default = 1]
  --bin BIN             Bin size for running the analysis [Default 50000]
  --bed BED             Regions to focus the analysis on. Three-column bed file
  --bed_out BED_OUT     Regions to ignore. Three-column bed
                        file
  --min_ac MIN_AC       Minimum alternative count to consider a genomic
                        site for further analysis. Default = 0
  --min_af MIN_AF       Minimum alternative allele fraction to consider a
                        genomic site for further analysis. Default = 0
  --min_dp MIN_DP       Minimum coverage to consider a genomic site for further analysis. Default
                        = 5
  --min_cc MIN_CC       Minimum number of cells required to consider a genomic
                        site for further analysis. Default = 5
  --min_bq MIN_BQ       Minimum base quality permited for the base counts.
                        Default = 20
  --min_mq MIN_MQ       Minimum mapping quality required to analyse read.
                        Default = 255
  --tmp_dir TMP_DIR     Temporary folder for tmp files
```

- **Example:** check [here](/docs/SComaticExample.md) to see how to run this step with an example sample.  

## Step 3: Merging  base count matrices
In Step 3, SComatic takes as input base count matrices computed in Step 2 for all cell types analysed to merge them into a single base count matrix, which is stored in tsv format. Individual base count matrices to be merged need to be stored in the same directory.

- List of parameters:
```python 
python scripts/MergeCounts/MergeBaseCellCounts.py --help
usage: MergeBaseCellCounts.py [-h] --tsv_folder TSV_FOLDER --outfile
                                       OUTFILE

Script to merge the cell/base counts tsv files per cell type in only one file

optional arguments:
  -h, --help            Show this help message and exit
  --tsv_folder TSV_FOLDER
                        Folder with cell/base count files (tsv) per cell type.
                        Avoid not desired tsv files in this folder
  --outfile OUTFILE   Out file name
```

- **Example:** check [here](/docs/SComaticExample.md) to see how to run this step with an example sample.  

## Step 4: Detection of somatic mutations
The last step consists of running two scripts to call somatic mutations.

### Step 4.1
SComatic applies a set of hard filters and Beta binomial tests to discount sites affected by recurrent technical artefacts as somatic mutations. 

- List of parameters:
```python
python scripts/BaseCellCalling/BaseCellCalling.step1.py --help
usage: BaseCellCalling.step1.py [-h] --infile INFILE --outfile
                                         OUTFILE --ref REF [--min_cov MIN_COV]
                                         [--min_cells MIN_CELLS]
                                         [--min_ac_cells MIN_AC_CELLS]
                                         [--min_ac_reads MIN_AC_READS]
                                         [--max_cell_types MAX_CELL_TYPES]
                                         [--min_cell_types MIN_CELL_TYPES]
                                         [--fisher_cutoff FISHER_CUTOFF]
                                         [--alpha1 ALPHA1] [--beta1 BETA1]
                                         [--alpha2 ALPHA2] [--beta2 BETA2]

Script to perform the scRNA somatic variant calling

optional arguments:
  -h, --help            Show this help message and exit
  --infile INFILE       Input file with all samples merged in a single tsv
  --outfile OUTFILE     Output file prefix
  --ref REF             Reference fasta file (*fai must exist)
  --min_cov MIN_COV     Minimum depth of coverage to consider a sample.
                        [Default: 5]
  --min_cells MIN_CELLS
                        Minimum number of cells with sequencing depth at a site to consider a
                        genomic site for further analysis. [Default: 5]
  --min_ac_cells MIN_AC_CELLS
                        Minimum number of cells supporting the alternative
                        allele to call a mutation. [Default: 2]
  --min_ac_reads MIN_AC_READS
                        Minimum number of reads supporting the alternative
                        allele to call a mutation. [Default: 3]
  --max_cell_types MAX_CELL_TYPES
                        Maximum number of cell types carrying a mutation to
                        make a somatic call. [Default: 1]
  --min_cell_types MIN_CELL_TYPES
                        Minimum number of cell types with enough coverage across enough
                        cells to consider a site as callable [Default: 2]
  --fisher_cutoff FISHER_CUTOFF
                        P value cutoff for the Fisher's exact test performed to
                        detect strand bias. A float value is expected, if applied,
                        we recommend 0.001. By default, this test is switched
                        off with a value of 1 [Default: 1]
  --alpha1 ALPHA1       Alpha parameter for Beta distribution of read counts.
                        [Default: 0.260288007167716]
  --beta1 BETA1         Beta parameter for Beta distribution of read counts.
                        [Default: 173.94711910763732]
  --alpha2 ALPHA2       Alpha parameter for Beta distribution of cell counts.
                        [Default: 0.08354121346569514]
  --beta2 BETA2         Beta parameter for Beta distribution of cell counts.
                        [Default: 103.47683488327257]
```

To estimate new Beta binomial parameters whenever required by the user, SComatic provides the following [scripts](/docs/betabinomialestimation.md).

- **Example:** check [here](/docs/SComaticExample.md) to see how to run this step with an example sample.  

### Step 4.2 
Scomatic takes the output of the previous step (4.1) and applies additional filters based on external datasets (RNA editing and Panel of Normals), and flags clustered mutations. High quality mutations are marked with the label “PASS” in the FILTER column of the output file. When using the provided [RNA editing](RNAediting) file, please cite the following articles: [Tan et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/29022589), [Kiran et al. 2012](https://www.ncbi.nlm.nih.gov/pubmed/23074185) and [Picardi et al. 2017](https://www.ncbi.nlm.nih.gov/pubmed/27587585).

- List of parameters: 
```python 
python scripts/BaseCellCalling/BaseCellCalling.step2.py --help
usage: BaseCellCalling.step2.py [-h] --infile INFILE --outfile
                                         OUTFILE [--editing EDITING]
                                         [--pon PON]
                                         [--min_distance MIN_DISTANCE]

Script to perform the scRNA somatic variant calling

optional arguments:
  -h, --help            Show this help message and exit
  --infile INFILE       Input file with all samples merged in a single tsv
  --outfile OUTFILE     Output file prefix
  --editing EDITING     RNA editing file to be used to remove RNA-diting sites
  --pon PON             Panel of normals (PoN) file to be used to remove
                        germline polymorphisms and recurrent artefacts
  --min_distance MIN_DISTANCE
                        Minimum distance allowed between potential somatic
                        variants [Default: 5]
```

* The *--pon* parameter permits to work with different types of formats/files depending on the availability and quantity of non-neoplastic samples. These are the main options: 

     * Using the [PoNs](/PoNs/) provided in this repository, which were computed using the data described in the main manuscript of SComatic.
     * Using your custom PoNs, which can be built using [SComatic modules](/docs/pon.md).
     * TSV file listing the mutations detected by running SComatic (output of *Detection of somatic mutations: Step 4.2*) on scRNA-seq data from e.g., matched non-neoplastic cells.
     * Custom PoN in VCF format (unzipped) generated by running a variant caller, such as *GATK-HaplotypeCaller*, on DNA sequencing (WES or WGS) data, such as a matched normal sample or a set of unmatched germline samples.
     * Using a PoN in VCF format (unzipped) generated by other tools like [*GATK*](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-)

- **Example:** check [here](/docs/SComaticExample.md) to see how to run this step with an example sample.  

## [Estimating new beta-binomial parameters](/docs/betabinomialestimation.md)

SComatic models the background error rate of the technology used to generate the single-cell data (e.g., single-cell RNA-seq) using a Beta binomial distribution. Specifically, non-reference allele counts at homozygous reference sites are modelled using a binomial distribution with parameter P (error rate), which is a random variable that follows a Beta distribution with parameters α and β. 
Default values for the Beta binomial tests used in  Step 4.1 are computed using the data sets described in the manuscript. However, we provide [scripts](/docs/betabinomialestimation.md) to allow the user to reparameterize the Beta binomial using other data sets. 

## [Generating a custom Panel of Normals](/docs/pon.md)

In Step 4.2, SComatic uses a Panel of Normals (PoN) to detect systematic errors and germline contamination in the somatic mutation callset. The PoN provided in this repository is computed using the data described in the manuscript (Hg38 reference genome). However, SComatic provides a [script](/docs/pon.md) to build a custom PoN using other data sets if required. 

## [Other SComatic functionalities](/docs/OtherFunctionalities.md)
SComatic provides the following additional functionalities, which are described in detail [here](/docs/OtherFunctionalities.md).

- Computing the number of callable sites per cell type
- Computing the number of callable sites per cell
- Computing the genotype for each cell at the variant sites
- Computing the trinucleotide context background
- Computing germline genotypes for known variants in single-cell datasets

## [FAQs - Frequently asked questions](/docs/faqs.md)
This section answers some of the users' most recurrent doubts when running SComatic.

1. [Are the SComatic parameters for scATAC-seq data the same as for scRNA-seq data?](/docs/faqs.md#1-are-the-scomatic-parameters-for-scatac-seq-data-the-same-as-for-scrna-seq-data)
2. [How can we perform the variant annotation with the SComatic output?](/docs/faqs.md#2-how-can-we-perform-the-variant-annotation-with-the-scomatic-output)
3. [Can SComatic work with other types of PoN files?](/docs/faqs.md#3-can-scomatic-work-with-other-types-of-pon-files)
4. [Can we use the calls from other callers to genotype unique cells using SComatic?](/docs/faqs.md#4-can-we-use-the-calls-from-other-callers-to-genotype-unique-cells-using-scomatic)
5. [How do different cell type labels (e.g. different levels of granularity) affect the SComatic performance?](docs/faqs.md#5-how-do-different-cell-type-labels-eg-different-levels-of-granularity-affect-the-scomatic-performance)
6. [What does it happen if CellRanger does not properly trim all non-genomic sequences (adapters) from the reads?](/docs/faqs.md#6-what-does-it-happen-if-cellranger-does-not-properly-trim-all-non-genomic-sequences-adapters-from-the-reads)
7. [How to interpret the SingleCellGenotype.py output?](https://github.com/cortes-ciriano-lab/SComatic/blob/main/docs/faqs.md#7-how-to-interpret-the-singlecellgenotypepy-output)

## Contact
If you have any comments or suggestions about SComatic please raise an issue or contact us: 

Francesc Muyas: fmuyas@ebi.ac.uk 

Isidro Cortes-Ciriano: icortes@ebi.ac.uk

## License
**SComatic is free for academic use only**. If you are not a member of a public funded academic and/or education and/or research institution you must obtain a commercial license from EMBL Enterprise Management GmbH (EMBLEM); please email EMBLEM (info@embl-em.de).
