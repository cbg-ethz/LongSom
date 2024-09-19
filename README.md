# LongSom
A [Snakemake](https://snakemake.github.io/) pipeline for calling somatic SNVs, fusions and CNAs in PacBio long-read single-cell RNA-seq cancer samples, using the [Trinity Cancer Transcriptome Analysis Toolkit (CTAT)](https://github.com/NCIP/Trinity_CTAT), and infer clones based on them.

LongSom takes a bam file and a barcodes file as input, and then uses [ctat-mutations](https://github.com/NCIP/ctat-mutations) to call SNVs, [ctat-LR-fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion) to call fusions. It then uses Bayesian non-parametric clustering [BnpC](https://github.com/cbg-ethz/BnpC) to cluster cells into subclones based on called SNVs and fusions. In parallel, LongSom uses [inferCNV](https://github.com/broadinstitute/infercnv) to call CNAs and cluster cells into subclones based on them.

## Contents
- [Installation](##Installation)
- [Usage](##Usage)
- [Cite](##Cite)

## Requirements
- Python 3.X
- Mamba/Conda 
- Singularity (https://sylabs.io/docs/)

## Installation

### Clone repository
First, download LongSom from github and change to the directory:
```bash
git clone https://github.com/cbg-ethz/LongSom
cd LongSom
```
### Create a conda environment with Snakemake
Install Snakemake: 
```bash
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n LongSom snakemake snakemake-executor-plugin-slurm
```
Using Mamba is highly recommended, for more information. visit [Snakemake's installation guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Then, activate the environment:
```bash
conda activate LongSom
```
This environment should be activated each time you want to use LongSom

### Install Subread (featurecount)
You can download [Subread](https://sourceforge.net/projects/subread/) and intall it this way:  
```bash
wget https://sourceforge.net/projects/subread/files/subread-2.0.6/subread-2.0.6-source.tar.gz
tar zxvf subread-2.0.6-source.tar.gz
cd subread-2.0.6-source/src/
make -f Makefile.Linux
```

### Install CTAT softwares
Download the simg of those three tools:
* ctat-LR-fusion (https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/ctat-LR-fusion/) (tested on V0.13.0)
* ctat InferCNV (https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/InferCNV/) (tested on V1.16.0)
* ctat-mutations (https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/CTAT_MUTATIONS/) (tested on V4.0.0)

Place all simg in the `bin` folder

### Install BnpC
Follow [BnpC installation instructions](https://github.com/cbg-ethz/BnpC?tab=readme-ov-file#Installation) (create a conda environment called BnpC). 

## Usage

File requirements:
* .bam file (with BC as barcode tag)
* barcodes .txt file 
* genome .fa file (hg38)
* transcriptome .gtf file  (https://www.gencodegenes.org/human/)

Before each usage, you should source the LongSom environment:

```bash
conda activate LongSom
```

The LongSom wrapper script `run_LongSom.py` can be run with the following shell command:
```bash
./run_LongSom
```

It should run for less than a day on HPC. Output files should be found in the `results` folder.


### Before running the pipeline


* **config file**
  * input directory
    Before running the pipeline, the `config/config.yaml` file needs to be adapted to contain the path to input bam files. It is provided in the first section (`specific`) of the config file.
  * resource information
    In addition to the input path, further resource information must be provided in the section `specific`. This information is primarily specifying
     the genomic reference used for the reads mapping and the transcriptomic reference required for isoform classification. An example `config.yaml` file ready for adaptation, as
    well as a brief description of the relevant config blocks, is provided in the directory `config/`.

* **reference files**
  * A genome fasta file (http://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38)
  * A GENCODE gene annotation gtf file (https://www.gencodegenes.org/human/)

* **sample map**
  * Provide a sample map file, i.e. a tab delimited text file listing all samples that should be analysed, and how many bam files are associated to it (see example below). ID will be used to name files and identify the sample throughout the pipeline.
  * Sample map example:
  ```
  sample     files
  SampleA     2
  SampleB     4
  SampleC     2
  ```
* **input data**
  * This pipeline take as input either concatenated or unconcatenated reads PacBio CCS bam files. I you use concatenated reads input, files should be named `SampleA_1.bam`, `SampleA_2.bam`, `SampleB_1.bam`, etc. (sample name should correspond to the sample map).  If you use unconcatenated reads as input, files should be named `SampleA_1.subreads.bam`, etc.

## Cite
 Arthur Dondi, Nico Borgsmüller, Pedro Ferreira, Brian Haas, Francis Jacob, Viola Heinzelmann-Schwarz, Tumor Profiler Consortium, Niko Beerenwinkel. De novo detection of somatic variants in long-read single-cell RNA sequencing data. Available on biorxiv soon