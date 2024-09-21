# FAQs - Frequently asked questions

## 1. Are the SComatic parameters for scATAC-seq data the same as for scRNA-seq data?
No, they are different. ATAC-seq data is a DNA-based approach. Therefore, there are differences in the way of processing the bam files. These are the parameters that you should change in comparison to the scRNA-seq workflow:
- **Step 1**: remove the _--max_nM_ and _--max_NH_ parameters, set the mapping quality filter to _--min_MQ_ 30 .
- **Step 2**: set the mapping quality filter to _--min_mq_ 30
- **Step 4.2**: Remove the _--editing_ parameter, and _--pon_ altered to point at the scATACseq PoN provided in this GitHub repo (or custom one)

## 2. How can we perform the variant annotation with the SComatic output?
In our manuscript, we have annotated all our variants using [annovar](https://annovar.openbioinformatics.org/en/latest/). With a couple of command lines, we can adapt our variant calling output for this annotation tool. Using our example data, we should run these lines:

**1.  Preparing SComatic output for annovar**

```
ANNOVAR=/path/to/annovar
hummandb=/path/to/annovar_humandb_hg38

grep -v '#' Example.calling.step2.pass.tsv |  tr '\t' '-' | awk -F'-' -v OFS='\t' '{print $1,$2,$3,$4,$5,$0}' > sample.variants.avinput
```

**2. Annotate variants using the annovar-ready variant file**

```
perl $ANNOVAR/table_annovar.pl sample.variants.avinput \
                $hummandb -buildver hg38 \
                -out sample.variants.annovar \
                -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f \
                -nastring . -csvout -polish --otherinfo
```

This step will generate a comma-separated file (.csv) called _sample.variants.annovar.hg38_multianno.csv_ , which can be easily opened and processed in your more desired software (R, excel...). The output should look like this:

```
Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,cytoBand,ExAC_ALL,ExAC_AFR,ExAC_AMR,ExAC_EAS,ExAC_FIN,ExAC_NFE,ExAC_OTH,ExAC_SAS,avsnp147,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,CADD_raw,CADD_phred,DANN_score,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,integrated_fitCons_score,integrated_confidence_value,GERP++_RS,phyloP7way_vertebrate,phyloP20way_mammalian,phastCons7way_vertebrate,phastCons20way_mammalian,SiPhy_29way_logOdds,gnomAD_genome_ALL,gnomAD_genome_AFR,gnomAD_genome_AMR,gnomAD_genome_ASJ,gnomAD_genome_EAS,gnomAD_genome_FIN,gnomAD_genome_NFE,gnomAD_genome_OTH,Otherinfo
chr10,29559501,29559501,A,T,"intronic","SVIL",.,.,.,"10p11.23",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-29559501-29559501-A-T-PASS-Epithelial_cells-ATTCA-TTGAT-1-12-5-4-2-0.3333-0.4-0.0-0.0001-2-2-0;38;1-0;16;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-12|5|3:0:2:0:0:0|8:0:4:0:0:0|306:0:164:0:0:0|0:0:0:0:0:0|8:0:4:0:0:0-30|13|13:0:0:0:0:0|30:0:0:0:0:0|1189:0:0:0:0:0|0:0:0:0:0:0|30:0:0:0:0:0"
chr10,73413109,73413109,G,A,"intronic","ANXA7",.,.,.,"10q22.2",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-73413109-73413109-G-A-PASS-Epithelial_cells-TGGAG-CATGG-1-12-8-4-3-0.3333-0.375-0.0-0.0-2-2-0;15;1-0;11;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-12|8|3:0:0:5:0:0|4:0:0:8:0:0|160:0:0:302:0:0|0:0:0:0:0:0|4:0:0:8:0:0-7|6|0:0:0:6:0:0|0:0:0:7:0:0|0:0:0:270:0:0|0:0:0:0:0:0|0:0:0:7:0:0"
chr10,89413978,89413978,G,A,"upstream","IFIT5","dist=590",.,.,"10q23.31",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-89413978-89413978-G-A-PASS-Epithelial_cells-GTATG-GTGTT-1-18-11-7-4-0.3889-0.3636-0.0-0.0-2-2-0;20;1-0;12;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-18|11|4:0:0:7:0:0|7:0:0:11:0:0|270:0:0:417:0:0|0:0:0:0:0:0|7:0:0:11:0:0-9|5|0:0:0:5:0:0|0:0:0:9:0:0|0:0:0:339:0:0|0:0:0:0:0:0|0:0:0:9:0:0"
chr10,97332670,97332670,G,A,"UTR3","FRAT2","NM_012083:c.*1201C>T",.,.,"10q24.1",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-97332670-97332670-G-A-PASS-Epithelial_cells-TTGGG-GTGGC-1-12-11-6-6-0.5-0.5455-0.0-0.0-2-2-0;12;1-0;11;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-6|6|0:0:0:6:0:0|0:0:0:6:0:0|0:0:0:237:0:0|0:0:0:0:0:0|0:0:0:6:0:0-12|11|6:0:0:5:0:0|6:0:0:6:0:0|229:0:0:225:0:0|0:0:0:0:0:0|6:0:0:6:0:0-NA"
chr10,126926205,126926205,G,A,"intronic","DOCK1",.,.,.,"10q26.2",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-126926205-126926205-G-A-PASS-Epithelial_cells-GATCC-GTGAC-1-35-17-15-6-0.4286-0.3529-0.0-0.0-2-2-0;35;1-0;18;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-35|17|6:0:0:12:0:0|15:0:0:20:0:0|603:0:0:786:0:0|15:0:0:20:0:0|0:0:0:0:0:0-15|7|0:0:0:7:0:0|0:0:0:15:0:0|0:0:0:590:0:0|0:0:0:15:0:0|0:0:0:0:0:0"
```

As you will notice, the first few columns of the file have the annotation information provided by _annovar_ (Gene, region, impact, Gnomad allele frequencies...). Importantly, the final column of the csv file (_Otherinfo_) shows all the info found in the original _Example.calling.step2.pass.tsv_, but separated by the symbol "-" . The column names of each one of these _Otherinfo_ items are the same as the ones found in the _Example.calling.step2.pass.tsv_ : 

```
> grep '^#CHROM' Example.calling.step2.pass.tsv

#CHROM	Start	End	REF	ALT	FILTER	Cell_types	Up_context	Down_context	N_ALT	Dp	Nc	Bc	Cc	VAF	CCF	BCp	CCp	Cell_types_min_BC	Cell_types_min_CC	Rest_BC	Rest_CC	Fisher_p	Cell_type_Filter	INFO	Myeloids	Epithelial_cells	Stromal_cells

```

Of course, depending on the columns you want to keep for the annotation process, you might slightly change the command line in the _step 1_ described above, as well as the parameters provided in the annovar computation. The parameters provided here are the ones used in our manuscript. 

## 3. Can SComatic work with other types of PoN files?
Yes, SComatic can deal with PoNs obtained by other tools.

The *--pon* parameter permits working with different types of formats/files depending on the availability and quantity of non-neoplastic samples. These are the main options: 

* Using the [PoNs](/PoNs/) provided in this repository, which were computed using the data described in the main manuscript of SComatic.
* Using your custom PoNs, which can be built using [SComatic modules](/docs/pon.md).
* TSV file listing the mutations detected by running SComatic (output of *Detection of somatic mutations: Step 4.2*) on scRNA-seq data from e.g., matched non-neoplastic cells.
* Custom PoN in VCF format (unzipped) generated by running a variant caller, such as *GATK-HaplotypeCaller*, on DNA sequencing (WES or WGS) data, such as a matched normal sample or a set of unmatched germline samples.
* Using a PoN in VCF format (unzipped) generated by other tools like [*GATK*](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-)

## 4. Can we use the calls from other callers to genotype unique cells using SComatic?
Yes, SComatic allows us to use calls obtained by other callers (in vcf format) to check the presence of each mutation at single-cell resolution. However, a minor edit is required to transform the vcf format into the SComatic format. Here is the required command line to be run in the vcf:

```
# Prepare SComatic input from an external vcf for the script SingleCellGenotype.py
VCF=/path/to/sample.vcf
grep -v '^#' $VCF | awk -F'\t' -v OFS='\t' '{print $1,$2,$2,$4,$5,".","Cell_type",".",".",".",".",".",".","."}' > Mutations.other_caller.tsv

SCOMATIC=/path/to/SComatic
cell_type_bam=/path/to/cell_type.bam
META=/path/to/cell_type_annotations.tsv
REF=/path/to/reference_genome.fasta
python $SCOMATIC/scripts/SingleCellGenotype/SingleCellGenotype.py --bam $cell_type_bam  \
    --infile Mutations.other_caller.tsv \
    --nprocs 16  \
    --meta $META   \
    --outfile Mutations.other_caller.single_cell_genotype.tsv  \
    --tmp_dir temp  \
    --ref $REF
```

## 5. How do different cell type labels (e.g. different levels of granularity) affect the SComatic performance?
SComatic can be run using different levels of granularity in terms of cell type annotations. 

As illustrated in the schematic below, in the case of non-neoplastic samples the relatedness between cell types from a development perspective determines which types of mutations can be detected. Using very granular cell type annotations that consider e.g. two cell types from the same differentiation hierarchy as different cell types, will remove any mutation present in the progenitor common to both cell types (assuming that sufficient sampling of all cell types in the data set analysed is achieved). As a result, only mutations acquired after clonal diversification will be detected (marked in green in the schematic). In contrast, mutations acquired in progenitor/stem cells and during early development will be discounted based on the fact that they will be present in all descendant cells (marked in red in the schematic). It follows from the preceding that mutations acquired very early during development will likely be considered germline based on their presence in multiple cell types when calling mutations using diverse tissue types from the same individual. By contrast, using broader cell type annotations (e.g., epithelial, immune, etc) permits the detection of mutations accumulated over longer periods of time in e.g., adult stem cells up to the point of lineage diversification. Overall, determining the granularity of the cell type annotations depends on the biological question of interest. We note that SComatic can be run using cell type annotations of variable granularity easily, which should enhance the applicability of our algorithm to diverse research areas.

<div align="center">
<img src="/docs/Granularity_plot.png" width="60%">
</div>


For the analysis of cancer samples, we note that the effect of cell type granularity is less relevant as compared to the analysis of non-neoplastic samples, in that somatic mutations are only expected to be detected in the cancer cells. Thus, unless cell state annotations are used for the cancer cells, mutations present in cancer cells would be readily detected when running SComatic if one cell type annotation is used for all cancer cells and non-neoplastic cell types are included to remove germline polymorphisms.

## 6. What does it happen if CellRanger does not properly trim all non-genomic sequences (adapters) from the reads?
This might be an issue and could increase the number of false positive calls in the final call set. If this is the case, you can use the parameters _--n_trim  NUMER_OF_BASES_TO_IGNORE_ in _Step1_ (_SplitBamCellTypes.py_) to ignore the first and last bases of each read. In addition, when running the [_CellRanger count_](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct) tool, provide the _–-chemistry_ parameter as precisely as possible. It will help to remove these non-desired regions. 

## 7. How to interpret the SingleCellGenotype.py output? 

This script provides base count and allele information for all cells with reads covering each variant site. By default, all cells (mutated and non-mutated)  will be present in the output file (see `--alt_flag` parameter to change it). However, it is important to understand the meaning of each one of the columns of this output file to properly interpret these results. Here a short description:

| Column | Description |
| --- | --- |
| #CHROM |  Chromosome carrying the mutation |
| Start |  Start genomic coordinate |
| End |  End genomic coordinate |
| REF |  Reference allele |
| ALT_expected |  Alternative allele as described in the input file (`--infile`) |
| Cell_type_expected |  Cell types harbouring the mutation as described in the input file (`--infile`) |
| Num_cells_expected |  Number of expected cells carrying the mutation as described in the input file (`--infile`) |
| CB | Unique cell barcode analysed |
| Cell_type_observed | Cell type attributed to the analysed _CB_ according to the input metadata file (`--meta`) |
| Base_observed | Allele observed in this _CB_ |
| Num_reads | Number of reads carrying the _Base_observed_ | 

<br>

Let's understand this table with an example. Looking at our [SComatic example data](https://github.com/cortes-ciriano-lab/SComatic/blob/main/docs/SComaticExample.md), we will focus on the variant site _chr10-29559501_ and the _SComatic/example_data/results/SingleCellAlleles/Epithelial_cells.single_cell_genotype.tsv_ file generated.

```
#CHROM	Start	End	REF	ALT_expected	Cell_type_expected	Num_cells_expected	CB	Cell_type_observed	Base_observed	Num_reads
chr10	29559501	29559501	A	T	Epithelial_cells	2	AGTCTTTGTGCATCTA	Epithelial_cells	A	5
chr10	29559501	29559501	A	T	Epithelial_cells	2	CCCTCCTAGGCTAGGT	Epithelial_cells	A	1
chr10	29559501	29559501	A	T	Epithelial_cells	2	GGGTCTGTCTTGAGGT	Epithelial_cells	T	2
chr10	29559501	29559501	A	T	Epithelial_cells	2	GTCCTCAAGGCTCATT	Epithelial_cells	T	2
chr10	29559501	29559501	A	T	Epithelial_cells	2	GAGTCCGAGGGTGTTG	Epithelial_cells	A	2
```

The columns `ALT_expected`, `Cell_type_expected`	and `Num_cells_expected` correspond to the values observed in the `--infile Example.calling.step2.pass.tsv`, so they represent the calls at cell type resolution. 

In contrast, the columns `CB`, 	`Cell_type_observed`, `Base_observed` and `Num_reads` correspond to the allele observations at unique cell resolution when interrogating the bam files.  

Each CB can be presented in the output file in as many rows as different alleles are found per cell, although in most cases, we only observed one allele per cell (so one row per unique CB). In order to find the alleles harbouring the called mutation, we have to look for those rows (unique CBs) where `ALT_expeced == Base_observed` and  `Cell_type_expected == Cell_type_observed`. In general terms, CBs not accomplishing these conditions can be understood as noise or non-mutated cells. 

Although this script is designed to be run with the SComatic calls, it can also be run with an [external vcf file](https://github.com/cortes-ciriano-lab/SComatic/edit/main/docs/faqs.md#4-can-we-use-the-calls-from-other-callers-to-genotype-unique-cells-using-scomatic). As in most of the vcf-based cases we do not know the `Cell_type_expected`, we cannot check the `Cell_type_expected == Cell_type_observed` condition. However, we still can check the `ALT_expeced == Base_observed` condition. 
