# Unravelling intratumor heterogeneity in relapsed/refractory multiple myeloma (RRMM) using somatic variants in mitochondrial DNA

----
***Author:*** Moritz Przybilla

----

***E-mail:*** przybilla@stud.uni-heidelberg.de

----

## Introduction (from ASH Abstract 2019)

Multiple myeloma (MM) is a heterogeneous malignancy of clonal plasma cells that accumulate in the bone marrow (BM). Despite new treatment approaches, in most patients treatment-resistant subclones are selected by therapy, resulting in the development of refractory disease. While the subclonal architecture in newly diagnosed patients has been investigated in great detail, intra-tumor heterogeneity in relapsed/refractory (RR) MM is poorly characterized. Recent technological advances provide the opportunity to analyze tumor samples at single-cell (sc) level with high accuracy and througput. Here, we present a pilot study for an integrative analysis of sc Assay for Transposase-Accessible Chromatin with high-throughput sequencing (scATAC-seq) and scRNA-seq with the aim to comprehensivly study the regulatory landscape, gene expression, and evolution of individual subclones in RRMM patients.

----
## Mitochondrial DNA

The foundation for the following approach was layed by two different studies, which investigated the possibility of using somatic mutation in the mitochondrial genome to infer lineages in cancer-unrelated and related context. *Ludwig et al.* (Cell, 2019) and *Xu et al.* (eLife, 2019) thereby leveraged a combination of bulk and single-cell (sc) sequencing in order to determine these variants. Recent studies have shown that mitochondrial DNA (mtDNA) has improtant properties, which allow the usage as a "tracking device" of the cancer evolution.

**Some fun facts about mtDNA:**

* Circular, intron-free, double-stranded genome that is only 16.5 kb in size, encompassing 13 protein coding genes OXPHOS – oxidative phosphorylation function
* 10x higher mutation rate compared to nuclear genome
* 100s – 1000s of mitochondria per cell
* Most mtDNA are likely to be passenger mutations in cancer
  * No association with known mutagens such as cigarette smoke or UV light (Signature 4 & 7)

Most of these facts were derived from a big mtDNA genomics study in the context of cancer from 2014 (*Ju et al.*, eLife).

----
## Idea

Ludwig et al., as well as others, made use of the fact, that the mitochondrial genome is captured within conventional assay for transposase accessible chromatin (ATAC)- and RNA-sequencing, thus delivering this information ‘for free’ 4,7. Importantly, each mitochondrion in a cell holds its own, circular, intron-free mtDNA, comprising only 16.5 kb in size, with the number of mitochondria per cell sometimes exceeding 10,000. Due to that, the coverage with which the mtDNA is captured, is sufficiently higher than the target sequencing depth of nuclear DNA 4,6,8. While the increased coverage in theory allows highly confident mutation calls, Ludwig et al.’s approach involved the joint bulk ATAC-seq and single cell RNA-seq (scRNA-seq) profiling of the mitochondrial genome. In general, bulk ATAC-seq was used first in order to determine a ‘pool’ of somatic mutations, representing the ‘ground truth’. Subsequently, scRNA-seq was used to infer different lineages carrying mutations present in this ‘pool’ of variants.

----
## Whole-genome sequencing (WGS)

The results from the study indicated above provide evidence that WGS provides a much deeper coverage of the mitochondrial genome than enriched sequencing methods like WES (*Ju et al.*, 2014, eLife). In theory, the coverage of the nuclear genome translates to the coverage of the mtDNA by approximately 100 times more (e.g. 30x coverage refers to approximately 3000x mtDNA coverage). In contrast to that, Ju et al. reported a low average coverage of mtDNA (92x) for a number of 971 samples which underwent WES (*Ju et al.*, 2014, eLife). Thus, we will use WGS to determine the somatic variants present in the *pool*, from which we can then infer lineages.

### Workflow

In order to call variants within the mtDNA, we are using a previously developed script from *Na Cai*, who kindly shared it with us. The script is based on the availability of standard bam files containing all sequencing reads from WGS. The adapted version of Na's script can be found here:

    /home/przybilm/bsub_commands/src/mtdna_preprocessing_variantcalling_MP.sh

For detailed information on the pre-processing and variant calling, please refer to the documentation within the script. In brief, extraction of MT reads from the original bam files is performed first (samtools), with subsequent remapping (picard) including splitting of paired end reads, map to an alternative mitochondrial reference genome (bwa) with trimming of 10 bases at the ends as well as sorting and indexing (samtools). For each step, an output path is specified, guiding the reader to the data which I produced.

##### ***Step 1***

To use the scripts a conda environment with the respective tools was installed on the odcf-cluster. To implement all steps of the pipeline, the availability of a tab-delimited file with <sampleID> and <path to bam> is mandatory. The file was generated using an R script present within the following folder.

    /home/przybilm/bsub_commands/WGS/step1/get_sampleSheet.R

***Output***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/bamlist.txt

##### ***Step 2***

Following this, ***step2*** of Na's script can be executed. The resulting bash script was submitted to the odcf-cluster using the following script.

    /home/przybilm/bsub_commands/WGS/step2/mtdna_bams.sh

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/rawMTbams/

***add information about the coverage plots here***

##### ***Step 3***

In ***step3***, remapping and trimming of the mtDNA reads is performed and implemented as described before. The respective script can be found using this path.

    /home/przybilm/bsub_commands/WGS/step3/picard.sh

At this point, the **alternative reference sequence** has to be used. All the respective files can be found following the following path.

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/REFs/

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/MTfastq/

##### ***Step 4***

In ***step4***, filtering of the remapped bam files containing mtDNA reads only, is performed using ngsutils and samtools.

    /home/przybilm/bsub_commands/WGS/step4/mtDNA_filter.sh

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/filtered/

##### ***Step 5***

The final pre-processing step before implementing the variant calling with GATK comprises the read group replacement. Therefore, picard is used within the following script (***step5***).

    /home/przybilm/bsub_commands/WGS/step5/ReplaceReads.sh

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/filtered_unique_rg/

After performance of ***step1*** - ***step5***, we are ready to perform variant calling of homoplasmic (over 90% VAF) and heteroplasmic variants (below 90% VAF).


***add information about the coverage plots here***


#### Homoplasmic variant calling (**step6**)

For calling homoplasmic variants, GATK is used. In this case, we leverage [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/article?id=11068) to determine the variants.

    /home/przybilm/bsub_commands/WGS/step6/HaplotypeCaller.sh

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/gatk/gvcfs/

After running the first part of the variant calling, we consolidate all variants present in the different samples. Therefore, we need to create a genomedb.list file, which is a tab delimited file with <sample ID> <full path to corresponding .g.vcf.gz>. The result is an allsamples.biallelicscnps.vcf.gz file.

    /home/przybilm/bsub_commands/WGS/step6/make_genomedbList.R

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/gatk/

The output will be in .vcf.gz format which can be investigated using ``` $ less myfile.vcf.gz```. For further information about the vcf format, please have a look at the respective page from the [BroadInstitute](https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it).

#### Heteroplasmic variant calling (**step7**)

##### mtdnaserver

The most interesting step of the WGS workflow comprises the determination of heteroplasmic variants in the samples of interest. To do so, we leverage a publicly available tool, called [mtdnaserver](https://github.com/seppinho/mutserve), which is designed to exactly determine somatic mutations in the mitchondrial genome. Although we know that the coverage was very high, we run mtdnaserver with different cutoffs, thus being able to capture most of the mutations. In the following scripts we use a conventional (1%), low (0.1%) and zero (0%) as argument.


    /home/przybilm/bsub_commands/WGS/step7/mtdnaserver/variantcalling.sh
    /home/przybilm/bsub_commands/WGS/step7/mtdnaserver/variantcalling_lowVAF.sh
    /home/przybilm/bsub_commands/WGS/step7/mtdnaserver/variantcalling_zeroVAF.sh

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/mtdnaserver/

After calling the somatic variants, the raw output was filtered according to the VAF in control and tumor samples. Therefore, a custom R script was used.

   /home/przybilm/bsub_commands/WGS/mtdnaserver/filtering_mtdnaserver_control-tumor-filter.R

To get a dataset to work with, 5% FDR were applied to mutations called for the tumor, while similarly 1% FDR were applied to the controls of the output from the low VAF variant calling. The respective dataframe was subsetted into patient-specific mutation files.

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/mtdnaserver/filtered/lowVAF/

These files were used as reference for cellSNP mode1.


***add information about the plots here***


##### mutect2

As an alternative to mtdnaserver, [mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php) was used for somatic variant calling. Therefore, a list of matched control and tumor samples was produced and mutect2 was run submitting the following script.

    /home/przybilm/bsub_commands/WGS/step7/mutect2/MuTect2.sh

***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/MuTect2/

The output of mutect2 is again in vcf file format.

##### Mpileup

In order to check the coverage of the called mutations with an additional method, *samtools mpileup* was used.

    /home/przybilm/bsub_commands/WGS/step7/run_mpileup.sh


***Output folder***:

    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/mpileup


----
## Single-cell Assay for Transposase Accessible Chromating sequencing (sc-ATAC-seq)

In the previously published literature, scATAC-seq is indicated as the method to use for the investigation of mtDNA. Thus, we first assessed mutations in this dataset. The tool EMBLEM, which was published by *Xu et al.* was made publicly available on github ([EMBLEM](https://github.com/ChangLab/EMBLEM)), which is why we preferentially use this.

### EMBLEM

Epigenome and Mitochondrial Barcode of Lineage from Endogeneous Mutations (EMBLEM) relies on joint bulk and scATAC-seq, which effectively enriches for mtDNA in comparison to WGS (17x higher coverage). The approach which the authors describe in their methods section is cited below:

**Quote from the method section of** ***Xu et al., 2019, eLife*** **:**

> **Single cell ATAC-seq data processing and mitochondrial DNA variant calling**
> Single cell ATAC-seq were processed similarly to the bulk ATAC-seq, taking each individual cell as one sample. After cleanin the alignment, files from every single cell were merged and heteroplasmic variants were first called with the merged bam and filtered using the same criteria as bulk data. Heteroplasmic variants called from merged data or from bulk data were re-counted in each individual cell using ``` samtools -q 20 -Q 20 ```. The non-reference allele had to match the variants detected in merged or bulk data.


> **Detection rate estimation**
> In every single cell, if the variant allele detected in merged or bulk data were supported by any reads, it was considered positive; otherwise it was counted as zero. A binary matrix was used to present the lineage relationship among single cells and plotted as a heatmap.

In general, they are using this method to infer intermediate populations in the clonal evolution of primary hematopoietic stem cells (pHSCs). They could show that a population of pHSC carry a clone with pHSC-specific mutations as well as another clone with shared mutations with leukemic stem cells (LSCs).

Consequently, we will try to use the same pipeline they published on Github, starting with their [scATAC-seq pipeline](https://github.com/ChangLab/ATAC_mito_sc) for alignment and mapping as well as filtering of fastq to bam files, which shall then serve as an input for the EMBLEM pipeline to determine somatic variants in the mtDNA.


### cellSNP

Since the EMBLEM pipeline is currently not available due to missing files on Github (correspondence with Jing Xu already started), we have to consider alternative ways to call the somatic variants. Therefore, we first use Yuanhua's [cellSNP](https://github.com/huangyh09/cellSNP) which is actually made for scRNA-sequencing data, but works for scDNA data as well.

We use the output from the cellranger v1.1 scATAC pipeline, which is essentially a position sorted bamfile.

From these bamfiles, the mtDNA reads are extracted and saved into sample specific bams. Subsequently, those are used for the


    /icgc/dkfzlsdf/analysis/B260/projects/przybilm/cellSNP/results/cellSNP/

    /icgc/dkfzlsdf/analysis/hipo2/hipo_K08K/cellRanger-ATAC_v1.1/GRCh38/


### monovar

[monovar](https://bitbucket.org/hamimzafar/monovar/src/master/)


## Single-cell RNA sequencing (sc-RNA-seq)

Cellranger output
    /icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3_GRCh38
    /icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/cellranger_results_v3_GRCh19

Script to generate cellranger output
    /icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/scripts/cellranger_pipeline/job.arr_cellranger.sh

Input file for the cellranger bash script

    /icgc/dkfzlsdf/analysis/hipo2/hipo_K43R/info_files/HIPO_IDs_Hana_change_final.txt

Read in seurat count matrix by

    SAMPLE/outs/filtered_feature_bc_matrix/


### cellSNP

[cellSNP](https://github.com/huangyh09/cellSNP)

