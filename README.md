# Unravelling intratumor heterogeneity in relapsed/refractory multiple myeloma (RRMM) using somatic variants in mitochondrial DNA

----
**Author:** Moritz Przybilla

----

**E-mail:** przybilla@stud.uni-heidelberg.de
**E-mail:** m.przybilla@dkfz-heidelberg.de

----

## Introduction (to be updated)

----
## Mitochondrial DNA (to be updated)

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

## Mitochondrial workflow

## Bulk Genotyping

In order to call variants within the mtDNA in bulk WGS, we are using a previously written script from *Na Cai*, who kindly shared it with us. The script is based on the availability of standard bam files containing all sequencing reads from WGS (both for nuclear and mitochondrial genome). The adapted version of Na's script can be found here:

```bash
# file path to bash script
/home/przybilm/bsub_commands/src/mtdna_preprocessing_variantcalling_MP.sh
```

For detailed information on the pre-processing and variant calling, please refer to the documentation within the script. In brief, extraction of MT reads from the original bam files is performed first (samtools), with subsequent remapping (picard) including splitting of paired end reads, map to an alternative mitochondrial reference genome (bwa) with trimming of 10 bases at the ends as well as sorting and indexing (samtools). For each step, an output path is specified, guiding the reader to the data which I produced.

### **Step 1**

To use the script, a conda environment with the respective tools was installed on the odcf-cluster. To implement all the following steps, the availability of a tab-delimited file with <sampleID> and <path to bam> is mandatory. The file was generated using an R script present within the following folder.

```bash
# file path to an example script for generating a samplesheet
~/mtdna_pipeline/R/get_sampleSheet.R
```

**Output**:

```bash
# the bamfile should be generated in your working directory
~/working_dir/bamlist.txt
```

This will look like this:

```bash
# example for a tab-delimited bamlist.txt file
K08K-1KAD5P_control1    /icgc/dkfzlsdf/project/hipo2/hipo_K08K/sequencing/whole_genome_sequencing/view-by-pid/K08K-1KAD5P/control1/paired/merged-alignment/control1_K08K-1KAD5P_merged.mdup.bam
K08K-1KAD5P_tumor1  /icgc/dkfzlsdf/project/hipo2/hipo_K08K/sequencing/whole_genome_sequencing/view-by-pid/K08K-1KAD5P/tumor1/paired/merged-alignment/tumor1_K08K-1KAD5P_merged.mdup.bam
K08K-1KAD5P_tumor2  /icgc/dkfzlsdf/project/hipo2/hipo_K08K/sequencing/whole_genome_sequencing/view-by-pid/K08K-1KAD5P/tumor2/paired/merged-alignment/tumor2_K08K-1KAD5P_merged.mdup.bam
K08K-1XZHPL_tumor1  /icgc/dkfzlsdf/project/hipo2/hipo_K08K/sequencing/whole_genome_sequencing/view-by-pid/K08K-1XZHPL/tumor1/paired/merged-alignment/tumor1_K08K-1XZHPL_merged.mdup.bam
```

### **Step 2**

Following this, **step 2** of Na's script can be executed. Within this step, all reads mapping to the mitochondrial genome will be extracted and stored as a per-sample bam file.

**Output folder**:

```bash
# example output folder
~/working_dir/rawMTbams/
```

The output folder will look like this:

```bash
# example output folder
-rw-r--r--. 1 przybilm B260  79619744 Aug 12  2019 K08K-1KAD5P_control1.mt.bam
-rw-r--r--. 1 przybilm B260      1224 Aug 12  2019 K08K-1KAD5P_control1.mt.bam.bai
-rw-r--r--. 1 przybilm B260 115155002 Aug 12  2019 K08K-1KAD5P_tumor1.mt.bam
-rw-r--r--. 1 przybilm B260      1352 Aug 12  2019 K08K-1KAD5P_tumor1.mt.bam.bai
-rw-r--r--. 1 przybilm B260 110552790 Aug 12  2019 K08K-1KAD5P_tumor2.mt.bam
-rw-r--r--. 1 przybilm B260      1352 Aug 12  2019 K08K-1KAD5P_tumor2.mt.bam.bai
```

### **Step 3**

In **step 3**, remapping and trimming of the mtDNA reads is performed and implemented.

At this point, the **alternative reference sequence** has to be used. All the respective reference files can be found following the following path.

```bash
# this is the folder where the reference sequence etc. is stored
~/mtdna_pipeline/REFs/
```

In contrast to the other steps, this one will create two output folders. One with the fastq files and a second one with sorted mitochondrial bam files.

**Output folder**:

```bash
# example output folder for the raw fastq files (MT reads only)
~/working_dir/MTfastq/

# example output folder for the sorted mt bam files
~/working_dir/MTremapped
```

The first output folder will look like this, containing distinct folders for each sample:

```bash
# example for MTfastq folder with an individual folder per sample
drwxr-sr-x. 2 przybilm B260 214 Aug 12  2019 K08K-1KAD5P_control1
drwxr-sr-x. 2 przybilm B260 204 Aug 12  2019 K08K-1KAD5P_tumor1
drwxr-sr-x. 2 przybilm B260 204 Aug 12  2019 K08K-1KAD5P_tumor2
```

While the second output folder will contain the mitochondrial bam files:

```bash
# example for MTremapped folder with a sorted bam file per sample
-rw-r--r--. 1 przybilm B260  75482061 Aug 12  2019 K08K-1KAD5P_control1.mt.bam
-rw-r--r--. 1 przybilm B260  72616735 Aug 12  2019 K08K-1KAD5P_control1.mt.sorted.bam
-rw-r--r--. 1 przybilm B260      1032 Aug 12  2019 K08K-1KAD5P_control1.mt.sorted.bam.bai
-rw-r--r--. 1 przybilm B260 110029564 Aug 12  2019 K08K-1KAD5P_tumor1.mt.bam
-rw-r--r--. 1 przybilm B260 105424750 Aug 12  2019 K08K-1KAD5P_tumor1.mt.sorted.bam
-rw-r--r--. 1 przybilm B260      1416 Aug 12  2019 K08K-1KAD5P_tumor1.mt.sorted.bam.bai
-rw-r--r--. 1 przybilm B260 105224513 Aug 12  2019 K08K-1KAD5P_tumor2.mt.bam
-rw-r--r--. 1 przybilm B260 101239248 Aug 12  2019 K08K-1KAD5P_tumor2.mt.sorted.bam
-rw-r--r--. 1 przybilm B260      1304 Aug 12  2019 K08K-1KAD5P_tumor2.mt.sorted.bam.bai
```

### **Step 4**

In **step 4**, filtering of the remapped bam files containing mtDNA reads only, is performed using ngsutils and samtools.

**Output folder**:

```bash
# example output folder for the filtered bam files
~/working_dir/filtered/
```

The output folder will look like this, containing distinct files for each sample:

```bash
# example for filtered folder with a filtered bam file per sample
-rw-r--r--. 1 przybilm B260  57627898 Aug 12  2019 K08K-1KAD5P_control1.filtered.bam
-rw-r--r--. 1 przybilm B260       152 Aug 12  2019 K08K-1KAD5P_control1.filtered.bam.bai
-rw-r--r--. 1 przybilm B260  70978143 Aug 12  2019 K08K-1KAD5P_control1.samfiltered.bam
-rw-r--r--. 1 przybilm B260       264 Aug 12  2019 K08K-1KAD5P_control1.samfiltered.bam.bai
-rw-r--r--. 1 przybilm B260  87235933 Aug 12  2019 K08K-1KAD5P_tumor1.filtered.bam
-rw-r--r--. 1 przybilm B260       152 Aug 12  2019 K08K-1KAD5P_tumor1.filtered.bam.bai
-rw-r--r--. 1 przybilm B260 103561583 Aug 12  2019 K08K-1KAD5P_tumor1.samfiltered.bam
-rw-r--r--. 1 przybilm B260       424 Aug 12  2019 K08K-1KAD5P_tumor1.samfiltered.bam.bai
-rw-r--r--. 1 przybilm B260  80028006 Aug 12  2019 K08K-1KAD5P_tumor2.filtered.bam
-rw-r--r--. 1 przybilm B260       152 Aug 12  2019 K08K-1KAD5P_tumor2.filtered.bam.bai
-rw-r--r--. 1 przybilm B260  99049175 Aug 12  2019 K08K-1KAD5P_tumor2.samfiltered.bam
-rw-r--r--. 1 przybilm B260       328 Aug 12  2019 K08K-1KAD5P_tumor2.samfiltered.bam.bai
```

### **Step 5**

The final pre-processing step before implementing the variant calling with GATK comprises the read group replacement. To do so, picard is used within the following script (**step 5**).

**Output folder**:

```bash
# example output folder for the filtered bam files with replaced read groups
~/working_dir/filtered_unique_rg/
```

```bash
# example for filtered_unique_rg folder with the final bam file per sample
-rw-r--r--. 1 przybilm B260       152 Aug 13  2019 K08K-1KAD5P_control1.bai
-rw-r--r--. 1 przybilm B260  59658246 Aug 13  2019 K08K-1KAD5P_control1.bam
-rw-r--r--. 1 przybilm B260       152 Aug 13  2019 K08K-1KAD5P_tumor1.bai
-rw-r--r--. 1 przybilm B260  90111766 Aug 13  2019 K08K-1KAD5P_tumor1.bam
-rw-r--r--. 1 przybilm B260       152 Aug 13  2019 K08K-1KAD5P_tumor2.bai
-rw-r--r--. 1 przybilm B260  82743741 Aug 13  2019 K08K-1KAD5P_tumor2.bam
```

After we implemented **step1** to **step5**, we are finally ready to perform variant calling of homoplasmic (over VAF ≥ 90%) and heteroplasmic variants (VAF < 90%) within the following steps of the pipeline.

### Homoplasmic variant calling (**Step 6**)

For calling homoplasmic variants, GATK is used. In this case, we leverage [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/article?id=11068) to determine the variants.

**Output folder**:

```bash
# example output folder for the homoplasmic variants
~/working_dir/gatk/gvcfs/
```

```bash
# example output folder for the homoplasmic variants containing an vcf file per sample
-rw-r--r--. 1 przybilm B260 5319 Aug 13  2019 K08K-1KAD5P_control1.g.vcf.gz
-rw-r--r--. 1 przybilm B260  124 Aug 13  2019 K08K-1KAD5P_control1.g.vcf.gz.tbi
-rw-r--r--. 1 przybilm B260 5407 Aug 13  2019 K08K-1KAD5P_tumor1.g.vcf.gz
-rw-r--r--. 1 przybilm B260  124 Aug 13  2019 K08K-1KAD5P_tumor1.g.vcf.gz.tbi
-rw-r--r--. 1 przybilm B260 5406 Aug 13  2019 K08K-1KAD5P_tumor2.g.vcf.gz
-rw-r--r--. 1 przybilm B260  123 Aug 13  2019 K08K-1KAD5P_tumor2.g.vcf.gz.tbi
```

After running the first part of the variant calling, we consolidate all variants present in the different samples. Therefore, we need to create a `genomedb.list` file, which is a tab delimited file with <sample ID> <full path to corresponding .g.vcf.gz>.

```bash
# file path to the make_genomedblist script
~/mtdna_pipeline/R/make_genomedbList.R
```

The output will be in .vcf.gz format which can be investigated using `less myfile.vcf.gz`. In particular, the resulting file will be called `allsamples.biallelicscnps.vcf.gz`. For further information about the vcf format, please have a look at the respective page from the [BroadInstitute](https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it).

### Heteroplasmic variant calling (**step 7**)

The most interesting step of the WGS workflow comprises the determination of heteroplasmic variants in the samples of interest. To do so, we leverage a publicly available tool, called [mtdnaserver](https://github.com/seppinho/mutserve), which is designed to exactly determine somatic mutations in the mitchondrial genome. As we know that the coverage in our multiple myeloma samples is very high, we run mtdnaserver with a low cuf-off of VAF = 0.1%.

**Output folder**:

```bash
# example output folder for the heteroplasmic variants
~/working_dir/mtdnaserver/
```

The output folder content will look like this:

```bash
# example output folder for the heteroplasmic variants containing an vcf, as well as txt files per sample
-rw-r--r--. 1 przybilm B260    3147 Aug 14  2019 K08K-1KAD5P_all.var.txt
-rw-r--r--. 1 przybilm B260    3914 Aug 13  2019 K08K-1KAD5P_control1.txt
-rw-r--r--. 1 przybilm B260    1088 Aug 13  2019 K08K-1KAD5P_control1.vcf.gz
-rw-r--r--. 1 przybilm B260 4539471 Aug 13  2019 K08K-1KAD5P_control1_raw.txt
-rw-r--r--. 1 przybilm B260    2501 Aug 13  2019 K08K-1KAD5P_tumor1.txt
-rw-r--r--. 1 przybilm B260     905 Aug 13  2019 K08K-1KAD5P_tumor1.vcf.gz
-rw-r--r--. 1 przybilm B260 4637307 Aug 13  2019 K08K-1KAD5P_tumor1_raw.txt
-rw-r--r--. 1 przybilm B260    2182 Aug 13  2019 K08K-1KAD5P_tumor2.txt
-rw-r--r--. 1 przybilm B260     869 Aug 13  2019 K08K-1KAD5P_tumor2.vcf.gz
-rw-r--r--. 1 przybilm B260 4288156 Aug 13  2019 K08K-1KAD5P_tumor2_raw.txt
```

Within the final step of the `mtdna_preprocessing_variantcalling.sh` script, we will again combine all variants together in a file. This will create an `all.var.txt` file, which will look like this:

```bash
ID	Pos	Ref	Variant	VariantLevel	MajorBase	MajorLevel	MinorBase	MinorLevel	Coverage	Type
K08K-1KAD5P_control1	2706	A	G	1.00	G	0.00	-	0.00	8802	1
K08K-1KAD5P_control1	1438	A	G	1.00	G	0.00	-	0.00	9158	1
K08K-1KAD5P_control1	720	T	C	0.009	T	0.99	C	0.009	6786	2
K08K-1KAD5P_control1	750	A	G	0.999	G	0.00	-	0.00	6801	1
K08K-1KAD5P_control1	72	T	A	0.999	A	0.00	-	0.00	1190	1
K08K-1KAD5P_control1	2539	A	C	0.001	A	0.999	C	0.001	8599	2
K08K-1KAD5P_control1	1245	T	C	0.002	T	0.998	C	0.002	8891	2
K08K-1KAD5P_control1	3777	T	C	0.002	T	0.998	C	0.002	7953	2
K08K-1KAD5P_control1	16189	T	C	0.001	T	0.998	C	0.001	6166	2
```

This is the end of the `mtdna_preprocessing_variantcalling.sh` script. From this point, we will need to prepare the output of the heteroplasmic variant calling for the single-cell genotyping. In order to perform this step, it is crucial to perform variant filtering. Importantly, we do not want to investigate uninformative germline variants.

### Filtering heteroplasmic variant calls

After calling the heteroplasmic (somatic) variants, we will filter the outputs according to the VAF in control and tumor samples. Therefore, a custom R script was used. It should be noted, that this step is very much dependent on the manual examination of your variant call results. In theory, you can write your own scripts to filter these outputs and customize it to your needs.

```bash
# file path to the filtering script with a normal reference
~/mtdna_pipeline/R/mtdna_var_filtering.R
```

To get a dataset to work with, we assume that variants with a VAF ≤ 5% in the control could still be wrongly detected, while similarly it could still be present in the tumor with VAF ≥ 0.5%. This is really much a way to increase the number fo variants, potentially being valuable for lineage tracing. This is neither very sophisticated, nor precise. However, it is very important anyway that you investigate your data, and eventually also check within IGV, whether your mutation calls are valid.

```bash
# example output folder for the filtered heteroplasmic variants
~/working_dir/mtdnaserver/filtered/
```

The output from this script will be a patient-specific file with all the distinct variants from individual samples. Note, that the filtering thresholds will also be indicated in the file name (control_VAF ≤ 5%; tumor_VAF ≥ 0.5%).

```bash
# example output folder for the filtered heteroplasmic variants with a variant txt file per sample
-rw-r--r--. 1 przybilm B260    346 Aug 14  2019 1KAD5P_0.05_0.005_filtered.var.txt
-rw-r--r--. 1 przybilm B260    598 Aug 14  2019 37HWC4_0.05_0.005_filtered.var.txt
-rw-r--r--. 1 przybilm B260    636 Aug 14  2019 41R2SE_0.05_0.005_filtered.var.txt
```

These files will have the following format (tab-delimited):

```bash
Mutation    K08K.1KAD5P_control1_VariantLevel   K08K.1KAD5P_control1_Coverage   K08K.1KAD5P_tumor1_VariantLevel K08K.1KAD5P_tumor1_Coverage K08K. 1KAD5P_tumor2_VariantLevel    K08K.1KAD5P_tumor2_Coverage
10946:A:C       0.005   3870    0.006   6567    0.005   4551
12705:C:T       0.005   9587    0.003   15404   0.005   13280
2701:G:A        0       0       0.993   12992   0.706   11921
8280:A:C        0       0       0       0       0.007   1209
```

In addition, the script mentioned above will create a `all_patient_0.05_0.005_filtered.var.txt`, that will contain

### Convert txt files to vcf files

In the very last step for the bulk sequencing analysis, we need to convert the `txt` files to `vcf` files, as they are the input for our single-cell genotyping in the second stage of this workflow. Thus, we provide a script which takes the `filtered.var.txt` generated above as input and converts them into `vcf` files. In addition, it also filters for common sequencing error resulting from changes in the mitochondrial reference genome. You can find more about these changes [here](https://www.mitomap.org/foswiki/bin/view/MITOMAP/MitoSeqs). The respective positions are found in the `REFs/chrM_REF_mutations.vcf` file and will be excluded from the final `vcf` file.

```bash
# file path to the filtering script with a normal reference
~/mtdna_pipeline/R/var_txt_to_vcf.R
```

Running this script will generate `filtered.var.vcf` files for each patient (samples get merged together). The output folder will comprise one file per patient.

```bash
# example output folder for the filtered heteroplasmic variants with a variant txt file per sample
-rw-r--r--. 1 przybilm B260  61 Feb  4 15:12 K08K-1KAD5P_0.05_0.01_filtered.var.vcf
-rw-r--r--. 1 przybilm B260 883 Feb  4 15:12 K08K-1XZHPL_0.05_0.01_filtered.var.vcf
-rw-r--r--. 1 przybilm B260 413 Feb  4 15:12 K08K-21Z2S9_0.05_0.01_filtered.var.vcf
```

Checking the files themselves, you will see this:

```bash
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chrM	2701	.	G	A	.	.	.
```

Ultimately, this is the result which will be used in the second stage of the workflow, where we will infer these *ground truth* variants in single cells.

----

## Single-cell Genotyping (to be updated)
