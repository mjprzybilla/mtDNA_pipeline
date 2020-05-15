# *********************************************************************
#
# Calling variants with the mtDNA reads from bam files containing sequencing reads
#
# Date: 1 July 2019
#
# Author: Na Cai, Moritz Przybilla
#
# Summary: This is a guide for pre-processing and calling variants from BAM files
#          containing sequencing reads for mtDNA. Pre-processing involves sectioning
#          out reads from the mtDNA from the BAMs, checking their alignment and read
#          groups to ensure a single read group per donor/cell (whichever is the unit),
#          realigning if necessary, pruning out read pairs that do not both map to the
#          mtDNA with high quality, pruning out reads with excess heteroplasmies (likely
#          contaminations). Variant calling will be performed with both GATK HaplotypeCaller
#          and local version of mtDNA-server. Afterwards, we will call variants from BAM files
#          containing sequencing reads for mtDNA from scATAC-sequencing. Pre-processing involves 
#          sectioning out reads from the mtDNA from the BAMs and Variant calling from 
#          single cell data using Yuanhua's cellSNP.
#
# Dependencies: To implement this script, a conda environment with the installed libraries
#               of samtools, picard, bwa, gatk4, ngsutils and locally installed
#               mtdnaserver, as well as cellSNp is required.
#
# *********************************************************************

## *********************************************************************
## PRE-PROCESSING
## *********************************************************************

## bash mtdna_preprocessing_variantcalling.sh /icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS /icgc/dkfzlsdf/analysis/B260/projects/przybilm/stanford/WGS/output chrM
#!/bin/bash

## load conda environment
source ~/.bashrc
conda activate snakes

## enable wildcard expansions in script
shopt -s extglob
shopt | grep extglob

## step1: prepare the mtDNA reference - currently the mitochondrial reference genome used is NC_012920, also known as "rCRS"
## to read more about the mtDNA reference, please see: https://www.mitomap.org/foswiki/bin/view/MITOMAP/MitoSeqs
## all files I use for mtDNA reference can be found in mtfasta.tar
## in mtfasta.tar, there's the fasta file hs37d5.mt.fa and the following indices that will be needed
## 1. hs37d5.mt.fa: samtools faidx hs37d5.mt.fa (samtools see: https://github.com/samtools/samtools)
## 2. hs37d5.mt.fa.sa, hs37d5.mt.fa.pac, hs37d5.mt.fa.bwt, hs37d5.mt.fa.ann, hs37d5.mt.fa.amb: bwa index hs37d5.mt.fa (bwa aligner see: http://bio-bwa.sourceforge.net/bwa.shtml)

## step2: extract mt reads from bams using samtools (see https://github.com/samtools/samtools)
## you may have to try MT, M or chrM, depending on which reference the bam was mapped to and what the labelling of mtDNA was
## can inspect by doing "samtools view $file MT | head" to see if there's anything, same for the other labellings
## specify an working directory that comprises the bam files which comprise the mitochondrial reads as well as the 
## standard reads. Additionally, give an output directory where everything should be saved in 
## Hipo 067
## wdir="/icgc/dkfzlsdf/project/hipo/hipo_067/sequencing/whole_genome_sequencing/view-by-pid"
## odir="/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_067"
## Hipo_K08K
## wdir="/icgc/dkfzlsdf/project/hipo2/hipo_K08K/sequencing/whole_genome_sequencing/view-by-pid"
## odir="/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K"
## bamdir="/icgc/dkfzlsdf/analysis/hipo2/hipo_K08K/CellRanger-ATAC_v1.1/GRCh38/" - scATAC data
## Medulloblastoma
## wdir="/icgc/dkfzlsdf/analysis/B060/share/B060_Stuttgart_case/bam_files/"
## odir="/icgc/dkfzlsdf/analysis/B260/projects/przybilm/medulloblastoma/"
## bamdir="/icgc/dkfzlsdf/analysis/B260/projects/chromothripsis_medulloblastoma/data/10xCNV" - for single-cell bam files from scRNA
## chrM="chrM"
wdir=$1
odir=$2
mkdir $odir
chrM=$3
echo "working directory: $wdir"
echo "output directory: $odir"
echo "reference ID: $chrM"

## generate the list of bams in a tab delimited file with <sample ID> <path to bam>
echo "Generating sample sheet..."
Rscript /home/przybilm/bsub_command/WGS/mtdna_pipeline/get_sampleSheet.R $wdir $odir
wait
echo "Done"

# make folders in output directory to refine output folder structure
mkdir $odir/rawMTbams 2>/dev/null
mkdir $odir/rawMTbams/src 2>/dev/null
mkdir $odir/rawMTbams/log 2>/dev/null
rm $odir/rawMTbams/src/getMT.s* 2>/dev/null
## make use of the sample sheet generated before
samplelist=$(cat $odir/bamlist.txt | cut -f 1)  # take sample id
for sample in $samplelist
do echo $sample # std out sample name
if [ ! -f $odir/rawMTbams/$sample.mt.bam ]
then file=$(grep $sample $odir/bamlist.txt | cut -f 2) # use the path to the respective bam file for the sample
y=${file%.bam} ; ## this gets the whole path before .bam
z=${y##*/} ; ## this removes the path leaving only the bam name
# calculate percentage of mitochondrial reads in WGS 
# echo -e "total_counts=$(samtools view -c $file) ;
# MT_counts=$(samtools view -c $file $chrM) ; \
# echo "Total number of read counts in the bam file: $total_counts" ; \
# echo "Total number of mitochondrial read counts in the bam file: $MT_counts" ; \
# echo "Percentage of mitochondrial reads in sample:" ; \
# echo "scale=3 ; $MT_counts / $total_counts" | bc ; \
echo -e "samtools view -b $file $chrM > $odir/rawMTbams/$sample.mt.bam ; \
samtools index $odir/rawMTbams/$sample.mt.bam" >> $odir/rawMTbams/src/getMT.sh
else echo "done"
fi
done

# submit the bash script to the cluster
echo "Extracting all reads mapped to the mitochondrial genome..."
split -l 10 $odir/rawMTbams/src/getMT.sh  $odir/rawMTbams/src/getMT.split 2>/dev/null
for file in $odir/rawMTbams/src/getMT.split*
do chmod u+rwx $file
z=${file##*/} ;
bsub -R "select[mem>=10000] rusage[mem=10000]" -M10000 -q verylong -o $odir/rawMTbams/log/$z.o  \
-e $odir/rawMTbams/log/$z.e $file
done
wait
echo "Done"

## step3: prepare the bams for use (optional/only if needed)
## there may be a need to remap them if they're on a difference reference
## use picardtools to remap by first converting mapped bams to fastqs, then remapping to the NC_012920 reference
## picardtools see https://broadinstitute.github.io/picard/
fromdir="$odir/rawMTbams"
mkdir $odir/MTremapped 2>/dev/null
mkdir $odir/MTfastq 2>/dev/null
refdir="/icgc/dkfzlsdf/analysis/B260/projects/przybilm/REFs"
bamdir="$odir/MTremapped"
todir="$odir/MTfastq"
mkdir $todir/src 2>/dev/null
mkdir $todir/log 2>/dev/null
rm $todir/src/outputFastqs.* 2>/dev/null
for file in $fromdir/*.mt.bam
do echo $file
y=${file%.mt.bam} ;
echo "$y" ;
z=${y##*/} ;
echo "$z" ;
mkdir $todir/${z} 2>/dev/null
## A. the picardtools command allows you to split paired end reads into .1.fq containing the first part of each pair, .2.fq that contains the second part of each pair, and .u which contains reads that are not paired
## validation stringency is by default STRICT, but many at times this won't allow you to run anything, so putting to SILENT is fine
## B. then the bwa aln comman allows you to remap all the mtDNA reads to the NC_012920 reference
## we use -q10 -t4 as options. -q10 trims the reads for the 10 bases at the ends of the reads - this is because sequencing error is usually much higher on the ends.
## -t4 means using 4 threads, can change this according to how multithreaded your processes can be
## C. bwa sampe takes the two mapped paired ends and combines them into a single bam file.
## D. samtools sort sorts the remapped bam file by coordinates, samtools index indexes the bam file allowing random access
## !!! please pay attention to the identifier of the mitochondrial genome - MT, chrM - choose the reference accordingly (two versions in the folder) !!!
echo -e "picard SamToFastq I=$file FASTQ=$todir/${z}/${z}.1.fq \
SECOND_END_FASTQ=$todir/${z}/${z}.2.fq UNPAIRED_FASTQ=$todir/${z}/${z}.u VALIDATION_STRINGENCY=SILENT; \
bwa aln -q10 -t4 $refdir/hs37d5.mt.fa $todir/${z}/${z}.1.fq > $todir/${z}/${z}.1.sai ;\
bwa aln -q10 -t4 $refdir/hs37d5.mt.fa $todir/${z}/${z}.2.fq > $todir/${z}/${z}.2.sai ;\
bwa sampe $refdir/hs37d5.mt.fa $todir/${z}/${z}.1.sai $todir/${z}/${z}.2.sai $todir/${z}/${z}.1.fq $todir/${z}/${z}.2.fq | \
samtools view -Sb - > $bamdir/${z}.mt.bam ;\
samtools sort $bamdir/${z}.mt.bam | samtools view -Sb > $bamdir/${z}.mt.sorted.bam; \
samtools index $bamdir/${z}.mt.sorted.bam" >> $todir/src/outputFastqs.sh
# fi
done

# submit the bash script to the cluster
echo "Splitting, trimming, realigning..."
rm $todir/src/outputFastqs.split* 2>/dev/null
split -l 5 $todir/src/outputFastqs.sh $todir/src/outputFastqs.split 2>/dev/null
for file in $todir/src/outputFastqs.split*
do chmod u+rwx $file
z=${file##*/} ;
bsub -R "select[mem>=10000] rusage[mem=10000]" -M10000 -q verylong -o $todir/log/$z.o  \
-e $todir/src/$z.e $file
done
echo "Done"

## step4 filtering the remapped bam files
## A. remove bad reads using samtools
## this step removes by addition of bit flags of the following types of reads:
## read unmapped (0X4), mate unmapped (0X8), not primary alignment (0X100), read fails platform/vendor quality checks (0X200), read is PCR or optical duplicate (0X400), supplemental alignment (0X800)
## to read about flags in SAM format, see https://broadinstitute.github.io/picard/explain-flags.html
## this site also helps you calculate what flag value to use for a different filtering approach
## B. remove reads with more than 5 bases mismatch (to remove contaminations) using bamutils (from ngsutils: https://github.com/ngsutils/ngsutils)
## this step removes reads with more than 5 mismatches and those reads with mapping quality smaller than Phred score of 50
## odir="/icgc/dkfzlsdf/analysis/B260/projects/przybilm/hipo_K08K"
rm $odir/filtered/src/Filter.* 2>/dev/null
mkdir $odir/filtered 2>/dev/null
mkdir $odir/filtered/src 2>/dev/null
mkdir $odir/filtered/log 2>/dev/null
for file in $odir/MTremapped/*.mt.sorted.bam
do
y=${file%.mt.sorted.bam} ; ## this gets the whole path before .bam
z=${y##*/} ;
echo -e "samtools view -b -F 3852 $file > $odir/filtered/$z.samfiltered.bam; \
samtools index $odir/filtered/$z.samfiltered.bam; \
bamutils filter $odir/filtered/$z.samfiltered.bam $odir/filtered/$z.filtered.bam -lte NM:i 5 -gte MAPQ:i 50; \
samtools index $odir/filtered/$z.filtered.bam" >> $odir/filtered/src/Filter.sh
done

# submit the bash script to the cluster
echo "Perform filtering of remapped reads..."
split -l 5 $odir/filtered/src/Filter.sh $odir/filtered/src/Filter.split 2>/dev/null
for file in $odir/filtered/src/Filter.split*
do echo $file
z=${file##*/} ;
chmod u+rwx $file
bsub -R "select[mem>=10000] rusage[mem=10000]" -M10000 -q verylong -o $odir/filtered/log/$z.o -e $odir/filtered/log/$z.e $file
done
wait
echo "Done"

## step 5 replace read groups (optional - but necessary for variant calling with GATK, see below)
## this is necessary if not all reads from the unit of analysis (cell or donor) are of a single read group. This happens when sequencing is multiplexed?
## to check this, for every bam file (which hopefully contains reads from a cell or a donor), do samtools view -H $bam | grep "@RG" and see whether there are files with more than 1 line
## samtools view -H allows you to see the header of the bams without seeing any of the alignments, and it contains all the information including the read groups (indexed by @RG)
## I strongly recommend (though it's very boring) to read through specification of the SAM format https://samtools.github.io/hts-specs/SAMv1.pdf
## use picardtools to replace readgroups, usually I just replace all the read groups with the cell/donor name
## some options: RGPL=ILLUMINA means this is from illumina sequencing, SORT_ORDER=coordinate gives a sorted bam, CREATE_INDEX=true creates the index too so there's no need to index using samtools
mkdir $odir/filtered_unique_rg 2>/dev/null
mkdir $odir/filtered_unique_rg/log 2>/dev/null
mkdir $odir/filtered_unique_rg/src 2>/dev/null
rm $odir/filtered_unique_rg/src/*replaceRG* 2>/dev/null
for file in $odir/filtered/*.filtered.bam
do echo $file
y=${file%.filtered.bam} ;
z=${y##*/} ;
echo -e "picard AddOrReplaceReadGroups I=$file O=$odir/filtered_unique_rg/$z.bam SORT_ORDER=coordinate RGPU=1 RGID=$z RGLB=$z RGPL=ILLUMINA RGSM=$z CREATE_INDEX=True">>$odir/filtered_unique_rg/src/replaceRG.sh
done

# submit the bash script to the cluster
echo "Replace read groups..."
split -l 5 $odir/filtered_unique_rg/src/replaceRG.sh $odir/filtered_unique_rg/src/replaceRG.split 2>/dev/null
for file in $odir/filtered_unique_rg/src/replaceRG.split*
do echo $file
z=${file##*/} ;
chmod u+rwx $file
bsub -R "select[mem>=10000] rusage[mem=10000]" -M10000 -q verylong -o $odir/filtered_unique_rg/log/$z.o -e $odir/filtered_unique_rg/log/$z.e $file
done
wait
echo "Done"


# *********************************************************************
## VARIANT CALLING (HOMOPLASMIC)
# *********************************************************************

## use GATK to call variants
## GATK 4 onwards (I think) uses a g.vcf approach where variants are called per sample, then consolidated.
## This allows each sample to be processed individually and easy addition of new samples to the same call set.
## the programme used is called HaplotypeCaller, see https://software.broadinstitute.org/gatk/documentation/article?id=11068
## for the mtDNA, we specify the region of the genome -L MT, and --sample-ploidy 1.
## because we specified --sample-ploidy 1, we'd only be good for calling homoplasmic (inherited) variations in the mtDNA
## this setting is NOT GOOD for calling heteroplasmic variations (somatic mutations)
## At this point, the header of the hs37d5.mt.fa file has to be changed so that it contain rCRS instead of MT only
mkdir $odir/gatk 2>/dev/null
mkdir $odir/gatk/log 2>/dev/null
mkdir $odir/gatk/src 2>/dev/null
mkdir $odir/gatk/gvcfs 2>/dev/null
rm $odir/gatk/src/*callGVCFs* 2>/dev/null
for bam in $odir/filtered_unique_rg/*.bam
do echo $bam
y=${bam%.bam} ;
z=${y##*/} ;
echo -e "gatk --java-options '-Xmx1g' HaplotypeCaller \
-R $refdir/hs37d5.mt.fa --sample-ploidy 1 --L MT \
-I $bam -O $odir/gatk/gvcfs/$z.g.vcf.gz -ERC GVCF">>$odir/gatk/src/callGVCFs.sh
done

# submit the bash script to the cluster
echo "Run GATK HaplotypeCaller to call homoplasmic variants..."
split -l 5 $odir/gatk/src/callGVCFs.sh $odir/gatk/src/callGVCFs.split 2>/dev/null
for file in $odir/gatk/src/callGVCFs.split*
do echo $file
z=${file##*/} ;
chmod u+rwx $file
bsub -R "select[mem>=10000] rusage[mem=10000]" -M10000 -q verylong -o $odir/gatk/log/$z.o -e $odir/gatk/log/$z.e $file
done
echo "Done"

## generate the genomedb.list file used afterwards - sample id and full path to corresponding .g.vcf
echo "Generate a genomedb.list..."
Rscript /home/przybilm/bsub_command/WGS/mtdna_pipeline/make_genomedbList.R $odir
echo "Done"

## consolidate all variants called in all bams
## $odir/gatk/gvcfs/genomedb.list is a file you need to create, it is a tab delimited file with <sample ID> <full path to corresponding .g.vcf.gz>
## usually this step is quite fast so I don't even run with cluster, just run the command and wait
gatk --java-options '-Xmx8g' GenomicsDBImport --intervals MT \
-R $refdir/hs37d5.mt.fa --sample-name-map $odir/gatk/gvcfs/genomedb.list \
--genomicsdb-workspace-path $odir/gatk/allsamples --batch-size 100 --consolidate

## genotype all bams at the consolidated variant sites
## This allows you to know when one doesn't have the minor/non-ref allele at a variant site, whether it's due to no data or no variant
## same for this step, usually not worth it to put on cluster to queue, can run interactively
gatk --java-options '-Xmx8g' GenotypeGVCFs --intervals MT --sample-ploidy 1 \
-R $refdir/hs37d5.mt.fa --max-alternate-alleles 2 -new-qual \
-V gendb://$odir/gatk/allsamples \
-O $odir/gatk/allsamples.vcf.gz

## filters variants from all the variants called and genotyped in all samples
## this step allows you to select only SNPs, restrict to only Biallelic SNPs etc, for details see GATK website:
## https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php
gatk --java-options '-Xmx8g' SelectVariants \
-R $refdir/hs37d5.mt.fa \
-V $odir/gatk/allsamples.vcf.gz \
-O $odir/gatk/allsamples.biallelicsnps.vcf.gz \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
--exclude-filtered \
--exclude-non-variants
# --exclude-sample-name $odir/gatk/excludesamples.notinp1p2.txt

## make VCF into plink file that is easier to handle/load into LIMIX etc using PLINK (https://www.cog-genomics.org/plink/2.0/)
## I use keep-allele-order so that the first allele is always ref, second allele is always alt.
## This helps combining several different datasets with variant calling done separately.
## if don't use this option, what happens is plink automatically makes the first allele the minor allele, the second allele the major one
## if two datasets have different major/minor alleles, this would make combining them confusing later on
plink --vcf $odir/gatk/allsamples.biallelicsnps.vcf.gz --keep-allele-order --make-bed --out $odir/gatk/allsamples.biallelicsnps

# *********************************************************************
## VARIANT CALLING (HETEROPLASMIC)
# *********************************************************************

## for calling homoplasmic and heteroplasmic variants, the best way to do this at the moment (other than writing something up ourselves) is using mtDNA-server
## see https://mtdna-server.uibk.ac.at/index.html, http://nar.oxfordjournals.org/lookup/doi/10.1093/nar/gkw247 for details.
## this is a web-based system, but there is a local version downloadable from https://github.com/seppinho/mutation-server
## run mtDNA-server, using --level to control the minimum level of heteroplasmy you're willing to call
## for example, --level 0.01 means I only want to count those sites with 1% of reads supporting a mutation or more
## any sites with a mutation supported by less than 1% of reads I'll regard as sequencing error
## this is dependent on coverage. If you have loads of coverage then this number can be lower.
## there are other filters in mtdna-server, but we don't normally need to specify them as we already filtered the bams earlier for mapping quality etc
## mtdna-server doesn't run for sites where coverage < 50 reads on each strand
mkdir $odir/mtdnaserver 2>/dev/null
mkdir $odir/mtdnaserver/log 2>/dev/null
mkdir $odir/mtdnaserver/src 2>/dev/null
wget https://github.com/seppinho/mutserve/releases/download/v1.3.0/mutserve-1.3.0.jar -P $odir/mtdnaserver/src 2>/dev/null
rm $odir/mtdnaserver/src/variantCalling.s* 2>/dev/null
for file in $odir/filtered_unique_rg/*.bam
do echo $file
y=${file%.bam} ;
z=${y##*/}
echo -e "java -Xmx100g -jar $odir/mtdnaserver/src/mutserve-1.3.0.jar analyse-local \
--input $file --output $odir/mtdnaserver/$z.vcf.gz --reference $refdir/rCRS.fasta --level 0.001 ">>$odir/mtdnaserver/src/variantCalling.sh
done

# submit the bash script to the cluster
echo "Run mtdnaserver to call heteroplasmic variants..."
split -l 5 $odir/mtdnaserver/src/variantCalling.sh $odir/mtdnaserver/src/variantCalling.split 2>/dev/null
for file in $odir/mtdnaserver/src/variantCalling.split*
do echo $file
z=${file##*/} ;
chmod u+rwx $file
bsub -R "select[mem>=20000] rusage[mem=20000]" -M20000 -q verylong -o $odir/mtdnaserver/log/$z.o -e $odir/mtdnaserver/log/$z.e $file
done
echo "Done"

## combine and get all variants
ls $odir/mtdnaserver | grep -v "_raw.txt" | grep -v "_all.var.txt" | grep ".txt" > $odir/txtfile.txt
txtfiles=$(cat $odir/txtfile.txt)
for file in $txtfiles
do echo $file
y=${file%.txt} ;
z=${y##*/}
sed 1d $odir/mtdnaserver/$z.txt >> $odir/mtdnaserver/all.var.txt
done
sed -n '1p' $odir/mtdnaserver/$file > $odir/mtdnaserver/varheader.txt
sed 's/\.bam//' $odir/mtdnaserver/all.var.txt > $odir/mtdnaserver/all.var.tmp
cat $odir/mtdnaserver/varheader.txt $odir/mtdnaserver/all.var.tmp > $odir/mtdnaserver/all.var.txt
