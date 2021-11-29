Processing SC BG4 CUT&Tag
================

Approximatly 10,000 nuclei were loaded into each lane. Two libraries for each cell line (MFC7, U2OS) and one library of mixed MCF and U2OS cells were prepared and sequenced on an Illumina NextSeq 500 instrument. Resulting fastq files were processed using the cellranger-atac pipeline calling the count function.

The main steps of the analysis included

-   quality control and alignment to hg38.

-   identification of valid barcodes (cells) and subsequent clustering specifically of mixed cell populations samples (both mixed and artifically mixed).

-   identification of clusters identity of cells based on unsupervised classification ofcell clusters' aggregated profile toghether with known samples (individual single cell replicates and bulk samples).

-   evaluation of G4 loci homogeneity across cell population by assessing the number of barcodes (cells) supporting them in each individual cell cluster.
-   G4 loci stratification based on homogeneity (all loci, top 25 percent loci) and subsequent structural (G4 topologies analysis) and spatial (fold enrichment at genomic loci) characterisation (TSS, exons, introns, etc.).

``` bash

# rename libraries
content_file=`ls *contents.csv`
echo $content_file
id_temp=${content_file#SLX*.}
id_final=${id_temp%%.s_1*}
echo $id_final
  # create a new content file for renaming
cut -d ',' -f 1,2,4 $content_file | sed 's/"//g' | awk -v var="$id_final" -F',' '{print $1"."$2"."var"\t" $3"_"$1"."$2"."var}'| sed 's/-N/_N/g' >  ${content_file%%.csv}_for_renaming.csv

while read -r old_name new_name; do
  echo "new_old: "$old_name
  #done < ${content_file%%.csv}_for_renaming.csv
  ls -lth ${old_name}.r_1.fq.gz
  #ls -lth ${old_name}.r_2.fq.gz
  echo "new_names: "$new_name
  mv ${old_name}.s_1.i_1.fq.gz ${new_name}_S1_L001_I1_001.fastq.gz
  mv ${old_name}.s_1.r_1.fq.gz ${new_name}_S1_L001_R1_001.fastq.gz
  mv ${old_name}.s_1.r_2.fq.gz ${new_name}_S1_L001_R2_001.fastq.gz
  mv ${old_name}.s_1.r_3.fq.gz ${new_name}_S1_L001_R3_001.fastq.gz

done < ${content_file%%.csv}_for_renaming.csv


#rename files
for file in *.s_1.i_1.fq.gz
do
mv $file ${file/.s_1.i_1.fq.gz/_S1_L001_I1_001.fastq.gz}

done

for file in *.s_1.r_1.fq.gz
do
mv $file ${file/.s_1.r_1.fq.gz/_S1_L001_R1_001.fastq.gz}
done

for file in *.s_1.r_2.fq.gz
do
mv $file ${file/.s_1.r_2.fq.gz/_S1_L001_R2_001.fastq.gz}
done

for file in *.s_1.r_3.fq.gz
do
mv $file ${file/.s_1.r_3.fq.gz/_S1_L001_R3_001.fastq.gz}
done


## rename fastq files so that they end with the following patterns in order to run the sc atac cellranger pipeline (count)
#*_S1_L001_I1_001.fastq.gz
#*_S1_L001_R1_001.fastq.gz"
#*_S1_L001_R2_001.fastq.gz"
#*_S1_L001_R3_001.fastq.gz"



# fastqc - check sequencign quality
cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521
mkdir fastqc
for file in *.fq.gz
  do
  
  # fastqc
  sbatch --time 02:00:00 --mem 4G --wrap "fastqc --noextract --nogroup -q -o fastqc/ $file"
done


# merge artifically fastq from different cell lines we concatenated the relative fastqs
sbatch --mem 12G --wrap "cat *I1_001.fastq.gz > unified_fq/K562_U2OS_S1_L001_I1_001.fastq.gz"
sbatch --mem 12G --wrap "cat *R1_001.fastq.gz > unified_fq/K562_U2O_S1_L001_R1_001.fastq.gz"
sbatch --mem 12G --wrap "cat *R2_001.fastq.gz > unified_fq/K562_U2O_S1_L001_R2_001.fastq.gz"
sbatch --mem 12G --wrap "cat *R3_001.fastq.gz > unified_fq/K562_U2O_S1_L001_R3_001.fastq.gz"


# align each library and treat artifically mixed library as an additional one

for file in *S1_L001_I1_001.fastq.gz
do
sample=${file%%_S1*}
echo $sample

id_to_use=${file%%_S1*}
id_to_use=`echo $id_to_use | tr "." "_"`
echo $id_to_use

cmd="cellranger-atac count --id $id_to_use \
--fastqs /scratchb/sblab/simeon01/20210525_winnie_secondSC_SLX-20233/SLX-20233 \
--sample $sample \
--reference ~/reference_genomes/sc_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
--jobmode=local \
--localcores=32 \
--localmem=115"

echo $cmd

sbatch -o $sample.out -e $sample.err --export=ALL --nodes=1 --ntasks-per-node=16 --signal=2 --no-requeue --mem=128G --wrap "$cmd"

echo '=========><============'
done


#Reanalised the data by requesting a different number of clusters



cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/

ls -./outs/*.gz
path_mcf7=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/
path_u2OS=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/

cd $path_mcf7
sbatch --mem 32G --wrap "cellranger-atac reanalyze --id=MCF7_SINAH10_reanalysis \
                            --peaks=outs/peaks.bed \
                            --params=MCF7_SINAH10_reanalysis.csv \
                            --reference=/scratcha/sblab/simeon01/reference_genomes/sc_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0\
                            --fragments=$path_mcf7/outs/fragments.tsv.gz"
cd $path_u2OS
sbatch --mem 32G --wrap "cellranger-atac reanalyze --id=U2OS_SINAG10_reanalysis \
                            --peaks=outs/peaks.bed \
                            --params=U2OS_SINAG10_reanalysis.csv \
                            --reference=/scratcha/sblab/simeon01/reference_genomes/sc_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0\
                            --fragments=$path_u2OS/outs/fragments.tsv.gz"
                            
path_mixed=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7
cd $path_mixed

sbatch --mem 32G --wrap "cellranger-atac reanalyze --id=U2OS_MCF7reanalysis \
                            --peaks=outs/peaks.bed \
                            --params=U2OS_MCF7_reanalysis.csv \
                            --reference=/scratcha/sblab/simeon01/reference_genomes/sc_ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0\
                            --fragments=$path_mixed/outs/fragments.tsv.gz"

num_comps,5
max_clusters,2
```

processing the output of the reanalysis

``` bash
# use sinto to extract reads corresponding to barcodes in each cluster

cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/outs/analysis/clustering/graphclust

awk  -F',' '{print $1"\t"$2}' clusters.csv > clusters_for_sinto.txt

# to extract the Reads for a subset of cells can be extracted from a BAM file
bam=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/outs/possorted_bam.bam
cells=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/outs/analysis/clustering/graphclust/clusters_for_sinto.txt
output_folder=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out
mkdir $output_folder
output_folder2=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out2
mkdir $output_folder2
cd $output_folder
sbatch --mem 64G --wrap "sinto filterbarcodes -b $bam -c $cells"
cd $output_folder2
sbatch --mem 64G --wrap "sinto filterbarcodes -b $bam -c $cells -p 8"

for file in *bam
do
sbatch --mem 7G --wrap "samtools view -c $file > ${file%%.bam}.stat5"
sbatch --mem 7G --wrap "samtools index $file"
done

# create tracks and call peaks
mkdir macs2_no_ctrl
for bam in *bam
do
tot_r_hg38=`cat ${bam%%.bam}.stat5`; 
scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }'); 
echo $scal_factor_hg38; 
echo sbatch --time 06:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"; 
#sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw" 
cmd_macs2_bam="macs2 callpeak --keep-dup all -t $bam -p 1e-5 -n macs2_no_ctrl/${bam%%.merged.markduplicates.bam}.p1e-5 --nomodel --extsize 147"
#sbatch --mem 12G --wrap "$cmd_macs2_bam"
echo "====="
done

grep "chr" 2.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | wc -l
#11956
grep "chr" 1.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | wc -l
#6958

for file in *narrowPeak 
do 
grep "chr" $file |awk '{if($9>=800) print $0}' | intersectBed -a - -b $oqs -wa | sort | uniq | wc -l
done

 for bed in $reference_cells/*bed; do  echo $bed;grep "chr" 2.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' 1.bam.p1e-5_peaks.narrowPeak | intersectBed -a - -b $bed -wa  -nonamecheck | sortBed -i - |wc -l; done
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/MCF7.multi2_bio2.bed
3011
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS.normoxia.bio2.bed
4355

for bed in $reference_cells/*bed; do  echo $bed;grep "chr" 1.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | intersectBed -a - -b $bed -wa  -nonamecheck | sortBed -i - |wc -l; done
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//MCF7.multi2_bio2.bed
3011
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_2L_1000.bio.bed
4919
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_bulk_1000.bio.bed
5391
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS.normoxia.bio2.bed
4355


for bed in $reference_cells/*bed; do  echo $bed;grep "chr" 2.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | intersectBed -a - -b $bed -wa  -nonamecheck | sortBed -i - |wc -l; done
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//MCF7.multi2_bio2.bed
3647  # **
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_2L_1000.bio.bed
3842
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_bulk_1000.bio.bed
4303
/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS.normoxia.bio2.bed
3286


# overlap with oqs
oqs=/scratcha/sblab/simeon01/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_annotations_hg38.bed
for file in *narrowPeak ; do echo $file;grep "chr" $file |awk '{if($9>=800) print $0}' | wc -l ; grep "chr" $file |awk '{if($9>=800) print $0}' | intersectBed -a - -b $oqs -wa  -nonamecheck | sortBed -i - |wc -l;  done
#1.bam.p1e-5_peaks.narrowPeak
#6958
#5412
#(5412/6958)
#2.bam.p1e-5_peaks.narrowPeak
#11956 2.bam.p1e-5_peaks.narrowPeak
#5906
#(5906/11956)


## overlap with OQs
oqs=/scratcha/sblab/simeon01/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_annotations_hg38.bed

mcf7_only=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/outs/filtered_peak_bc_matrix/peaks.bed
u2os_only=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/outs/filtered_peak_bc_matrix/peaks.bed

reference_cells=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/ # here there are bed files

wc -l $mcf7_only
# 5513
wc -l $u2os_only
# 19269
# overlap with OQs
grep "chr" $mcf7_only | intersectBed -a - -b $oqs -wa  -nonamecheck | sort| uniq |wc -l
grep "chr" $u2os_only | intersectBed -a - -b $oqs -wa  -nonamecheck | sort| uniq |wc -l

## overlap with individual cell types
for bed in $reference_cells*bed; do  echo $bed;grep "chr" $mcf7_only | intersectBed -a - -b $bed -wa  -nonamecheck | sort| uniq |wc -l; done
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/MCF7.multi2_bio2.bed
#3898
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_100k_bg4_2L_1000.bio.bed
#4100
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_100k_bg4_bulk_1000.bio.bed
#4315
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS.normoxia.bio2.bed
#3401

for bed in $reference_cells*bed; do  echo $bed;grep "chr" $u2os_only | intersectBed -a - -b $bed -wa  -nonamecheck | sort| uniq |wc -l; done
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/MCF7.multi2_bio2.bed
#6338
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_100k_bg4_2L_1000.bio.bed
#14985
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_100k_bg4_bulk_1000.bio.bed
#17830
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS.normoxia.bio2.bed
#9571


## for the library containing the mixed cell population extract the bam relative to each cluster using sinto

cd /scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/outs/analysis/clustering/graphclust
awk  -F',' '{print $1"\t"$2}' clusters.csv > clusters_for_sinto.txt

# to extract the Reads for a subset of cells can be extracted from a BAM file
bam=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/outs/possorted_bam.bam
cells=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/outs/analysis/clustering/graphclust/clusters_for_sinto.txt
output_folder=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out
mkdir $output_folder
cd $output_folder
sbatch --mem 64G --wrap "sinto filterbarcodes -b $bam -c $cells -p 8"


for file in *bam
do
sbatch --mem 7G --wrap "samtools view -c $file > ${file%%.bam}.stat5"
sbatch --mem 7G --wrap "samtools index $file"
done

# create tracks and call peaks
mkdir macs2_no_ctrl
for bam in *bam
do
tot_r_hg38=`cat ${bam%%.bam}.stat5`; 
scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }'); 
echo $scal_factor_hg38; 
echo sbatch --time 06:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"; 
sbatch --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw" 
cmd_macs2_bam="macs2 callpeak --keep-dup all -t $bam -p 1e-5 -n macs2_no_ctrl/${bam%%.merged.markduplicates.bam}.p1e-5 --nomodel --extsize 147"
sbatch --mem 12G --wrap "$cmd_macs2_bam"
echo "====="
done

grep "chr" 2.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | wc -l
#2972
grep "chr" 1.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | wc -l
#10724

reference_cells=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/ 
for bed in $reference_cells/*bed; 
do  echo $bed;grep "chr" 2.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | intersectBed -a - -b $bed -wa  -nonamecheck | sortBed -i - |wc -l; done
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//MCF7.multi2_bio2.bed
#989 **
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_2L_1000.bio.bed
#977
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_bulk_1000.bio.bed
#1087
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS.normoxia.bio2.bed
#924


for bed in $reference_cells/*bed; do  echo $bed;grep "chr" 1.bam.p1e-5_peaks.narrowPeak |awk '{if($9>=800) print $0}' | intersectBed -a - -b $bed -wa  -nonamecheck | sortBed -i - |wc -l; done
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//MCF7.multi2_bio2.bed
#2335
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_2L_1000.bio.bed
#4138
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS_100k_bg4_bulk_1000.bio.bed
#4872
#/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7//U2OS.normoxia.bio2.bed
#3312 ***


# overlap with oqs
oqs=/scratcha/sblab/simeon01/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_annotations_hg38.bed
for file in *narrowPeak ; do echo $file;grep "chr" $file |awk '{if($9>=800) print $0}' | wc -l ; grep "chr" $file |awk '{if($9>=800) print $0}' | intersectBed -a - -b $oqs -wa  -nonamecheck | sortBed -i - |wc -l;  done
#1.bam.p1e-5_peaks.narrowPeak
#10724
#5080
#(5080/10724)
#2.bam.p1e-5_peaks.narrowPeak
#2972
#1511
#(1511/2972)

for file in *narrowPeak ; do echo $file;grep "chr" $file |awk '{if($9>=800) print $0}' > selected_${file} ;  
done
cd 
cd /scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out/macs2_no_ctrl/overlap_clusters/sets
intervene upset -i *bed /scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/*bed --filenames -o pairwise_overlaps
intervene pairwise -i *bed /scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/*bed --filenames --compute frac --htype color

# generate bams by chromomsomes (to speed up computations on large bam)

# ==================================
# U2OS 1st and 2nd
# ==================================

bams_U2OS_1st=/scratchb/sblab/simeon01/20210525_winnie_secondSC_SLX-20233/SLX-20233/SINAF8/U2OS_SINAF8/outs/possorted_bam.bam
bams_U2OS_2st=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/outs/possorted_bam.bam
# generate chromosomes bams
chrS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
bam=possorted_bam.bam
cd /scratchb/sblab/simeon01/20210525_winnie_secondSC_SLX-20233/SLX-20233/SINAF8/U2OS_SINAF8/outs
cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/outs
chrS=(1)
for chr in ${chrS[@]}
do
echo $chr
echo $chr
sbatch --mem 12G --wrap "samtools view -b $bam chr${chr} > ${bam%%.bam}.chr${chr}.bam"
done

#bams_U2OS_1st_chr1=/scratchb/sblab/simeon01/20210525_winnie_secondSC_SLX-20233/SLX-20233/SINAF8/U2OS_SINAF8/outs/possorted_bam.chr1.bam
#bams_U2OS_2st_chr1=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/outs/possorted_bam.chr1.bam

# ==================================
# MCF7 1st and 2nd
# ==================================

bams_MCF1_1st=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/outs/possorted_bam.bam
bams_MCF1_2st=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/MCF7_SINAH1/outs/possorted_bam.bam
chrS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
bam=possorted_bam.bam
#cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/outs
cd /scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/MCF7_SINAH1/outs
for chr in ${chrS[@]}
do
echo $chr
echo $chr
sbatch --mem 12G --wrap "samtools view -b $bam chr${chr} > ${bam%%.bam}.chr${chr}.bam"
done

#bams_MCF1_1st_chr1=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/outs/possorted_bam.chr1.bam
#bams_MCF1_2st_chr1=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/MCF7_SINAH1/outs/possorted_bam.chr1.bam

# ==================================
# MCF7 and U2OS artificially mixed #
# ==================================

bams_cluster1_art_mixed=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/1.bam
bams_cluster2_art_mixed=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/2.bam
cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/
chrS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

for bam in *bam
do
for chr in ${chrS[@]}
do
echo $chr
echo $chr
sbatch --mem 12G --wrap "samtools view -b $bam chr${chr} > ${bam%%.bam}.chr${chr}.bam"
done
done

bams_cluster2_art_mixed_chr1=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/2.chr1.bam
# /scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.BG4ref.bed -b /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/

# ==================================
# MCF7 and U2OS mixed              #
# ==================================

bams_cluster1_mixed=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out/1.bam
bams_cluster2_mixed=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out/1.bam
cd /scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out
chrS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

for bam in *bam
do
for chr in ${chrS[@]}
do
echo $chr
echo $chr
sbatch --mem 12G --wrap "samtools view -b $bam chr${chr} > ${bam%%.bam}.chr${chr}.bam"
done
done
# select only chr1 from all of the bam
bams=($bams_U2OS_1st $bams_U2OS_2st $bams_MCF1_1st $bams_MCF1_2st $bams_cluster1_art_mixed $bams_cluster2_art_mixed $bams_cluster1_mixed $bams_cluster2_mixed)
for bam in ${bams[@]}
do
echo $bam
sbatch --mem 18G --wrap "samtools view -b $bam chr1 > ${bam%%.bam}.chr1.bam"
done

## compute total number of reads in bam (.stat5 file)
cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/outs
sbatch --mem 12G --wrap "samtools view -c possorted_bam.bam > SLX-20521_SINAH10_MCF7_SINAH10.stat5"
cd /scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/MCF7_SINAH1/outs
sbatch --mem 12G --wrap "samtools view -c possorted_bam.bam > SLX-20523_MCF7_SINAH1_outs.stat5"
cd /scratchb/sblab/simeon01/20210525_winnie_secondSC_SLX-20233/SLX-20233/SINAF8/U2OS_SINAF8/outs
sbatch --mem 12G --wrap "samtools view -c possorted_bam.bam > SLX-20233_SINAF8_U2OS_SINAF8.stat5"
cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/outs
sbatch --mem 12G --wrap "samtools view -c possorted_bam.bam > SLX-20521_SINAG10_U2OS_SINAG10.stat5"
cd /scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/
sbatch --mem 12G --wrap "samtools view -c 1.bam  > SLX-20521_U2OS_MCF7_U2OS_MCF7reanalysis_sinto_out__1.stat5"
sbatch --mem 12G --wrap "samtools view -c 2.bam > SLX-20521_U2OS_MCF7_U2OS_MCF7reanalysis_sinto_out__2.stat5"
cd /scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out
sbatch --mem 12G --wrap "samtools view -c 1.bam  > SLX-20523_U2OS_MCF7_SINAG1_sinto_out__1.stat5"
sbatch --mem 12G --wrap "samtools view -c 2.bam > SLX-20523_U2OS_MCF7_SINAG1_sinto_out__2.stat5"


# create output folder
output_coverage_SC=/scratchb/sblab/simeon01/20210709_cov_SC
mkdir $output_coverage_SC


#create reference = based on BG4 (ChIP-seq--> external file) 
cd /scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7
cat *bio2.bed | sortBed -i - | mergeBed -i -> U2OS_MCF7.BG4ref.bed


cd $output_coverage_SC
ref=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.BG4ref.bed
ref1=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.BG4ref.chr1.bed
# coverageCommand
# bams_U2OS_1st
# bams_U2OS_2st
# bams_MCF1_1st
# bams_MCF1_2st
# bams_cluster1_art_mixed
# bams_cluster2_art_mixed
# bams_cluster1_mixed
# bams_cluster2_mixed
sbatch --mem 32G -e U2OS_1st.err --wrap "bedtools coverage -a $ref -b $bams_U2OS_1st -counts > $output_coverage_SC/U2OS_1st.U2OS_MCF7.BG4ref"
sbatch --mem 32G -e U2OS_2st.err --wrap "bedtools coverage -a $ref -b $bams_U2OS_2st -counts > $output_coverage_SC/U2OS_2st.U2OS_MCF7.BG4ref"
sbatch --mem 32G --wrap "bedtools coverage -a $ref -b $bams_MCF1_1st -counts > $output_coverage_SC/MCF1_1st.U2OS_MCF7.BG4ref"
sbatch --mem 32G --wrap "bedtools coverage -a $ref -b $bams_MCF1_2st -counts > $output_coverage_SC/MCF1_2st.U2OS_MCF7.BG4ref"
sbatch --mem 16G --wrap "bedtools coverage -a $ref -b $bams_cluster1_art_mixed -counts > $output_coverage_SC/cluster1_art_mixed.U2OS_MCF7.BG4ref"
sbatch --mem 64G --wrap "bedtools coverage -a $ref1 -b $bams_cluster2_art_mixed -counts > $output_coverage_SC/cluster2_art_mixed.U2OS_MCF7_chr1.BG4ref"
sbatch --mem 16G --wrap "bedtools coverage -a $ref -b $bams_cluster1_mixed -counts > $output_coverage_SC/cluster1_mixed.U2OS_MCF7.BG4ref"
sbatch --mem 16G --wrap "bedtools coverage -a $ref -b $bams_cluster2_mixed -counts > $output_coverage_SC/cluster2_mixed.U2OS_MCF7.BG4ref"

sbatch --mem 32G --wrap "bedtools coverage -a $ref -b $bams_U2OS_1st_chr1 -counts > $output_coverage_SC/U2OS_1st.U2OS_MCF7_chr1.BG4ref"
sbatch --mem 32G --wrap "bedtools coverage -a $ref -b $bams_U2OS_2st_chr1 -counts > $output_coverage_SC/U2OS_2st.U2OS_MCF7_chr1.BG4ref"
sbatch --mem 32G --wrap "bedtools coverage -a $ref -b $bams_MCF1_1st_chr1 -counts > $output_coverage_SC/MCF1_1st.U2OS_MCF7_chr1.BG4ref"
sbatch --mem 32G --wrap "bedtools coverage -a $ref -b $bams_MCF1_2st_chr1 -counts > $output_coverage_SC/MCF1_2st.U2OS_MCF7_chr1.BG4ref"

# count tot reads
sbatch --mem 10G --wrap "samtools view -c $bams_U2OS_1st > ${bams_U2OS_1st%%.bam}.stat5"
sbatch --mem 10G --wrap "samtools view -c $bams_U2OS_2st > ${bams_U2OS_2st%%.bam}.stat5"
sbatch --mem 10G --wrap "samtools view -c $bams_MCF1_1st > ${bams_MCF1_1st%%.bam}.stat5"
sbatch --mem 10G --wrap "samtools view -c $bams_MCF1_2st > ${bams_MCF1_2st%%.bam}.stat5"
sbatch --mem 10G --wrap "samtools view -c $bams_cluster1_art_mixed > ${bams_cluster1_art_mixed%%.bam}.stat5"
sbatch --mem 10G --wrap "samtools view -c $bams_cluster2_art_mixed > ${bams_cluster2_art_mixed%%.bam}.stat5"
sbatch --mem 10G --wrap "samtools view -c $bams_cluster1_mixed > ${bams_cluster1_mixed%%.bam}.stat5"
sbatch --mem 10G --wrap "samtools view -c $bams_cluster2_mixed > ${bams_cluster2_mixed%%.bam}.stat5"

#### alternative: compute coverage by chromosomes and then combine data during import
bam_mcf7_1st=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/outs
bam_mcf7_2nd=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/MCF7_SINAH1/outs
bam_u2OS_1st=/scratchb/sblab/simeon01/20210525_winnie_secondSC_SLX-20233/SLX-20233/SINAF8/U2OS_SINAF8/outs
bam_u2OS_2nd=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/outs
bam_mcf7_u2OS_art=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/
bam_mcf7_u2OS_mix=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out
#all_bams=($bam_mcf7_1st $bam_mcf7_2nd $bam_u2OS_1st $bam_u2OS_2nd $bam_mcf7_u2OS_art $bam_mcf7_u2OS_mix)
all_bams=($bam_mcf7_u2OS_art $bam_mcf7_u2OS_mix)

output_coverage_SC=/scratchb/sblab/simeon01/20210719_cov_SC
#mkdir $output_coverage_SC

chrS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
ref=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.BG4ref.bed

for chr in ${chrS[@]} # chromosome to screen
do
for bam_dir in ${all_bams[@]} # path
do
echo $bam_dir
#echo ${bam_dir##*SLX}
label=`echo ${bam_dir##*SLX} | sed 's/\//_/'g`
echo SLX${label}
cd $bam_dir
bam_to_use=`ls *chr${chr}.bam | grep -v ^Clu*`
for bam in ${bam_to_use[@]}
do
echo $bam_to_use
ls -lh $bam_to_use
echo '-----'
sbatch --mem 3G --wrap "bedtools coverage -a $ref -b $bam -counts > $output_coverage_SC/SLX${label}_${bam%%.bam}.BG4ref"
echo "====="
done
done
done

## calculate coverage using the bulk samples (individual cells type only i.e. MCF7, U2OS)
cd /scratchb/sblab/simeon01/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned
bams=(A672_2_SLX-20522_U2OS_10k_bg4_10X_5th_ref_rep1_tec1.merged.markduplicates.bam A672_2_SLX-20522_U2OS_10k_bg4_10X_5th_ref_rep1_tec1.merged.markduplicates.bam A702_A_SLX-20548_U2OS_10k_bg4_10X_6th_ref_rep2_tec1.merged.markduplicates.bam A703_A_SLX-20550_MCF7_10k_bg4_10X_6th_ref_rep2_tec1.merged.markduplicates.bam)
output_coverage_SC=/scratchb/sblab/simeon01/20210719_cov_SC
for bam in ${bams[@]}
do
echo $bam
ls -lh $bam
echo '-----'
sbatch --mem 3G --wrap "bedtools coverage -a $ref -b $bam -counts > $output_coverage_SC/${bam%%.bam}.BG4ref"
echo "====="
done


### do the same as before on the reference obtained from bulk CUT&Tag ( merge MCF7 100k, U2OS 100k)
# create a reference based on the consensus of the C&T
consensus_MCF7_ct=/scratchb/sblab/simeon01/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned/seacr_no_ctrl/consensus_dir_min8/MCF7_100k_bg4_high_salt_1000.bed
consensus_U2OS_ct=/scratchb/sblab/simeon01/20210519_winnie_bulk_trimmed/seacr_no_ctrl/consensus_dir_min8/U2OS_100k_bg4_bulk_1000.bio.bed
cat $consensus_MCF7_ct $consensus_U2OS_ct| sortBed -i - | mergeBed -i - >  /scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.100k.CTref.bed

refCT=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.100k.CTref.bed

output_coverage_SC_CTref=/scratchb/sblab/simeon01/20210719_cov_SC_ctref

mkdir $output_coverage_SC_CTref

#### alternative: compute coverage by chromosomes and then combine data during import
bam_mcf7_1st=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAH10/MCF7_SINAH10/outs
bam_mcf7_2nd=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/MCF7_SINAH1/outs
bam_u2OS_1st=/scratchb/sblab/simeon01/20210525_winnie_secondSC_SLX-20233/SLX-20233/SINAF8/U2OS_SINAF8/outs
bam_u2OS_2nd=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/SINAG10/U2OS_SINAG10/outs
bam_mcf7_u2OS_art=/scratchb/sblab/simeon01/20210525_winnie_SC_SLX-20521/SLX-20521/U2OS_MCF7/U2OS_MCF7reanalysis/sinto_out/
bam_mcf7_u2OS_mix=/scratchb/sblab/simeon01/20210701_winnie_SC_SLX-20523/SLX-20523/U2OS_MCF7_SINAG1/sinto_out
all_bams=($bam_mcf7_1st $bam_mcf7_2nd $bam_u2OS_1st $bam_u2OS_2nd $bam_mcf7_u2OS_art $bam_mcf7_u2OS_mix)

chrS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
for chr in ${chrS[@]} # chromosome to screen
do
for bam_dir in ${all_bams[@]} # path
do
echo $bam_dir
#echo ${bam_dir##*SLX}
label=`echo ${bam_dir##*SLX} | sed 's/\//_/'g`
echo SLX${label}
cd $bam_dir
bam_to_use=`ls *chr${chr}.bam | grep -v ^Clu*`
for bam in ${bam_to_use[@]}
do
echo $bam_to_use
ls -lh $bam_to_use
echo '-----'
sbatch --mem 3G --wrap "bedtools coverage -a $refCT -b $bam -counts > $output_coverage_SC_CTref/SLX${label}_${bam%%.bam}.CTref"
echo "====="
done
done
done

## calculate coverage for the bulk sample (individual cells only)
cd /scratchb/sblab/simeon01/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned
bams=(A672_2_SLX-20522_U2OS_10k_bg4_10X_5th_ref_rep1_tec1.merged.markduplicates.bam A673_2_SLX-20523_MCF7_10k_bg4_10X_5th_ref_rep1_tec1.merged.markduplicates.bam A702_A_SLX-20548_U2OS_10k_bg4_10X_6th_ref_rep2_tec1.merged.markduplicates.bam A703_A_SLX-20550_MCF7_10k_bg4_10X_6th_ref_rep2_tec1.merged.markduplicates.bam)

refCT=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.100k.CTref.bed
output_coverage_SC_CTref=/scratchb/sblab/simeon01/20210719_cov_SC_ctref
for bam in ${bams[@]}
do
echo $bam
ls -lh $bam
echo '-----'
sbatch --mem 3G --wrap "bedtools coverage -a $refCT -b $bam -counts > $output_coverage_SC_CTref/${bam%%.bam}.CTref"
echo "====="
done


refCT=/scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.100k.CTref.bed
oqs=/scratcha/sblab/simeon01/reference_genomes/OQs/OQ_hits.lifted_hg19_to_hg38_no_annotations_hg38.bed

annotateBed -i $refCT -files $oqs > /scratcha/sblab/simeon01/Data/Reference_For_U2OS_MCF7/U2OS_MCF7.100k.CTref.annotatedForOQs.bed
```

## Prepare data to explore homogeneity/heterogeneity of G4s across indivual cells of clusters

``` bash
# select barcode to use (clustered barcodes)
cut -f 1 clusters_for_sinto.txt| sort | uniq > barcode_to_use.txt

awk '{if ($2==1) print $1}' clusters_for_sinto.txt | sort | uniq > barcodes_cluster1.txt
awk '{if ($2==2) print $1}' clusters_for_sinto.txt | sort | uniq > barcodes_cluster2.txt

grep -f barcodes_cluster1.txt fragments.tsv  > selected_fragments_by_barcodes_cluster1.tsv
grep -f barcodes_cluster2.txt fragments.tsv  > selected_fragments_by_barcodes_cluster2.tsv

# filter peaks that overlap with fragments and count distinct barcodes 
# --> generate a summarised list of G4 loci with number of supporting barcodes
intersectBed -a selected_1.bam.p1e-5_peaks.narrowPeak -b selected_fragments_by_barcodes_cluster1.tsv -wa -wb| groupBy -g 1,2,3,4 -c 14 -o count_distinct > summarise_cluster1.bed
intersectBed -a selected_2.bam.p1e-5_peaks.narrowPeak -b selected_fragments_by_barcodes_cluster2.tsv -wa -wb | groupBy -g 1,2,3,4 -c 14 -o count_distinct > summarise_cluster2.bed

# identify common and specific G4 regions
#G4s in common between cluster1 and cluster2
intersectBed -b summarise_cluster2.bed -a summarise_cluster1.bed -wa -wb > common_summarise_cluster1_cluster2.bed
#G4 regions specific to cluster1
intersectBed -b common_summarise_cluster1_cluster2.bed -a summarise_cluster1.bed -wa -v> specific_summarise_cluster1.bed
#G4 regions specific to cluster2
intersectBed -a summarise_cluster2.bed -b common_summarise_cluster1_cluster2.bed -wa -v> specific_summarise_cluster2.bed

# generate simple 3 columns bed (to try also GREAT annotation tool or gprofiler tool
cut -f 1,2,3 common_summarise_cluster1_cluster2.bed > common_summarise_cluster1_cluster2_forGreat.bed
cut -f 1,2,3 specific_summarise_cluster1.bed > specific_summarise_cluster1_forGreat.bed
cut -f 1,2,3 specific_summarise_cluster2.bed > specific_summarise_cluster2_forGreat.bed

# get promoter coordinates and annotate the G4 sites from above according to promoters (+-1kb)
cd /Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs
#prom=/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS_minus1000.bed
prom=/Users/simeon01/Documents/genomes/hg38/gencode.v28.annotation.gene_1kbaroundTSS.bed
for file in `ls *_summarise_* | grep -v forGreat | grep -v prom`
do
intersectBed -a $prom -b $file -wa -wb | cut -f 1,2,3,4,11 > ${file%%.bed}.prom.bed
cut -f 4 ${file%%.bed}.prom.bed | sed 's/\..*//g' |sed 's/"//g' | sort | uniq > ${file%%.bed}.prom.txt
done

#prom=/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS_minus1000.bed
prom=/Users/simeon01/Documents/genomes/hg38/gencode.v28.annotation.gene_1kbaroundTSS.bed
for file in `ls summarise_cluster*bed`
do
intersectBed -a $prom -b $file -wa -wb | cut -f 1,2,3,4,11 > ${file%%.bed}.prom.bed
cut -f 4 ${file%%.bed}.prom.bed | sed 's/\..*//g' |sed 's/"//g' | sort | uniq > ${file%%.bed}.prom.txt
done
```

``` bash

## for all the bed intervals of the annotated peaks prepare fasta for structural analysis
mkdir bed_to_use
cp common_summarise_cluster1_cluster2.bed ./bed_to_use
cp specific_summarise_cluster1.bed ./bed_to_use
cp specific_summarise_cluster2.bed ./bed_to_use
cp summarise_cluster1.bed ./bed_to_use
cp summarise_cluster2.bed ./bed_to_use


cd /scratchb/sblab/simeon01/20210707_final_sets_regions/bed_to_use
 wc -l *bed
#    491 common_summarise_cluster1_cluster2.bed
#   9810 specific_summarise_cluster1.bed
#   2361 specific_summarise_cluster2.bed
#  10301 summarise_cluster1.bed
#   2850 summarise_cluster2.bed

 
 
# select top quartile (this is different for each individual case)
# for common regions O use the most conservative threshold (the max of the distr)

mkdir bed_top_quartile
awk '{if($4>=35) print $0"\t"$1"\t"$2"\t"$3}' common_summarise_cluster1_cluster2.bed > ./bed_top_quartile/common_summarise_cluster1_cluster2_top_quartile_cl1.bed
awk '{if($10>=28) print $0"\t"$1"\t"$2"\t"$3}' common_summarise_cluster1_cluster2.bed > ./bed_top_quartile/common_summarise_cluster1_cluster2_top_quartile_cl2.bed

awk '{if($5>=10) print $0"\t"$1"\t"$2"\t"$3}' specific_summarise_cluster1.bed > ./bed_top_quartile/specific_summarise_cluster1_top_quartile.bed

awk '{if($5>=6) print $0"\t"$1"\t"$2"\t"$3}' specific_summarise_cluster2.bed > ./bed_top_quartile/specific_summarise_cluster2_top_quartile.bed


for f in `ls *summa*bed | grep -v prom | grep -v Great | grep -v 1000`
do
output_name=`basename $f`
echo $output_name
sbatch --mem 6G --wrap="bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed $f > ${output_name%%.bed}.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome -seed 10 | bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_10.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome -seed 20 | bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_20.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome -seed 30 | bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_30.fa"
echo  '+++++ done ---->'
done

cd /scratchb/sblab/simeon01/20210707_final_sets_regions/bed_to_use/bed_top_quartile
for f in `ls *summa*bed | grep -v prom | grep -v Great | grep -v 1000`
do
output_name=`basename $f`
echo $output_name
sbatch --mem 6G --wrap="bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed $f > ${output_name%%.bed}.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome -seed 10 | bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_10.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome -seed 20 | bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_20.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.sorted.genome -seed 30 | bedtools getfasta -fi /scratcha/sblab/simeon01/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_30.fa"
echo  '+++++ done ---->'
done


# on complete cluster set
Rscript ~/sequence_hits_analysis_ChIP_C.R summarise_cluster2.fa summarise_cluster2.shuf_10.fa summarise_cluster2.shuf_20.fa summarise_cluster2.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R summarise_cluster1.fa summarise_cluster1.shuf_10.fa summarise_cluster1.shuf_20.fa summarise_cluster1.shuf_30.fa

# on set specific and common
Rscript ~/sequence_hits_analysis_ChIP_C.R common_summarise_cluster1_cluster2.fa common_summarise_cluster1_cluster2.shuf_10.fa common_summarise_cluster1_cluster2.shuf_20.fa common_summarise_cluster1_cluster2.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R specific_summarise_cluster1.fa specific_summarise_cluster1.shuf_10.fa specific_summarise_cluster1.shuf_20.fa specific_summarise_cluster1.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R specific_summarise_cluster2.fa specific_summarise_cluster2.shuf_10.fa specific_summarise_cluster2.shuf_20.fa specific_summarise_cluster2.shuf_30.fa



## on the top quartile
Rscript ~/sequence_hits_analysis_ChIP_C.R common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl1.fa common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl1.shuf_10.fa common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl1.shuf_20.fa common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl1.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl2.fa common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl2.shuf_10.fa common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl2.shuf_20.fa common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2_top_quartile_cl2.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R common_summarise_cluster1_cluster2_top_quartile_cl1.fa common_summarise_cluster1_cluster2_top_quartile_cl1.shuf_10.fa common_summarise_cluster1_cluster2_top_quartile_cl1.shuf_20.fa common_summarise_cluster1_cluster2_top_quartile_cl1.shuf_30.fa 
Rscript ~/sequence_hits_analysis_ChIP_C.R common_summarise_cluster1_cluster2_top_quartile_cl2.fa common_summarise_cluster1_cluster2_top_quartile_cl2.shuf_10.fa common_summarise_cluster1_cluster2_top_quartile_cl2.shuf_20.fa common_summarise_cluster1_cluster2_top_quartile_cl2.shuf_30.fa

Rscript ~/sequence_hits_analysis_ChIP_C.R specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_top_quartile.fa specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_top_quartile.shuf_10.fa specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_top_quartile.shuf_20.fa specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_top_quartile.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster2_top_quartile.fa specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster2_top_quartile.shuf_10.fa specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster2_top_quartile.shuf_20.fa specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster2_top_quartile.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R specific_summarise_cluster1_top_quartile.fa specific_summarise_cluster1_top_quartile.shuf_10.fa specific_summarise_cluster1_top_quartile.shuf_20.fa specific_summarise_cluster1_top_quartile.shuf_30.fa
Rscript ~/sequence_hits_analysis_ChIP_C.R specific_summarise_cluster2_top_quartile.fa specific_summarise_cluster2_top_quartile.shuf_10.fa specific_summarise_cluster2_top_quartile.shuf_20.fa specific_summarise_cluster2_top_quartile.shuf_30.fa
```

``` r
## in R ################ ======================================= #########################
setwd('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/bed_to_use')

piePlot <- function(count, categories) {
    dat <- data.frame(count = count, category = categories)
    dat$fraction <- dat$count / sum(dat$count)
    dat$ymax <- cumsum(dat$fraction)
    dat$ymin <- c(0, head(dat$ymax, n = -1))
    dat$label <- factor(paste(dat$category, dat$count), levels = paste(dat$category, dat$count))
    plot <-
        ggplot(dat, aes(
            fill = label, # fill by label not category
            ymax = ymax,
            ymin = ymin,
            xmin = 0,
            xmax = 1
        )) +
        geom_rect() +
        coord_polar(theta = "y") +
        theme(legend.position="top") + theme_void() # no need for labels anymore
    plot
}

generate_structure_plots <- function(dataframe_peaks,outuput_file,outuput_file_pie,outuput_file_fe, label_to_print,labels_legend){
  
  #pie chart
  pie_all_peaks <- ggplot(dataframe_peaks, aes(x = "", y = actual_count, fill = motif_name)) +
  geom_bar(width = 1, stat = "identity", 
           color = "white") +
  coord_polar("y", start = 0) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=labels_legend)+
  #geom_text(aes(label = actual_count), position = position_stack(vjust = 0.5)) +
  ggtitle(label_to_print) + 
  theme(legend.position ="rigth") +
  theme_void()
  
  bar_fold_enrichments <- ggplot(dataframe_peaks, 
       aes(x = motif_name, y = fold_change_over_random, fill = motif_name)) +
  geom_bar(stat = "identity", color = "white",
           width = 0.7) + 
  scale_fill_manual(values=c("skyblue4","tomato3","darkseagreen3","mediumpurple1","lightskyblue","coral","azure3"),
                    label=label_to_print)+
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=0.8) +
  ggtitle(label_to_print) + ylab('fold_enrichm') + 
  theme(legend.position = "none",
        text = element_text(size=12),
        #labels=labels_legend,
        axis.text.x = element_text(angle = 90))
  
  Final_plot <- grid.arrange(pie_all_peaks,bar_fold_enrichments)
  ggsave(file=outuput_file,Final_plot)
  ggsave(file=outuput_file_pie,pie_all_peaks)
  ggsave(file=outuput_file_fe,bar_fold_enrichments)
  rm(Final_plot)
}
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

assign('macs2_cluster1_peaks_motifs', get(load('summarise_cluster1.results_motifs.Rdata')))
macs2_cluster1_peaks_motifs <- macs2_cluster1_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_cluster1_peaks_motifs <- paste0(macs2_cluster1_peaks_motifs$actual_count,' (',macs2_cluster1_peaks_motifs$percentage,'%)')

assign('macs2_cluster2_peaks_motifs', get(load('summarise_cluster2.results_motifs.Rdata')))
macs2_cluster2_peaks_motifs <- macs2_cluster1_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_cluster2_peaks_motifs <- paste0(macs2_cluster2_peaks_motifs$actual_count,' (',macs2_cluster2_peaks_motifs$percentage,'%)')

# common and specific
assign('macs2_common_peaks_motifs', get(load('common_summarise_cluster1_cluster2.results_motifs.Rdata')))
macs2_common_peaks_motifs <- macs2_cluster1_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_common_peaks_motifs <- paste0(macs2_common_peaks_motifs$actual_count,' (',macs2_common_peaks_motifs$percentage,'%)')

assign('macs2_cluster1_spec_peaks_motifs', get(load('specific_summarise_cluster1.results_motifs.Rdata')))
macs2_cluster1_spec_peaks_motifs <- macs2_cluster1_spec_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_cluster1_spec_peaks_motifs <- paste0(macs2_cluster1_spec_peaks_motifs$actual_count,' (',macs2_cluster1_spec_peaks_motifs$percentage,'%)')

assign('macs2_cluster2_spec_peaks_motifs', get(load('specific_summarise_cluster2.results_motifs.Rdata')))
macs2_cluster2_spec_peaks_motifs <- macs2_cluster2_spec_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_cluster2_spec_peaks_motifs <- paste0(macs2_cluster2_spec_peaks_motifs$actual_count,' (',macs2_cluster2_spec_peaks_motifs$percentage,'%)')

generate_structure_plots(macs2_cluster1_peaks_motifs,
                         'macs2_cluster1_peaks_motifs_fin.pdf',
                         'macs2_cluster1_peaks_motifs_pie.pdf',
                         'macs2_cluster1_peaks_motifs_fe.pdf',
                         'macs2_cluster1_peaks',labels_legend_macs2_cluster1_peaks_motifs)
                         
generate_structure_plots(macs2_cluster2_peaks_motifs,
                         'macs2_cluster2_peaks_motifs_fin.pdf',
                         'macs2_cluster2_peaks_motifs_pie.pdf',
                         'macs2_cluster2_peaks_motifs_fe.pdf',
                         'macs2_cluster2_peaks',labels_legend_macs2_cluster2_peaks_motifs)

generate_structure_plots(macs2_common_peaks_motifs,
                         'macs2_common_peaks_motifs_fin.pdf',
                         'macs2_common_peaks_motifs_pie.pdf',
                         'macs2_common_peaks_motifs_fe.pdf',
                         'macs2_common_peaks',labels_legend_macs2_common_peaks_motifs)
                         
generate_structure_plots(macs2_cluster1_spec_peaks_motifs,
                         'macs2_cluster1_spec_peaks_motifs_fin.pdf',
                         'macs2_cluster1_spec_peaks_motifs_pie.pdf',
                         'macs2_cluster1_spec_peaks_motifs_fe.pdf',
                         'macs2_cluster1_spec_peaks',labels_legend_macs2_cluster1_spec_peaks_motifs)
generate_structure_plots(macs2_cluster2_spec_peaks_motifs,
                         'macs2_cluster2_spec_peaks_motifs_fin.pdf',
                         'macs2_cluster2_spec_peaks_motifs_pie.pdf',
                         'macs2_cluster2_spec_peaks_motifs_fe.pdf',
                         'macs2_cluster2_spec_peaks',labels_legend_macs2_cluster2_spec_peaks_motifs)

######## top quartile
# common and specific top quartile

setwd('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/bed_to_use/bed_top_quartile')
assign('macs2_common_top_peaks_motifs', get(load('common_summarise_cluster1_cluster2_top_quartile_cl1.results_motifs.Rdata')))
macs2_common_top_peaks_motifs <- macs2_common_top_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_common_top_peaks_motifs <- paste0(macs2_common_top_peaks_motifs$actual_count,' (',macs2_common_top_peaks_motifs$percentage,'%)')

assign('macs2_cluster1_spec_top_peaks_motifs', get(load('specific_summarise_cluster1_top_quartile.results_motifs.Rdata')))
macs2_cluster1_spec_top_peaks_motifs <- macs2_cluster1_spec_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_cluster1_spec_top_peaks_motifs <- paste0(macs2_cluster1_spec_top_peaks_motifs$actual_count,' (',macs2_cluster1_spec_top_peaks_motifs$percentage,'%)')

assign('macs2_cluster2_spec_top_peaks_motifs', get(load('specific_summarise_cluster2_top_quartile.results_motifs.Rdata')))
macs2_cluster2_spec_top_peaks_motifs <- macs2_cluster2_spec_top_peaks_motifs %>% mutate(percentage = round(100*actual_count/sum(actual_count),1))
labels_legend_macs2_cluster2_spec_top_peaks_motifs <- paste0(macs2_cluster2_spec_top_peaks_motifs$actual_count,' (',macs2_cluster2_spec_top_peaks_motifs$percentage,'%)')

generate_structure_plots(macs2_common_top_peaks_motifs,
                         'macs2_common_top_peaks_motifs_fin.pdf',
                         'macs2_common_top_peaks_motifs_pie.pdf',
                         'macs2_common_top_peaks_motifs_fe.pdf',
                         'macs2_common_top_peaks_peaks',labels_legend_macs2_common_top_peaks_motifs)
                         
generate_structure_plots(macs2_cluster1_spec_top_peaks_motifs,
                         'macs2_cluster1_spec_top_peaks_motifs_fin.pdf',
                         'macs2_cluster1_spec_top_peaks_motifs_pie.pdf',
                         'macs2_cluster1_spec_top_peaks_motifs_fe.pdf',
                         'macs2_cluster1_spec_top',labels_legend_macs2_cluster1_spec_top_peaks_motifs)

generate_structure_plots(macs2_cluster2_spec_top_peaks_motifs,
                         'macs2_cluster2_spec_top_peaks_motifs_fin.pdf',
                         'macs2_cluster2_spec_top_peaks_motifs_pie.pdf',
                         'macs2_cluster2_spec_top_peaks_motifs_fe.pdf',
                         'macs2_cluster2_spec_top',labels_legend_macs2_cluster2_spec_top_peaks_motifs)
```

Functional analysis with fGSEA of the promoter with G4 regions (ranked by number of supporting cells)

``` r
######################### GSEA ######################################

## == fuctions  and packages ==
library(dendextend)
library(dplyr)
library(ggdendro)
library(ggplot2)
library(biomaRt)
library(fgsea)
library(cluster)
library(ComplexHeatmap)
library(reshape)

pathways.hallmark <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c2.cp.v7.0.entrez.gmt.txt")
go_anno <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c5.all.v7.0.entrez.gmt.txt")
short_sequence_motifs <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c3.all.v7.0.entrez.gmt.txt")
oncogenic_signatures <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c6.all.v7.0.entrez.gmt.txt")
immunologic_signatures <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/c7.all.v7.0.entrez.gmt.txt")
gsea.hallmark <- gmtPathways("/Users/simeon01/Documents/genomes/hg38/h.all.v7.1.entrez.gmt.txt")

tss_file <- '/Users/simeon01/Documents/genomes/hg38/gencode.v28.TSS.bed'
source('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/perform_annotation_and_enrichment_analysis.R')

imported_common <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/common_summarise_cluster1_cluster2.prom.bed')


imported_clust1_spec <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/specific_summarise_cluster1.prom.bed', stringsAsFactors = F)
imported_clust2_spec <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/specific_summarise_cluster2.prom.bed', stringsAsFactors = F)

imported_clust1_all <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/summarise_cluster1_1000.prom.bed', stringsAsFactors = F)
imported_clust2_all <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/summarise_cluster2_1000.prom.bed', stringsAsFactors = F)

Table_GSEA_common <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_common$V4),
                                      metric=imported_common$V6)
Table_GSEA_cluster1 <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust1_spec$V4),
                                      metric=imported_clust1_spec$V6)
Table_GSEA_cluster2 <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust2_spec$V4),
                                      metric=imported_clust2_spec$V6)
                                      
Table_GSEA_cluster1_all <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust1_all$V4),
                                      metric=imported_clust1_all$V6)
Table_GSEA_cluster2_all <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust2_all$V4),
                                      metric=imported_clust2_all$V6)



Res_fgsea_GSEA_common <- perform_annotation_and_enrichment_analysis(Table_GSEA_common)
Res_fgsea_GSEA_cluster1 <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1)
Res_fgsea_GSEA_cluster2 <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2)

Res_fgsea_GSEA_cluster1_all <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1_all)
Res_fgsea_GSEA_cluster2_all <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2_all)

explort_fsga_tables_to_txt('common_regions',Res_fgsea_GSEA_common)
explort_fsga_tables_to_txt('cluster1_regions',Res_fgsea_GSEA_cluster1)
explort_fsga_tables_to_txt('cluster2_regions',Res_fgsea_GSEA_cluster2)



## import data obtained using peaks (from cell ranger)
imported_common_bulk <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/common_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1_cluster2.prom.bed')
imported_clust1_spec_bulk <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster1.prom.bed', stringsAsFactors = F)
imported_clust2_spec_bulk <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/specific_summarise_SLX-20523_U2OS_MCF7_SINAG_cluster2.prom.bed', stringsAsFactors = F)


Table_GSEA_common_bulk <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_common_bulk$V4),
                                      metric=imported_common_bulk$V6)
Table_GSEA_cluster1_bulk <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust1_spec_bulk$V4),
                                      metric=imported_clust1_spec_bulk$V6)
Table_GSEA_cluster2_bulk <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust2_spec_bulk$V4),
                                      metric=imported_clust2_spec_bulk$V6)
Res_fgsea_GSEA_common_bulk <- perform_annotation_and_enrichment_analysis(Table_GSEA_common_bulk)
Res_fgsea_GSEA_cluster1_bulk <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1_bulk)
Res_fgsea_GSEA_cluster2_bulk <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2_bulk)

head(Res_fgsea_GSEA_common_bulk$fgseaRes_GO[order(Res_fgsea_GSEA_common_bulk$fgseaRes_GO$pval)],20)
 plotEnrichment(gsea.hallmark[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]],Res_fgsea_GSEA_cluster1_bulk$ranks)
 plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]],Res_fgsea_GSEA_common_bulk$ranks)
 
 
topPathwaysUp <- Res_fgsea_GSEA_common_bulk$fgseaRes_hallmark$pathway[which(Res_fgsea_GSEA_common_bulk$fgseaRes_hallmark$ES > 0)[head(order(Res_fgsea_GSEA_common_bulk$fgseaRes_hallmark$pval), n=10)]]
topPathwaysDown <- Res_fgsea_GSEA_common_bulk$fgseaRes_hallmark$pathway[which(Res_fgsea_GSEA_common_bulk$fgseaRes_hallmark$ES < 0)[head(order(Res_fgsea_GSEA_common_bulk$fgseaRes_hallmark$pval), n=10)]]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gpathways.hallmark[topPathways], Res_fgsea_GSEA_common_bulk$ranks, Res_fgsea_GSEA_common_bulk$fgseaRes_hallmark,gseaParam=0.5)


topPathwaysUp <- Res_fgsea_GSEA_cluster1_bulk$fgseaRes_gseahallmark$pathway[which(Res_fgsea_GSEA_cluster1_bulk$fgseaRes_gseahallmark$ES > 0)[head(order(Res_fgsea_GSEA_cluster1_bulk$fgseaRes_gseahallmark$pval), n=5)]]
topPathwaysDown <- Res_fgsea_GSEA_cluster1_bulk$fgseaRes_gseahallmark$pathway[which(Res_fgsea_GSEA_cluster1_bulk$fgseaRes_gseahallmark$ES < 0)[head(order(Res_fgsea_GSEA_cluster1_bulk$fgseaRes_gseahallmark$pval), n=5)]]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gsea.hallmark[topPathways], Res_fgsea_GSEA_cluster1_bulk$ranks, Res_fgsea_GSEA_cluster1_bulk$fgseaRes_gseahallmark,gseaParam=0.5)


topPathwaysUp <- Res_fgsea_GSEA_cluster2_bulk$fgseaRes_gseahallmark$pathway[which(Res_fgsea_GSEA_cluster2_bulk$fgseaRes_gseahallmark$ES > 0)[head(order(Res_fgsea_GSEA_cluster2_bulk$fgseaRes_gseahallmark$pval), n=5)]]
topPathwaysDown <- Res_fgsea_GSEA_cluster2_bulk$fgseaRes_gseahallmark$pathway[which(Res_fgsea_GSEA_cluster2_bulk$fgseaRes_gseahallmark$ES < 0)[head(order(Res_fgsea_GSEA_cluster2_bulk$fgseaRes_gseahallmark$pval), n=5)]]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gsea.hallmark[topPathways], Res_fgsea_GSEA_cluster2_bulk$ranks, Res_fgsea_GSEA_cluster2_bulk$fgseaRes_gseahallmark,gseaParam=0.5)



plotEnrichment(gsea.hallmark[["REACTOME_FGFR1C_LIGAND_BINDING_AND_ACTIVATION"]],Res_fgsea_GSEA_common_bulk$ranks)


## import data obtained using seacr
imported_common_seacr <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/common_summarise_cluster1_cluster2_seacr.prom.bed')
imported_clust1_spec_seacr <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/specific_summarise_cluster1_seacr.prom.bed', stringsAsFactors = F)
imported_clust2_spec_seacr <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/specific_summarise_cluster2_seacr.prom.bed', stringsAsFactors = F)

Table_GSEA_common_seacr <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_common_seacr$V4),
                                      metric=imported_common_seacr$V6)
Table_GSEA_cluster1_seacr <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust1_spec_seacr$V4),
                                      metric=imported_clust1_spec_seacr$V6)
Table_GSEA_cluster2_seacr <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_clust2_spec_seacr$V4),
                                      metric=imported_clust2_spec_seacr$V6)
Res_fgsea_GSEA_common_seacr <- perform_annotation_and_enrichment_analysis(Table_GSEA_common_seacr)
Res_fgsea_GSEA_cluster1_seacr <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1_seacr)
Res_fgsea_GSEA_cluster2_seacr <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2_seacr)

plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]],Res_fgsea_GSEA_common_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]],Res_fgsea_GSEA_cluster2_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_LATE"]],Res_fgsea_GSEA_cluster1_seacr$ranks)

plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_common_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster2_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster1_seacr$ranks)

HALLMARK_PI3K_AKT_MTOR_SIGNALING
plotEnrichment(gsea.hallmark[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]],Res_fgsea_GSEA_common_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]],Res_fgsea_GSEA_cluster2_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]],Res_fgsea_GSEA_cluster1_seacr$ranks)

HALLMARK_MTORC1_SIGNALING

plotEnrichment(gsea.hallmark[["HALLMARK_MTORC1_SIGNALING"]],Res_fgsea_GSEA_common_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_MTORC1_SIGNALING"]],Res_fgsea_GSEA_cluster2_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_MTORC1_SIGNALING"]],Res_fgsea_GSEA_cluster1_seacr$ranks)

HALLMARK_NOTCH_SIGNALING

plotEnrichment(gsea.hallmark[["HALLMARK_NOTCH_SIGNALING"]],Res_fgsea_GSEA_common_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_NOTCH_SIGNALING"]],Res_fgsea_GSEA_cluster2_seacr$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_NOTCH_SIGNALING"]],Res_fgsea_GSEA_cluster1_seacr$ranks)

## load expression levels 
RNA_cell_atlas <- read.table('/Users/simeon01/Documents/cell_atlas_data/rna_celline.tsv', header = T, sep = "\t", stringsAsFactors = F)
library(dplyr)
tpm_mcf7 <- RNA_cell_atlas[grep('MCF7',RNA_cell_atlas$Cell.line),]%>% dplyr::select(Gene,Gene.name,TPM)
tpm_U2OS <- RNA_cell_atlas[grep('U-2 OS',RNA_cell_atlas$Cell.line),] %>% dplyr::select(Gene,Gene.name,TPM)

Merged_mcf7_u2os <- left_join(tpm_mcf7,tpm_U2OS,by='Gene') %>% mutate(log2FC=log2(TPM.x/TPM.y)) %>% dplyr::filter(!is.na(log2FC)) %>% dplyr::filter(!is.infinite(log2FC)) %>% mutate(log2FC_rank =rank(log2(TPM.x/TPM.y)))
colnames(Merged_mcf7_u2os)[1] <- 'ensembl_gene_id'

Table_GSEA_cluster1_exp <-  Table_GSEA_cluster1 %>% filter(metric > 3) %>% left_join(.,Merged_mcf7_u2os,by='ensembl_gene_id') %>% mutate(metric2=metric*TPM.y) %>% dplyr::select(ensembl_gene_id,metric2) %>% dplyr::filter(!is.na(metric2)) %>% dplyr::rename(metric = metric2)
Res_fgsea_GSEA_cluster1_exp <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1_exp)

Table_GSEA_cluster2_exp <- Table_GSEA_cluster2 %>% filter(metric > 3) %>% left_join(. ,Merged_mcf7_u2os,by='ensembl_gene_id') %>% mutate(metric2=metric*TPM.x) %>% dplyr::select(ensembl_gene_id,metric2) %>% dplyr::filter(!is.na(metric2))%>% dplyr::rename(metric = metric2)
Res_fgsea_GSEA_cluster2_exp <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2_exp)

plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster2_exp$ranks) # *****************
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster1_exp$ranks) # *****************
HALLMARK_OXIDATIVE_PHOSPHORYLATION
plotEnrichment(gsea.hallmark[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],Res_fgsea_GSEA_cluster2_exp$ranks) # *****************

temp1<- Res_fgsea_GSEA_cluster1_exp$fgseaRes_gseahallmark
temp2<- Res_fgsea_GSEA_cluster2_exp$fgseaRes_gseahallmark
 temp_all <- left_join(temp1,temp2,by="pathway")
 A <- cbind(temp_all$pathway,temp_all$ES.x,temp_all$ES.y, temp_all$pval.x, temp_all$pval.y,temp_all$size.x,temp_all$size.y)
 write.table(A,file='temp_all_pathway.csv',sep = ",",quote = F, col.names=NA)

HALLMARK_TNFA_SIGNALING_VIA_NFKB
plotEnrichment(gsea.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],Res_fgsea_GSEA_cluster2_exp$ranks)

HALLMARK_G2M_CHECKPOINT
plotEnrichment(gsea.hallmark[["HALLMARK_G2M_CHECKPOINT"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_G2M_CHECKPOINT"]],Res_fgsea_GSEA_cluster2_exp$ranks)

HALLMARK_KRAS_SIGNALING_DN
plotEnrichment(gsea.hallmark[["HALLMARK_KRAS_SIGNALING_DN"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_KRAS_SIGNALING_DN"]],Res_fgsea_GSEA_cluster2_exp$ranks)


temp1 <- Res_fgsea_GSEA_cluster1_exp$fgseaRes_GO
temp2 <- Res_fgsea_GSEA_cluster2_exp$fgseaRes_GO
temp_all <- left_join(temp1,temp2,by="pathway")

temp_all$pathway[which(temp_all$ES.x>0 & temp_all$ES.y<0 & size.x>=10 & size.y>=10)]
GO_TELOMERASE_HOLOENZYME_COMPLEX_ASSEMBLY
plotEnrichment(go_anno[["GO_TELOMERASE_HOLOENZYME_COMPLEX_ASSEMBLY"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_TELOMERASE_HOLOENZYME_COMPLEX_ASSEMBLY"]],Res_fgsea_GSEA_cluster2_exp$ranks)
GO_SECONDARY_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY
plotEnrichment(go_anno[["GO_SECONDARY_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_SECONDARY_ACTIVE_TRANSMEMBRANE_TRANSPORTER_ACTIVITY"]],Res_fgsea_GSEA_cluster2_exp$ranks)
GO_IRON_ION_BINDING
plotEnrichment(go_anno[["GO_ORGANIC_CYCLIC_COMPOUND_CATABOLIC_PROCESS"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_ORGANIC_CYCLIC_COMPOUND_CATABOLIC_PROCESS"]],Res_fgsea_GSEA_cluster2_exp$ranks)

GO_CHROMATIN_BINDING
plotEnrichment(go_anno[["GO_CHROMATIN_BINDING"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_CHROMATIN_BINDING"]],Res_fgsea_GSEA_cluster2_exp$ranks)

GO_CELL_CYCLE_G2_M_PHASE_TRANSITION
plotEnrichment(go_anno[["GO_CELL_CYCLE_G2_M_PHASE_TRANSITION"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_CELL_CYCLE_G2_M_PHASE_TRANSITION"]],Res_fgsea_GSEA_cluster2_exp$ranks)
GO_ORGANIC_CYCLIC_COMPOUND_CATABOLIC_PROCESS
plotEnrichment(go_anno[["GO_ORGANIC_CYCLIC_COMPOUND_CATABOLIC_PROCESS"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_ORGANIC_CYCLIC_COMPOUND_CATABOLIC_PROCESS"]],Res_fgsea_GSEA_cluster2_exp$ranks)
GO_UBIQUITIN_LIGASE_COMPLEX
plotEnrichment(go_anno[["GO_UBIQUITIN_LIGASE_COMPLEX"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_UBIQUITIN_LIGASE_COMPLEX"]],Res_fgsea_GSEA_cluster2_exp$ranks)
GO_CELL_CYCLE_PHASE_TRANSITION
plotEnrichment(go_anno[["GO_CELL_CYCLE_PHASE_TRANSITION"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_CELL_CYCLE_PHASE_TRANSITION"]],Res_fgsea_GSEA_cluster2_exp$ranks)
GO_INTRACELLULAR_PROTEIN_TRANSPORT
plotEnrichment(go_anno[["GO_INTRACELLULAR_PROTEIN_TRANSPORT"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_INTRACELLULAR_PROTEIN_TRANSPORT"]],Res_fgsea_GSEA_cluster2_exp$ranks)
GO_POSITIVE_REGULATION_OF_CATABOLIC_PROCESS
plotEnrichment(go_anno[["GO_POSITIVE_REGULATION_OF_CATABOLIC_PROCESS"]],Res_fgsea_GSEA_cluster1_exp$ranks)
plotEnrichment(go_anno[["GO_POSITIVE_REGULATION_OF_CATABOLIC_PROCESS"]],Res_fgsea_GSEA_cluster2_exp$ranks)

#imported_common_bulk  imported_clust1_spec_bulk imported_clust2_spec_bulk
Table_GSEA_imported_common_bulk_exp <- Table_GSEA_common_bulk %>% dplyr::filter(metric>3) %>% left_join(.,Merged_mcf7_u2os,by='ensembl_gene_id') %>% mutate(metric2=metric*TPM.y) %>% dplyr::select(ensembl_gene_id,metric2) %>% dplyr::filter(!is.na(metric2)) %>% dplyr::rename(metric = metric2)
Table_GSEA_cluster1_bulk_exp <- Table_GSEA_cluster1_bulk %>% dplyr::filter(metric>3) %>% left_join(.,Merged_mcf7_u2os,by='ensembl_gene_id') %>% mutate(metric2=metric*TPM.y) %>% dplyr::select(ensembl_gene_id,metric2) %>% dplyr::filter(!is.na(metric2)) %>% dplyr::rename(metric = metric2)
Table_GSEA_cluster2_bulk_exp <- Table_GSEA_cluster2_bulk %>% dplyr::filter(metric>3) %>% left_join(.,Merged_mcf7_u2os,by='ensembl_gene_id') %>% mutate(metric2=metric*TPM.y) %>% dplyr::select(ensembl_gene_id,metric2) %>% dplyr::filter(!is.na(metric2)) %>% dplyr::rename(metric = metric2)

Res_fgsea_GSEA_common_bulk_exp <- perform_annotation_and_enrichment_analysis(Table_GSEA_imported_common_bulk_exp)
Res_fgsea_GSEA_cluster1_bulk_exp <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1_bulk_exp)
Res_fgsea_GSEA_cluster2_bulk_exp <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2_bulk_exp)

plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster1_bulk_exp$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster2_bulk_exp$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_common_bulk_exp$ranks)

## load all clusters based on cellranger peak calling
imported_cr_bulk_cluster1 <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/SLX-20523_U2OS_MCF7_SINAG1_peaks_summarise_cluster1.prom.bed')
imported_cr_bulk_cluster2 <- read.table('/Users/simeon01/Documents/Winnie/20210701_winnie_SC_SLX-20523/outs/SLX-20523_U2OS_MCF7_SINAG1_peaks_summarise_cluster2.prom.bed')

boxplot(imported_cr_bulk_cluster1$V11,imported_cr_bulk_cluster2$V11,names=c('cluster1','cluster2'))
Table_GSEA_cluster1_all_bulk <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_cr_bulk_cluster1$V4),
                                      metric=imported_cr_bulk_cluster1$V11)
Table_GSEA_cluster2_all_bulk <- data.frame( ensembl_gene_id= gsub('\\..*','',imported_cr_bulk_cluster2$V4),
                                      metric=imported_cr_bulk_cluster2$V11)
       
Table_GSEA_cluster1_all_bulk_exp <- Table_GSEA_cluster1_all_bulk %>% dplyr::filter(metric>5) %>% left_join(.,Merged_mcf7_u2os,by='ensembl_gene_id') %>% mutate(metric2=metric*TPM.y) %>% dplyr::select(ensembl_gene_id,metric2) %>% dplyr::filter(!is.na(metric2)) %>% dplyr::rename(metric = metric2)
Table_GSEA_cluster2_all_bulk_exp <- Table_GSEA_cluster2_all_bulk %>% dplyr::filter(metric>5) %>% left_join(.,Merged_mcf7_u2os,by='ensembl_gene_id') %>% mutate(metric2=metric*TPM.y) %>% dplyr::select(ensembl_gene_id,metric2) %>% dplyr::filter(!is.na(metric2)) %>% dplyr::rename(metric = metric2)
dim(Table_GSEA_cluster1_all_bulk_exp)
dim(Table_GSEA_cluster2_all_bulk_exp)
                                      
Res_fgsea_GSEA_cluster1_all_bulk_exp <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1_all_bulk_exp)
Res_fgsea_GSEA_cluster2_all_bulk_exp <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2_all_bulk_exp)

plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster1_all_bulk_exp$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster2_all_bulk_exp$ranks)

# select only genes specific
Table_GSEA_cluster1_all_bulk_exp_spec <- anti_join(Table_GSEA_cluster1_all_bulk_exp,Table_GSEA_cluster2_all_bulk_exp,by='ensembl_gene_id')
Table_GSEA_cluster2_all_bulk_exp_spec <- anti_join(Table_GSEA_cluster2_all_bulk_exp,Table_GSEA_cluster1_all_bulk_exp,by='ensembl_gene_id')

Res_fgsea_GSEA_cluster1_all_bulk_exp_spec <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster1_all_bulk_exp_spec)
Res_fgsea_GSEA_cluster2_all_bulk_exp_spec <- perform_annotation_and_enrichment_analysis(Table_GSEA_cluster2_all_bulk_exp_spec)

plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster1_all_bulk_exp_spec$ranks)
plotEnrichment(gsea.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]],Res_fgsea_GSEA_cluster2_all_bulk_exp_spec$ranks)



### load the summarise files
macs2_clust1 <- read.table('summarise_cluster1.bed',stringsAsFactors = F)
macs2_clust2 <- read.table('summarise_cluster2.bed',stringsAsFactors = F)


#     671 barcodes_cluster1.txt
#     467 barcodes_cluster2.txt
```

Fold enrichments of peaks at genomic regions

``` bash


conda activate py2_interv
# enrichments at genomic location
# perform fold enrichment

workspace_gen=/scratcha/sblab/simeon01/reference_genomes/hg38/hg38.whitelist.bed
genomic_annotation=/scratcha/sblab/simeon01/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.bed


conda activate py2_interv
out_dir=(/scratchb/sblab/simeon01/20210707_final_sets_regions/bed_to_use/bed_top_quartile /scratchb/sblab/simeon01/20210707_final_sets_regions/bed_to_use)
for dir in ${out_dir[@]}
do
  cd $dir
  pwd
  mkdir FE_genomic_annotations


  for file in *bed
    do
    echo $file
  base_annotation=`basename $genomic_annotation`
  
    cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$genomic_annotation --segments=${file} -S ./FE_genomic_annotations/${file%%.bed}_intervals_${base_annotation%%.bed}.dat -t 4"
    echo ${file%%.bed}_intervals_${base_annotation%%.bed}.dat
    echo $cmd_gat
    sbatch --mem 4G --wrap "$cmd_gat"

  done
done



for file in *dat; do cut -f 1,2,3,8 $file > ${file%%.dat}.reduced.dat; done
#tar -czf FoldEnrichment.fdr.tar.gz *.reduced.dat
```
