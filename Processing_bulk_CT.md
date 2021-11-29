Processing bulk BG4 CUT&Tag
================

Bulk CUT&Tag data have been processed as illustrated below. The main steps of the processing include:

-   low quality bases trimming
-   alignemnt
-   du-duplication, collection of stats

-   fragment-size election and generation of tracks (bigwig and bedgraph files)
-   peak calling and peak filtering
-   identification of consensus regions
-   higher level analysis

## Basic processing

``` bash
## ## ## ## ## ## ## ## ## ## ## ## ##
## ## low base quality trimming  ## ##
## ## ## ## ## ## ## ## ## ## ## ## ##

mkdir trimmed
# paired-end cut and tag
for file in *.r1*.fastq
do
fq1=$file
fq2=${fq1/r1/r2}
ls $fq1
ls $fq2
echo "----"
sbatch -o %j.out -e %j.err --mem 16000 --wrap "cutadapt -q 20  -o trimmed/${fq1%%.fastq}.trimmed.fq.gz -p trimmed/${fq2%%.fastq}.trimmed.fq.gz $fq1 $fq2 "
done

## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##  alignment to hg38+ecoliK12   ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ##


cd ./trimmed
mkdir aligned
g='~/reference_genomes/hg38_ecoliK12/hg38_selected_ecoli_K12.fa'
w='~/reference_genomes/hg38/hg38.whitelist.sorted.bed'
  
#aln with bwa
for file in *.r1*trimmed.fq.gz 
  do
  f1=$file
  f2=${file/r1/r2}
  sbatch --time 12:00:00 --mem 16G --wrap "bwa mem -M $g $f1 $f2  | samtools view -Sb -F780 -q 10 -L $w - | samtools sort -@ 8 -   > aligned/${f1%%.trimmed.fq.gz}.hg38.sort.bam"
done


## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##  mark and remove duplicates   ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ##

#remove duplicates
for bam in *merged.bam
do
sbatch --time 12:00:00 --mem 8G --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt | samtools index ${bam%%.bam}.markduplicates.bam"
done


## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##         recover stats         ## ##
## ##       stat2: tot mapped       ## ##
## ##    stat5: nodup mapped        ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ##

for file in *merged.bam
do
bam_hg38=$file
bam_hg38_nodup=${file%%.bam}.markduplicates.bam

sbatch --time 00:05:00 --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%.bam}.stat2"
sbatch --time 00:05:00 --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done

# collect stats by chr
for file in *.merged*bam
do
sbatch --mem 4G --wrap "samtools flagstats $file -O tsv > ${file%%.bam}.flagstat"
sbatch --mem 4G --wrap "samtools idxstats $file > ${file%%.bam}.idxstats"
done

# index bam
for bam in *.markduplicates.bam; do sbatch --mem 4G --wrap "samtools index $bam"; done


## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##        generate tracks        ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ##

# generate bw
 for bam in *.markduplicates.bam; 
 do 
 tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`; 
 scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }'); 
 echo $scal_factor_hg38; 
 echo sbatch --time 06:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"; 
 sbatch -e bw_gen%j.${bam%%.markduplicates.bam}.out --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw" 
 done
 
```

## identify enrichments

Find enrichments using

-   seacr on fragments file &lt; 1000.

Before using seacr, prepare files.

``` bash
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##        sort paired-end by name      ## ##
## ##        select 1k fragments          ## ##
## ##        generate bedgraphs           ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# sort by name
# convert to bedpe
for file in *.merged.markduplicates.bam
do
label=${file%%.merged.markduplicates.bam}
cmd1="samtools sort -n $file > $label.sortName.bam &&\
  bedtools bamtobed -bedpe -i $label.sortName.bam > $label.sortName.bed"
echo $cmd1
sbatch --mem 4G --wrap "$cmd1"
done


# fragment size selection
genome=~/reference_genomes/hg38/hg38_selected.sorted.genome
for file in *sortName.bed
do
cmd_1000="awk '{if (\$1==\$4 && \$6-\$2 < 1000) print \$0}' $file > ${file%%.bed}.1000.clean.bed && cut -f 1,2,6 ${file%%.bed}.1000.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${file%%.bed}.1000.clean.fragments.bed && bedtools genomecov -bg -i ${file%%.bed}.1000.clean.fragments.bed -g $genome > ${file%%.bed}.1000.clean.fragments.bedgraph"
echo $cmd_1000
echo " >> =============== ============= << "
sbatch --mem 1G --wrap "$cmd_1000"
done

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##          call peaks with seacr      ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

#call regions
#macs2
#seacr
cd ~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned
mkdir seacr_no_ctrl

for file in *sortName.bed
do
bdg_1000=${file%%.bed}.1000.clean.fragments.bedgraph
cmd_s_1000="~/applications/SEACR_1.3.sh $bdg_1000 0.05 non stringent seacr_no_ctrl/${bdg_1000%%clean}.0.05fdr"
echo $cmd_s_1000
sbatch --mem 1G --wrap "$cmd_s_1000"
done

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##          filter peaks with          ## ##
## ##          min 8 reads support        ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


# bash script:
./bash_scripts/eval_read_support_seacr.sh


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##     generate consensus regions      ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# before that: 
#- technical consensus
#- biological consensus

#simpler explicit way for MCF7 (less conditions)
mkdir consensus_dir_min8

#50k rep1
multiIntersectBed -i *50k_bg4_MCF7_1st_rep*min8.reduced.bed | awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_50k_bg4_rep1_1000.multi2.bed

multiIntersectBed -i *50k_bg4_MCF7_2nd*min8.reduced.bed | awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_50k_bg4_rep2_1000.multi2.bed

multiIntersectBed -i ./consensus_dir_min8/MCF7_50k_bg4_rep*multi2.bed| awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_50k_bg4_1000.bio2.bed
#10k
multiIntersectBed -i  *MCF7_10k_bg4_MCF7*min8.reduced.bed | awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_10k_bg4.bio2_1000.bed
#5k
multiIntersectBed -i  *MCF7_5k_bg4_MCF7*min8.reduced.bed | awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_5k_bg4.bio2_1000.bed
#1k
multiIntersectBed -i  *MCF7_1k_bg4_MCF7*min8.reduced.bed | awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_1k_bg4.bio2_1000.bed

multiIntersectBed -i *high_salt*min8.reduced.bed | awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_100k_bg4_high_salt_1000.bed
multiIntersectBed -i *low_salt*min8.reduced.bed | awk '{if ($4>=2) print $0}' |  sortBed -i | mergeBed -i - > ./consensus_dir_min8/MCF7_100k_bg4_low_salt_1000.bed


# more complex for U2OS and K562 with more conditions explored
######################################## for seacr ###########################################

cd ~/20210519_winnie_bulk_trimmed
path=~/20210519_winnie_bulk_trimmed
#peaks_folders=(seacr_no_ctrl seacr_with_ctrl)
peaks_folders=(seacr_with_ctrl)
fragm_size=(1000 800 400 200)
cell_n=(100k 50k)
cell=(U2OS K562)
#cell=(K562)
exp_sel=(2L bulk cellNO)

for dir in ${peaks_folders[@]} #1
do

cd $path/$peaks_folders
mkdir consensus_dir

for cell_type in ${cell[@]} #2
do

for N in ${cell_n[@]} #3
do

for exp in ${exp_sel[@]} #4
do

for frag in ${fragm_size[@]} #5
do
biorep=`ls *${cell_type}*${N}*bg4*${exp}*${frag}*0.05fdr.stringent.bed | grep -o 'rep.' | sort | uniq`

for rep in ${biorep[@]} #loop over biol #6
do

bed_to_use=`ls *${cell_type}*${N}*bg4*${exp}*${rep}*${frag}*0.05fdr.stringent.bed`
echo "============================"
wc -l $bed_to_use
echo "============================"

multiIntersectBed -i $bed_to_use | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > ./consensus_dir/${cell_type}_${N}_bg4_${exp}_${rep}_${frag}.multi2.bed # consensus by biorep

done #1

cd consensus_dir
echo " *** consensus_dir -- START***"
wc -l ${cell_type}*${N}*bg4*${exp}*${frag}*.multi2.bed
echo " *** consensus_dir -- END ***"
multiIntersectBed -i ${cell_type}*${N}*bg4*${exp}*${frag}*.multi2.bed | awk '{if($4>=2) print $0}' | sortBed -i | mergeBed -i - > ${cell_type}_${N}_bg4_${exp}_${frag}.bio.bed # consensus by condition (cell-Ncells-fragment)
cd $path/$peaks_folders
echo " ================ . ================"

done #2

done #3

done #4

done #5

done #6



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## ##         assess reproducibility      ## ##
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# extract overalp with OQs and BG4 ChIP peaks
./bash_scripts/eval_overl_seacr_min8.sh 
```

### Comparison across different cell numbers and cell types of G4 signal

``` bash
mkdir ~/20210707_final_sets_regions
cd ~/20210707_final_sets_regions

# copy regions to use
cp ~/20210519_winnie_bulk_trimmed/seacr_no_ctrl/consensus_dir_min8/*100k*bio*bed .
cp ~/20210519_winnie_bulk_trimmed/seacr_no_ctrl/consensus_dir_min8/*50k*bio*bed .
cp ~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned/seacr_no_ctrl/consensus_dir_min8/MCF7_50k_bg4_1000.bio2.bed .
cp ~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned/seacr_no_ctrl/consensus_dir_min8/MCF7_100k_bg4_high_salt_1000.bed .

#create a unified set of regions
cat *bed | sortBed -i -| mergeBed -i - |sortBed -i - > multiCell.100k_50k.sort.bed

## check overlap between 3 cell lines 100k cells
cd ~/20210707_final_sets_regions

intervene venn -i MCF7_100k_bg4_high_salt_1000.bed K562_100k_bg4_bulk_1000.bio.bed U2OS_100k_bg4_bulk_1000.bio.bed -o ./overlap3cells_MCF7_U2OS_K562


## compute coverages looping on reference and bams
dir_for_coverage_output=~/20210707_final_sets_regions/coverage_data_out_3cells
mkdir $dir_for_coverage_output
ref_bed=~/20210707_final_sets_regions/multiCell.100k_50k.sort.bed
cd ~/20210519_winnie_bulk
dir_set=(~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned/ ~/20210519_winnie_bulk_trimmed/)
#for file in ${dir_for_coverage}/*bed
for dir in ${dir_set[@]}
do
cd $dir
for bam in *.merged.markduplicates.bam
do
cmd_coverage="bedtools coverage -a $ref_bed -b ${bam} -counts > $dir_for_coverage_output/${bam%%.bam}.3cellmerge.bed"
echo $cmd_coverage
echo " = = = "
sbatch --time 01:00:00 --mem 4G --wrap "$cmd_coverage"
done
done
```

### Structural analysis prepare data (of the consensus only)

``` bash

# use consensus regions only, generate fasta before and after randomising the bed intervals

for f in ${list_bed[@]}
do
output_name=`basename $f`
echo $output_name
  sbatch --mem 6G --wrap="bedtools getfasta -fi ~/reference_genomes/hg38/hg38_selected.fa -bed $f > ${output_name%%..bed}.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g ~/reference_genomes/hg38/hg38_selected.sorted.genome -seed 10 | bedtools getfasta -fi ~/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_10.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g ~/reference_genomes/hg38/hg38_selected.sorted.genome -seed 20 | bedtools getfasta -fi ~/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_20.fa"
  sbatch --mem 6G --wrap="shuffleBed -i $f -g ~/reference_genomes/hg38/hg38_selected.sorted.genome -seed 30 | bedtools getfasta -fi ~/reference_genomes/hg38/hg38_selected.fa -bed - > ${output_name%%.bed}.shuf_30.fa"

done

# for each consensus
Rscript ~/sequence_hits_analysis_ChIP_C.R consensus.fa consensus.shuf_10.fa consensus.shuf_20.fa consensus.shuf_30.fa
```

### Extract FRIP stats (Fraction of reads in peaks) and FROP (Fraction of reads outside peaks).

For FRIP and FROP use the [sh script](./bash_scripts/stat_extrac_winnie_seacr.sh) that call the [python](./bash_scripts/extract_frip.py) that performs the computations.

The output of those is then processed in R to generate a final table of results using an [R script](./R_scripts/compute_save_frip_other_STATS_winnie.R).

``` bash
### seacr on 1kb fragments
Rscript compute_save_frip_other_STATS_winnie.R /Users/simeon01/Documents/Winnie/202105518_winnie_bulk_SLX-20232/FRIP_output/frip_seacr_no_no_ctrl_all 1000.frip_stats
```

### Calculate fold enrichments of consensus at different genomic loci (TSS, exon, introns ...)

``` bash

conda activate py2_interv
# enrichments at genomic location
# perform fold enrichment
cd ~/20210519_winnie_bulk
workspace_gen=~/reference_genomes/hg38/hg38.whitelist.bed
genomic_annotation=~/reference_genomes/hg38/annotations_genomic_regions_gencode_v28.anno2.bed

out_dir=(~/20210519_winnie_bulk_trimmed/seacr_no_ctrl/consensus_dir_min8 )

#/Users/simeon01/applications/gat
for dir in ${out_dir[@]}
do
  cd $dir
  pwd
  mkdir FE_genomic_annotations_consensus


  for file in *bed
    do
    echo $file
  base_annotation=`basename $genomic_annotation`
  
    cmd_gat="gat-run.py --ignore-segment-tracks --workspace=${workspace_gen}  --annotations=$genomic_annotation --segments=${file} -S ./FE_genomic_annotations_consensus/${file%%.bed}_intervals_${base_annotation%%.bed}.dat -t 4"
    echo ${file%%.bed}_intervals_${base_annotation%%.bed}.dat
    echo $cmd_gat
    sbatch --mem 4G --wrap "$cmd_gat"

  done
done



for file in *dat; do cut -f 1,2,3,8 $file > ${file%%.dat}.reduced.dat; done
#tar -czf FoldEnrichment.fdr.tar.gz *.reduced.dat
```

final coverage - for trimmed only

``` bash
cd ~/20210519_winnie_bulk_trimmed/seacr_no_ctrl/consensus_dir_min8
dir_for_coverage=~/20210519_winnie_bulk_trimmed/coverage_data_final
mkdir $dir_for_coverage
cat **bio** | sortBed -i - | mergeBed -i - >$dir_for_coverage/merged_bio.seacr.min8.consensus.bed

cd $dir_for_coverage

## compute coverages looping on reference and bams
dir_for_coverage=~/20210519_winnie_bulk_trimmed/coverage_data_final
dir_for_coverage_output=~/20210519_winnie_bulk_trimmed/coverage_data_final_out
mkdir $dir_for_coverage_output

cd ~/20210519_winnie_bulk_trimmed

#for file in ${dir_for_coverage}/*bed
for file in ${dir_for_coverage}/*consensus.bed
do
ref_bed=$file
label_to_use=`basename $ref_bed`
for bam in *.merged.markduplicates.bam
do
cmd_coverage="bedtools coverage -a $ref_bed -b ${bam} -counts > $dir_for_coverage_output/${bam%%.bam}.$label_to_use"

echo $cmd_coverage
echo " = = = "

sbatch --time 01:00:00 --mem 4G --wrap "$cmd_coverage"

done
done
```

#################################### extract FRIP

``` bash
cd ~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned/
conda activate py2
bam_file_dir="~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned"
genome_file="~/reference_genomes/hg38/hg38_selected.sorted.genome"
bed_file_dir="~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned/macs2_no_ctrl"
for bam_file in `ls ${bam_file_dir}/*bg4*.merged.markduplicates.bam`
    do
    
    echo "bam file  ====="
    echo $bam_file
    ls -lh $bam_file

    base_bam_name=`basename ${bam_file}`
    temp_base=${base_bam_name%%.merged.markduplicates.bam}
    temp_base2=${temp_base}
    bed_file=${bed_file_dir}/${temp_base2}.sortName.1000._peaks.min8.reduced.bed
    
    echo "peak file ====="
    echo $bed_file
    ls -lh $bed_file

    output_file=${bam_file_dir}/${base_bam_name%%.merged.markduplicates.bam}.frip_stats
    echo "--"
    echo "outputfile ====="
    echo $output_file
    echo "========== .. ========== "

    # call python function that performs the stats extractions
    sbatch --mem 4G --wrap "python2.7 ~/applications/extract_frip.py -i $bam_file -b $bed_file -g $genome_file -o $output_file"
done

conda activate py2
bam_file_dir="~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned"
genome_file="~/reference_genomes/hg38/hg38_selected.sorted.genome"
bed_file_dir="~/20210628_winnie_bulk_trimmed_SLX-20522/trimmed/aligned/seacr_no_ctrl"
for bam_file in `ls ${bam_file_dir}/*bg4*.merged.markduplicates.bam`
    do
    
    echo "bam file  ====="
    echo $bam_file
    ls -lh $bam_file

    base_bam_name=`basename ${bam_file}`
    temp_base=${base_bam_name%%.merged.markduplicates.bam}
    temp_base2=${temp_base}
    bed_file=${bed_file_dir}/${temp_base2}.sortName.1000.clean.fragments.bedgraph.0.05fdr.stringent.min8.reduced.bed
    
    echo "peak file ====="
    echo $bed_file
    ls -lh $bed_file

    output_file=${bam_file_dir}/${base_bam_name%%.merged.markduplicates.bam}.frip_stats_seacr
    echo "--"
    echo "outputfile ====="
    echo $output_file
    echo "========== .. ========== "

    # call python function that performs the stats extractions
    sbatch --mem 4G --wrap "python2.7 ~/applications/extract_frip.py -i $bam_file -b $bed_file -g $genome_file -o $output_file"
done


## local machine
cd /Users/simeon01/Documents/Winnie/20210628_winnie_second_bulk_SLX-20522/macs2_no_ctrl

Rscript /Users/simeon01/Documents/Winnie/process_frip_stats.R /Users/simeon01/Documents/Winnie/20210628_winnie_second_bulk_SLX-20522/macs2_no_ctrl frip_stats

cd /Users/simeon01/Documents/Winnie/20210628_winnie_second_bulk_SLX-20522/seacr_no_ctrl
Rscript /Users/simeon01/Documents/Winnie/process_frip_stats.R /Users/simeon01/Documents/Winnie/20210628_winnie_second_bulk_SLX-20522/seacr_no_ctrl frip_stats_seacr
```

### not sure we used any of this

``` r
library(dplyr)
library(ggplot2)
M_stats <- read.delim('/Users/simeon01/Documents/Winnie/20210707_final_sets_regions/summary_mapping_both_bulk.csv',sep = ",",stringsAsFactors = F)
M_stats$name <- gsub('.merged.stat5','',M_stats$X)
M_stats <- M_stats[,-c(5,6)]
M_stats <- data.frame(M_stats)
M_stats <- M_stats %>% dplyr::filter(.,!grepl('ct4',name)) %>% dplyr::filter(.,!grepl('10X',name)) 

OQs <- read.delim('/Users/simeon01/Documents/Winnie/20210707_final_sets_regions/summary_peak_OQs.csv',sep = ",",stringsAsFactors = F)
OQs$name <- gsub('.sortName.1000._peaks.min8.reduced.bed','',OQs$X)
OQs$name <- gsub('.sortName.*.','',OQs$X)
OQs$perc <- OQs$N_peaks_oqs/OQs$N_peaks


frip_bulk1 <- read.delim('/Users/simeon01/Documents/Winnie/20210707_final_sets_regions/bulk1_seacr_no_ctrl_1000.frip_stats_summary_table.csv',sep = ",",stringsAsFactors = F)
frip_bulk2 <- read.delim('/Users/simeon01/Documents/Winnie/20210707_final_sets_regions/bulk2_seact_no_ctrl_1000.frip_stats_seacr_summary_table.csv',sep = ",",stringsAsFactors = F)
frip_bulk1$name <- gsub('.sortName.1000.*','',frip_bulk1$X)
frip_bulk2$name <- gsub('.frip_stats_seacr','',frip_bulk2$X)

frip_bulk_all <- as.data.frame(rbind(frip_bulk1,frip_bulk2))

#merge all together

Merge <- dplyr::left_join(M_stats,frip_bulk_all,by="name")
Merge <- dplyr::left_join(Merge,OQs,by="name")
head(Merge)
Merge$cell_system[grep('K562',Merge$name)] <- 'K562'
Merge$cell_system[grep('MCF7',Merge$name)] <- 'MCF7'
Merge$cell_system[grep('U2OS',Merge$name)] <- 'U2OS'
Merge$cell_number <- Merge$cell_system
Merge$cell_number[grep('_100k',Merge$name)] <- '100'
Merge$cell_number[grep('_50k',Merge$name)] <- '50'
Merge$cell_number[grep('_10k',Merge$name)] <- '10'
Merge$cell_number[grep('_5k',Merge$name)] <- '5'
Merge$cell_number[grep('_1k',Merge$name)] <- '1'
Merge$cell_number 


Merge %>% select(trimmed_mapped_no_dup,name,cell_system) %>% ggplot(aes(x=name,y=trimmed_mapped_no_dup,fill=cell_system)) + geom_bar(stat = "identity") +
  labs(title = "noDup reads",
       x = "libs",
       y = "map.reads") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))

Merge %>% select(trimmed_mapped_no_dup,tot_reads,name,cell_system,cell_number) %>% ggplot(aes(x=name,y=(tot_reads-trimmed_mapped_no_dup)/tot_reads,fill=cell_system)) + geom_bar(stat = "identity") +
  labs(title = "dup_rate",
       x = "libs",
       y = "map.reads") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))

Merge %>% select(trimmed_mapped_no_dup,tot_reads,name,cell_system,cell_number) %>% ggplot(aes(x=name,y=(tot_reads-trimmed_mapped_no_dup)/tot_reads,fill=cell_number)) + geom_bar(stat = "identity") +
  labs(title = "dup_rate",
       x = "libs",
       y = "map.reads") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))

Merge %>% select(trimmed_mapped_no_dup,tot_reads,name,FRIP,cell_system,cell_number) %>% ggplot(aes(x=name,y=FRIP,fill=cell_system)) + geom_bar(stat = "identity") +
  labs(title = "frip",
       x = "libs",
       y = "frip") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))

Merge %>% select(trimmed_mapped_no_dup,tot_reads,name,FRIP,cell_system,perc) %>% ggplot(aes(x=name,y=perc,fill=cell_system)) + geom_bar(stat = "identity") +
  labs(title = "overl_OQs",
       x = "libs",
       y = "%OQs") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))


Merge %>% select(N_peaks,cell_number,cell_system) %>% ggplot(aes(x=cell_number,y=N_peaks)) + geom_boxplot() + 
  geom_point() + 
  labs(title = "N_peaks_vs_N_cells",
       x = "N_cells",
       y = "N_peaks") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))

Merge %>% select(N_peaks,cell_number,cell_system,trimmed_mapped_no_dup) %>% ggplot(aes(x=cell_number,y=trimmed_mapped_no_dup)) + geom_boxplot() + 
  geom_point() + 
  labs(title = "N_nodupReads_vs_N_cells",
       x = "N_cells",
       y = "N_noDupReads") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))

Merge %>% select(N_peaks,trimmed_mapped_no_dup,cell_number) %>% ggplot(aes(y=N_peaks,x=trimmed_mapped_no_dup,col=cell_number)) + geom_point()


Merge %>%  dplyr::filter(.,grepl('_bg4_',name))%>%  dplyr::filter(.,grepl('U2OS_100k',name)) %>%  dplyr::filter(.,!grepl('bg4_low_salt',name)) %>% select(N_peaks,cell_number,cell_system,trimmed_mapped_no_dup,perc) %>% ggplot(aes(x=perc,y=trimmed_mapped_no_dup,col=cell_number)) + geom_point() + 
  labs(title = "N_nodupReads_vs_N_cells",
       x = "N_cells",
       y = "N_noDupReads") +
  theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10), 
        legend.title = element_text(size = 10))
```
