#!/bin/bash


####### Note that bedtools assume you use "\t"  as the deliminator!!!!!!!

A=$1
fn_name=$2
promoter_window=$3
enhancer_window=$4

sorted_peak_file=${fn_name}_extra_sorted.bed



cd GSE86797_RAW


###################################################################
##    2). genes whole promoter region overlap with the extra peaks
###################################################################


## bedtools closest
## Requie pre-sort
## -D: report the distance with refective of the reference genome
##     upstream peaks have distance < 0
## -t: for genes with the same distance closeness, use the first one


## find the cloest gene to the peaks' region
promoter_f=${fn_name}_promoter.bed
sorted_gene_file=mm10_sorted.bed
bedtools closest -D ref -t first -a ${sorted_peak_file} -b ${sorted_gene_file} > ${promoter_f}

## restrict to genes whose +- 1kb region has overlap with the peaks
awk -v var="$promoter_window" '{if($11<=var && $11>=-var && $11!=-1) print $0}' ${promoter_f} > temp
mv temp ${promoter_f}
echo "Number of genes with the extra peaks in promoter regions: "
awk '{print $9}' ${promoter_f} | uniq | wc -l



###################################################################
##  3.1). obtain the flanking region of these genes
##  3.2). get the peaks that overlap with these flanking regions
###################################################################

## flanking region of +- 20kb
## record the 11th column: distance upstream of promoter

flanking_f=${fn_name}_flanking.bed
awk -v var="$enhancer_window" 'function abs(x){return ((x < 0.0) ? -x : x)} {$7-=var;$8+=var;$11=abs($11);print $0}' ${promoter_f} | cut -d" " -f6-9,11 > ${flanking_f}

## intersect with the original peak file
cat ${flanking_f} | tr ' ' '\t' | bedtools intersect -wao -a - -b ${A} > temp  ## 11th column in number of overlap
mv temp ${flanking_f}


# echo "Confirm the number of genes interact with flanking_f >= number of genes interacting with promoters"
# awk '{print $4}' ${flanking_f} | uniq | wc -l
# awk '{print $9}' ${promoter_f} | uniq |  wc -l




################################################################################
##   4). enhance = flanking - promoter
##       get the genes that have extra peaks in both promoter and enhancer
################################################################################

Rscript ../2_2_peaks_in_enhancer.R ${fn_name}
echo "Number of genes with the extra ATAC peaks in promoter and enhancer regions"
cat ${fn_name}_genes_prom+enhac.bed | awk '{print $1}' | uniq | wc -l


## flanking regions
## cat ${fn_name}_genes_prom+enhac.bed | sed "1d" | awk -v var="$enhancer_window" '{$2-=var;$3+=var;print $0}'  | sed 's/"//g'  | tr ' ' '\t' > ${fn_name}_genes_flanking.bed




################################################################################
##   5). find overlap of promoter/enhancer of these gene with myc's peaks
################################################################################


cat ${fn_name}_genes_promoter.bed | sed "1d" | sed 's/"//g'  | tr ' ' '\t' > temp

bedtools intersect -wo -a temp -b myb_ENCSR000ETR.bed > ${fn_name}_gene_promoter_myc_overlap.bed
echo "Number of promoters regions (whose genes are open in both promoter and enhancer only in "${A%.bed}")"
awk '{print $4}' temp | uniq | wc -l
echo "overlap with myb peaks within "$enhancer_window" window"
awk '{print $4}' ${fn_name}_gene_promoter_myc_overlap.bed | uniq | wc -l 

echo "Number of genes whose promoter has unique openness"
cat ${fn_name}_gene_promoter_myc_overlap.bed | cut -f6-9 | uniq > myb/${fn_name}_genes_promoter.txt
awk '{print $4}' myb/${fn_name}_genes_promoter.txt | uniq | wc -l



cat ${fn_name}_genes_enhancer.bed | sed "1d" | sed 's/"//g'  | tr ' ' '\t' > temp

bedtools intersect -wo -a temp -b myb_ENCSR000ETR.bed > ${fn_name}_gene_enhancer_myc_overlap.bed
echo "Number of enhancers regions (whose genes are open in both promoter and enhancer only in "${A%.bed}")"
awk '{print $4}' temp | uniq | wc -l
echo "overlap with myb peaks within "$enhancer_window" window"
awk '{print $4}' ${fn_name}_gene_enhancer_myc_overlap.bed | uniq | wc -l 

echo "Number of genes whose enhancer has unique openness"
cat ${fn_name}_gene_enhancer_myc_overlap.bed | cut -f6-9 | uniq > myb/${fn_name}_genes_enhancer.txt
awk '{print $4}' myb/${fn_name}_genes_enhancer.txt | uniq | wc -l










