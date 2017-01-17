#!/bin/bash


####### Note that bedtools assume you use "\t"  as the deliminator!!!!!!!

A=$1   # <- only for output purpose
fn_name=$2
promoter_window=$3
enhancer_window=$4
genebase=$5

cd GSE86797_RAW



cat myb/${genebase}.txt | sort -k1,1 -k2,2n | tr ' ' '\t' > temp
mv temp myb/${genebase}.txt


###################################################################
##    2). genes whole promoter region overlap with the extra peaks
###################################################################


## find the cloest gene to the peaks' region
sorted_gene_file=${genebase}.txt
sorted_peak_file=${fn_name}_extra_sorted.bed
promoter_f=${fn_name}_promoter.bed

bedtools closest -D ref -t first -a ${sorted_peak_file} -b myb/${sorted_gene_file} > ${promoter_f}

## restrict to genes whose +- 1kb upstream region has overlap with the peaks
awk -v var="$promoter_window" '{if($10<=var && $10>=-var && $10!=-1) print $0}' ${promoter_f} > temp
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
awk -v var="$enhancer_window" 'function abs(x){return ((x < 0.0) ? -x : x)} {$7-=var;$8+=var;$10=abs($10);print $0}' ${promoter_f} | cut -d" " -f6-10 > ${flanking_f}

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

Rscript ../peaks_in_enhancer2.R ${fn_name}
echo "Number of genes with the extra ATAC peaks in promoter and enhancer regions"
cat ${fn_name}_genes_prom+enhac.bed | sed "1d" | awk '{print $1}' | uniq | wc -l


## flanking regions
## cat ${fn_name}_genes_prom+enhac.bed | sed "1d" | awk -v var="$enhancer_window" '{$2-=var;$3+=var;print $0}'  | sed 's/"//g'  | tr ' ' '\t' > ${fn_name}_genes_flanking.bed




################################################################################
##   5). find overlap of promoter/enhancer of these gene with cebpb's peaks
################################################################################


cat ${fn_name}_genes_promoter.bed | sed "1d" | sed 's/"//g'  | tr ' ' '\t' > temp

bedtools intersect -wo -a temp -b cebpb_ENCSR000AIB.bed > ${fn_name}_gene_promoter_cebpb_overlap.bed
echo "Number of promoters regions (whose genes are open in both promoter and enhancer only in "${A%.bed}")"
awk '{print $4}' temp | uniq | wc -l
echo "overlap with cebpb peaks within "$enhancer_window" window"
awk '{print $4}' ${fn_name}_gene_promoter_cebpb_overlap.bed | uniq | wc -l 

echo "Number of genes whose promoter has unique openness"
cat ${fn_name}_gene_promoter_cebpb_overlap.bed | cut -f6-9 | uniq > cebpb/${genebase}_promoter.txt
awk '{print $4}' cebpb/${genebase}_promoter.txt | uniq | wc -l



cat ${fn_name}_genes_enhancer.bed | sed "1d" | sed 's/"//g'  | tr ' ' '\t' > temp

bedtools intersect -wo -a temp -b cebpb_ENCSR000AIB.bed > ${fn_name}_gene_enhancer_cebpb_overlap.bed
echo "Number of enhancers regions (whose genes are open in both promoter and enhancer only in "${A%.bed}")"
awk '{print $4}' temp | uniq | wc -l
echo "overlap with cebpb peaks within "$enhancer_window" window"
awk '{print $4}' ${fn_name}_gene_enhancer_cebpb_overlap.bed | uniq | wc -l 

echo "Number of genes whose enhancer has unique openness"
cat ${fn_name}_gene_enhancer_cebpb_overlap.bed | cut -f6-9 | uniq > cebpb/${genebase}_enhancer.txt
awk '{print $4}' cebpb/${genebase}_enhancer.txt | uniq | wc -l





