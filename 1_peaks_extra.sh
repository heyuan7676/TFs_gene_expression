#!/bin/bash


####### Note that bedtools assume you use "\t"  as the deliminator!!!!!!!

A=$1
B=$2
C=$3
D=$4

fn_name=$5
fromBEGINING=$6
mergeAB=$7


cd GSE86797_RAW

if ${fromBEGINING}
	then
        cp original_data/* .
        ## merge two replicates
        cat GSM2308860_Memory-1_PairedAlignment.bam_MACS_peaks.bed  GSM2308861_Memory-2_PairedAlignment.bam_MACS_peaks.bed > GSM_Memory.bed
        cat GSM2308856_oNaive-1_PairedAlignment.bam_MACS_peaks.bed  GSM2308857_oNaive-2_PairedAlignment.bam_MACS_peaks.bed > GSM_oNaive.bed
        cat GSM2308858_Effector-1_PairedAlignment.bam_MACS_peaks.bed  GSM2308859_Effector-2_PairedAlignment.bam_MACS_peaks.bed > GSM_Effector.bed
        cat GSM2308862_NT-1_PairedAlignment.sam.bam_MACS_peaks.bed GSM2308863_NT-2_PairedAlignment.bam_MACS_peaks.bed > GSM_NT.bed
        cat GSM2308864_aPDL1-1_PairedAlignment.sam.bam_MACS_peaks.bed GSM2308865_aPDL1-2_PairedAlignment.bam_MACS_peaks.bed > GSM_aPDL1.bed
        mkdir original_data
        mv GSM*MACS_peaks.bed original_data/

        #####  sort all files and merge
        for x in GSM*.bed
        do
            sort -k1,1 -k2,2n ${x} > temp
            cat temp  | tr ' ' '\t' | bedtools merge -c 5 -o mean -i - | awk '{print $1,$2,$3,"MACS_PEAK_"NR,$4}' | tr ' ' '\t' > ${x}
        done
        rm temp


        ### genes on the mm10 genome 
        genefile=mm10.bed
        sorted_gene_file=mm10_sorted.bed

        #####  annotate the gene symbols
        cp original_data/mm10.bed .
        awk 'BEGIN {FS = OFS = "\t"}
        NR == FNR {
          # while reading the 1st file
          # store its records in the array f
          f[$4] = $0
          next
          }
        $1 in f {
          # when match is found
          # print all values
          print f[$1], $0
          }' mm10.bed mm10_gene_symbol.txt > temp


        awk '{print $1,$2,$3,$14, $5}' temp | tr ' ' '\t' > mm10.bed
        cat mm10.bed | grep -v "random" | grep -v "chrUn" | sort -k1,1 -k2,2n -k3,3n | sortBed |  tr ' ' '\t' > ${sorted_gene_file}


        #### myb peaks
        cp original_data/ENCFF316XIC.bed .
        cp original_data/ENCFF345FTL.bed .
        cat ENCFF316XIC.bed ENCFF345FTL.bed > myb_ENCSR000ETR.bed
        rm ENCFF316XIC.bed ENCFF345FTL.bed

        ### sort and merge two replicates
        sort -k1,1 -k2,2n myb_ENCSR000ETR.bed > temp
        cat temp  | tr ' ' '\t' | bedtools merge -c 5 -o mean -i - | awk '{print $1,$2,$3,"MYB_PEAK_"NR,$4}' | tr ' ' '\t' > myb_ENCSR000ETR.bed


        #### cebpb peaks
        cp original_data/ENCFF255COK.bed .
        cp original_data/ENCFF980ZHP.bed .
        cat ENCFF255COK.bed ENCFF980ZHP.bed > cebpb_ENCSR000AIB.bed
        rm ENCFF255COK.bed ENCFF980ZHP.bed

        ### sort and merge two replicates
        sort -k1,1 -k2,2n cebpb_ENCSR000AIB.bed > temp
        cat temp  | tr ' ' '\t' | bedtools merge -c 5 -o mean -i - | awk '{print $1,$2,$3,"CEBPB_PEAK_"NR,$4}' | tr ' ' '\t' > cebpb_ENCSR000AIB.bed



fi




########################################################
##    1.1). extra peaks in A but not in B/C/D
########################################################

extrafile=${fn_name}_extra.bed
sorted_peak_file=${fn_name}_extra_sorted.bed

bedtools intersect -v -a ${A} -b ${B} ${C} ${D} > ${extrafile}   ### -v: find the *not* overlapping regions in B
cat ${extrafile} | tr ' ' '\t' | sortBed -i - > ${sorted_peak_file}
rm ${extrafile}

# echo "Number of extra ATAC peaks in memory cells:"
# wc -l ${sorted_peak_file}




########################################################
##    1.2). extra peaks in A or B but not in C/D
########################################################


if ${mergeAB}
    then
        cat ${A} ${B} > ${A%.bed}${B#GSM}
        sort -k1,1 -k2,2n ${A%.bed}${B#GSM} > temp
        cat temp  | tr ' ' '\t' | bedtools merge -c 5 -o mean -i - | awk '{print $1,$2,$3,"MACS_PEAK_"NR,$4}' | tr ' ' '\t' > ${A%.bed}${B#GSM}
        rm temp

        newA=${A%.bed}${B#GSM}
        extrafile=${fn_name}_extra.bed
        sorted_peak_file=${fn_name}_extra_sorted.bed

        bedtools intersect -v -a ${newA} -b ${C} ${D} > ${extrafile}   ### -v: find the *not* overlapping regions in B
        cat ${extrafile} | tr ' ' '\t' | sortBed -i - > ${sorted_peak_file}
        rm ${extrafile}
fi

# echo "Number of extra ATAC peaks in memory cells:"
# wc -l ${sorted_peak_file}


