#!/bin/bash


####### Note that bedtools assume you use "\t"  as the deliminator!!!!!!!

B=GSM_Memory.bed
A=GSM_NT.bed

C=GSM_oNaive.bed
D=GSM_Effector.bed


fn_name=NT__Memory_oNaive_Effector
fromBEGINING=false
mergeAB=false

promoter_window=1000
enhancer_window=50000


mkdir GSE86797_RAW/cebpb
mkdir GSE86797_RAW/myb

bash 1_peaks_extra.sh ${A} ${B} ${C} ${D} ${fn_name} ${fromBEGINING} ${mergeAB}


newA=${A}
if ${mergeAB}
	then
	newA=${A%.bed}${B#GSM}
fi


echo   > ${fn_name}.output.txt
echo "========================================================================"   >> ${fn_name}.output.txt
echo "## For all genes in mm10"    >> ${fn_name}.output.txt
echo "========================================================================"   >> ${fn_name}.output.txt


bash 2_1_peaks_myb.sh ${newA} ${fn_name} ${promoter_window} ${enhancer_window}  >> ${fn_name}.output.txt



echo   >> ${fn_name}.output.txt
echo "========================================================================"   >> ${fn_name}.output.txt
echo "## Among the genes whose promoters are open and have myb peaks:"   >> ${fn_name}.output.txt
echo "========================================================================"   >> ${fn_name}.output.txt

genebase=${fn_name}_genes_promoter
bash 3_1_peaks_second_filter_cpbeb.sh ${newA} ${fn_name} ${promoter_window} ${enhancer_window} ${genebase}  >> ${fn_name}.output.txt



echo   >> ${fn_name}.output.txt
echo "========================================================================"   >> ${fn_name}.output.txt
echo "## Among the genes whose enhancer are open and have myb peaks:"   >> ${fn_name}.output.txt
echo "========================================================================"   >> ${fn_name}.output.txt

genebase=${fn_name}_genes_enhancer
bash 3_1_peaks_second_filter_cpbeb.sh ${newA} ${fn_name} ${promoter_window} ${enhancer_window} ${genebase}  >> ${fn_name}.output.txt


mkdir GSE86797_RAW/${fn_name}
mv GSE86797_RAW/${fn_name}*bed GSE86797_RAW/${fn_name}/
mv ${fn_name}*txt GSE86797_RAW/${fn_name}/

rm GSE86797_RAW/temp