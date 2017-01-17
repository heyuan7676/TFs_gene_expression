

args <- commandArgs(TRUE)
fn = args[1]

all_peaks_in_flanking = read.table(paste0(fn,"_flanking.bed"), header=F,stringsAsFactors=F)
peaks_in_promoters = read.table(paste0(fn,"_promoter.bed"), header=F,stringsAsFactors=F)
colnames(all_peaks_in_flanking) = c("chr", "start", "end", "gene", "distance_to_peak", "peak_chr", "peak_start", "peak_end", "peak", "peakReads","overlap_peak_gene")
colnames(peaks_in_promoters) = c("peak_chr", "peak_start", "peak_end", "peak", "peakReads", "chr", "start", "end", "gene", "NA","distance_to_peak")


peaks_for_flanking_regions = by(all_peaks_in_flanking,all_peaks_in_flanking$gene,function(x) x[9])
peaks_for_promoter_regions = by(peaks_in_promoters,peaks_in_promoters$gene, function(x) x[4])


### enhancer = all - promoter

peaks_for_enhancer_regions = list()  

for(n in names(peaks_for_flanking_regions)){
	a = peaks_for_flanking_regions[[n]]
	b = peaks_for_promoter_regions[[n]]
	peaks_in_enhancer = a[! a %in% b]
	if(length(peaks_in_enhancer) > 0){
		peaks_for_enhancer_regions[[n]] = peaks_in_enhancer
	}
}


colnames_save = c("peak_chr", "peak_start", "peak_end", "peak", "peakReads", "chr", "start", "end", "gene")

### promoter for these genes

genes_for_promoter = peaks_in_promoters$gene
idx = which(genes_for_promoter %in% names(peaks_for_enhancer_regions))

peaks_in_promoters = peaks_in_promoters[idx,]
write.table(peaks_in_promoters[,colnames_save], paste0(fn,"_genes_promoter.bed"), row.names=F,sep='\t')



### enhancer for these genes

genes_for_enhance = all_peaks_in_flanking$gene
idx = which(genes_for_enhance %in% names(peaks_for_enhancer_regions))

peaks_in_enhancer = all_peaks_in_flanking[idx,]
write.table(peaks_in_enhancer[,colnames_save], paste0(fn,"_genes_enhancer.bed"), row.names=F,sep='\t')



### genes

genes = unique(all_peaks_in_flanking[,c("gene","start","end")])
idx = which(genes$gene %in% names(peaks_for_enhancer_regions))

peaks_in_enhancer = genes[idx,]
write.table(genes, paste0(fn,"_genes_prom+enhac.bed"), row.names=F,sep='\t')




