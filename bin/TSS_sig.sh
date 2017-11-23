
######### get wg segmentation bed
cat ~/scratch/vision/ERY_fl.h3k4me1rep.116.input.bamtobed5endintersect | awk -F '\t' -v OFS='\t' '{print $1,$2,$3}'> 200_noblack.11_22_2017.bed

#########
### get tss input signal bed
paste 200_noblack.11_22_2017.bed /storage/home/gzx103/group/projects/vision/input/merged_normed_input.rounding.txt > 200_noblack.11_22_2017.bed.merged_normed_input.bed

### get tss bed
### sort bed # and X Y remove M
cat gencode_pc_sort.TSSexp10kb.bed | awk -F 'chr' -v OFS='\t' '{print $2}' | sort -k1,1n -k2,2n | awk -F '\t' -v OFS='\t' '{if ($1!="M" && $1!="Y" && $1!="X") print "chr"$1,$2,$3,$4,$5,$6}' > gencode_pc.TSSexp10kb.sort.chr_num.bed
cat gencode_pc_sort.TSSexp10kb.bed | awk -F 'chr' -v OFS='\t' '{print $2}' | sort -k1,1n -k2,2n | awk -F '\t' -v OFS='\t' '{if ($1=="M") print "chr"$1,$2,$3,$4,$5,$6}' > gencode_pc.TSSexp10kb.sort.chr_M.bed
cat gencode_pc_sort.TSSexp10kb.bed | awk -F 'chr' -v OFS='\t' '{print $2}' | sort -k1,1n -k2,2n | awk -F '\t' -v OFS='\t' '{if ($1=="X") print "chr"$1,$2,$3,$4,$5,$6}' > gencode_pc.TSSexp10kb.sort.chr_X.bed
cat gencode_pc_sort.TSSexp10kb.bed | awk -F 'chr' -v OFS='\t' '{print $2}' | sort -k1,1n -k2,2n | awk -F '\t' -v OFS='\t' '{if ($1=="Y") print "chr"$1,$2,$3,$4,$5,$6}' > gencode_pc.TSSexp10kb.sort.chr_Y.bed
### merge bed # and X Y
cat gencode_pc.TSSexp10kb.sort.chr_num.bed gencode_pc.TSSexp10kb.sort.chr_X.bed gencode_pc.TSSexp10kb.sort.chr_Y.bed > gencode_pc.TSSexp10kb.sort.chr_num_X_Y.bed

### get signal matrix
bedtools map -c 4 -o mean -F 0.9 -a gencode_pc.TSSexp10kb.sort.chr_num_X_Y.bed -b 200_noblack.11_22_2017.bed.merged_normed_input.bed > TSS.200_noblack.11_22_2017.bed.merged_normed_input.bed

for marker in $(cat marker_list.txt)
do
	#########
	### get atac input signal bed (signorm)
	sample_marker_folder='/storage/home/gzx103/group/projects/vision/signorm/'$marker'/'

	ls $sample_marker_folder > $marker'_signorm_list.txt'

	cp gencode_pc.TSSexp10kb.sort.chr_num_X_Y.bed 'gencode_pc.TSSexp10kb.chr_num_X_Y.signorm.'$marker'.txt'
	### get tss sample signal bed
	for file in $(cat $marker'_signorm_list.txt')
	do
		echo $file
		### get tss sample signal bed
		paste 200_noblack.11_22_2017.bed $sample_marker_folder$file > 200_noblack.11_22_2017.bed.tmp_sig.bed
		### get signal matrix
		bedtools map -c 4 -o mean -F 0.9 -a gencode_pc.TSSexp10kb.sort.chr_num_X_Y.bed -b 200_noblack.11_22_2017.bed.tmp_sig.bed > 'TSS.200_noblack.11_22_2017.'$file'.bed'
		### merge to signal matrix
		cat 'TSS.200_noblack.11_22_2017.'$file'.bed' | cut -f 7 > 'TSS.200_noblack.11_22_2017.'$file'.bed.sig'
		paste 'gencode_pc.TSSexp10kb.chr_num_X_Y.signorm.'$marker'.txt' 'TSS.200_noblack.11_22_2017.'$file'.bed.sig' > 'gencode_pc.TSSexp10kb.chr_num_X_Y.signorm.'$marker'.txt.tmp'
		mv 'gencode_pc.TSSexp10kb.chr_num_X_Y.signorm.'$marker'.txt.tmp' 'gencode_pc.TSSexp10kb.chr_num_X_Y.signorm.'$marker'.txt'
		### rm tmp files
		rm 200_noblack.11_22_2017.bed.tmp_sig.bed
		rm 'TSS.200_noblack.11_22_2017.'$file'.bed'
		rm 'TSS.200_noblack.11_22_2017.'$file'.bed.sig'
	done
	#########
	### get atac input signal bed (signorm)
	sample_marker_folder='/storage/home/gzx103/group/projects/vision/totalmean/'$marker'/'

	ls $sample_marker_folder > $marker'_totalmean_list.txt'

	cp gencode_pc.TSSexp10kb.sort.chr_num_X_Y.bed 'gencode_pc.TSSexp10kb.chr_num_X_Y.totalmean.'$marker'.txt'
	### get tss sample signal bed
	for file in $(cat $marker'_totalmean_list.txt')
	do
		echo $file
		### get tss sample signal bed
		paste 200_noblack.11_22_2017.bed $sample_marker_folder$file > 200_noblack.11_22_2017.bed.tmp_sig.bed
		### get signal matrix
		bedtools map -c 4 -o mean -F 0.9 -a gencode_pc.TSSexp10kb.sort.chr_num_X_Y.bed -b 200_noblack.11_22_2017.bed.tmp_sig.bed > 'TSS.200_noblack.11_22_2017.'$file'.bed'
		### merge to signal matrix
		cat 'TSS.200_noblack.11_22_2017.'$file'.bed' | cut -f 7 > 'TSS.200_noblack.11_22_2017.'$file'.bed.sig'
		paste 'gencode_pc.TSSexp10kb.chr_num_X_Y.totalmean.'$marker'.txt' 'TSS.200_noblack.11_22_2017.'$file'.bed.sig' > 'gencode_pc.TSSexp10kb.chr_num_X_Y.totalmean.'$marker'.txt.tmp'
		mv 'gencode_pc.TSSexp10kb.chr_num_X_Y.totalmean.'$marker'.txt.tmp' 'gencode_pc.TSSexp10kb.chr_num_X_Y.totalmean.'$marker'.txt'
		### rm tmp files
		rm 200_noblack.11_22_2017.bed.tmp_sig.bed
		rm 'TSS.200_noblack.11_22_2017.'$file'.bed'
		rm 'TSS.200_noblack.11_22_2017.'$file'.bed.sig'
	done
	#########
done





