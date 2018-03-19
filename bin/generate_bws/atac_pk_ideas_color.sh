

ls /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_label/20ct.*.bed.sort.bed > ideas_file_list.txt

for file in $(cat ideas_file_list.txt)
do
	filename=$(echo "$file" | awk -F '/' '{print $10}')
	ct=$(echo "$filename" | awk -F '.' '{print $2}')
	echo $ct
	time python ~/group/software/signorm/bin/generate_bws/ideas_state_color_by_atac_sig.py -a '/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_label/20ct.'$ct'.bed' -b 5 -c '/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/atac_sig/'$ct'.atac.pkn16.txt' -d 1 -e ideas_state_id_color_name_list.txt -u 16.0 -l 0.0 -s /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes -o 'ideas.v15.8.'$ct'.bed'
	sort -k1,1 -k2,2n 'ideas.v15.8.'$ct'.bed' > 'ideas.v15.8.'$ct'.sort.bed'
	./bedToBigBed 'ideas.v15.8.'$ct'.sort.bed' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes 'ideas.v15.8.'$ct'.bb'
done


time python ~/group/software/signorm/bin/generate_bws/ideas_state_color_by_atac_sig.py -a /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_label/20ct.B_SPL.bed -b 5 -c /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/atac_sig/B_SPL.atac.pkn16.txt -d 1 -e ideas_state_id_color_name_list.txt -u 16.0 -l 0.0 -o ideas.v15.8.B_SPL.bed

done



wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed

chmod 777 bedToBigBed

./bedToBigBed ideas.v15.8.B_SPL.bed chrom.sizes ideas.v15.8.B_SPL.bb



