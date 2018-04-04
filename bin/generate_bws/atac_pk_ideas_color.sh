
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod 777 bedToBigBed


time python ~/group/software/signorm/bin/generate_bws/atac_pk_ideas_color.py -a /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/snapshot18/atac_20cell.index.matrix.txt -b 5 -c /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/peak_list.txt -d /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/snapshot18/atac_20cell.function.matrix.txt -e 5 -f /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_list.txt -g /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_atac_sig_color/ideas_state_id_color_name_list.txt -o atac_pk

ls atac_pk.*.bed > atac_pk_ideas_state_color_file_list.txt


for filename in $(cat atac_pk_ideas_state_color_file_list.txt)
do
	ct=$(echo "$filename" | awk -F '.' '{print $2}')
	echo $ct
	sort -k1,1 -k2,2n 'atac_pk.'$ct'.bed' > 'atac_pk.'$ct'.sort.bed'
	./bedToBigBed 'atac_pk.'$ct'.sort.bed' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes 'atac_pk.'$ct'.bb'
	rm 'atac_pk.'$ct'.sort.bed'
done

time python ~/group/software/signorm/bin/generate_bws/check_ccRE_types.py

