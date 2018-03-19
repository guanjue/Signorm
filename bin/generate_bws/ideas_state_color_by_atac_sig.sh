#PBS -l nodes=2:ppn=8
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -A open

module load gcc
module load python/2.7
module load r/3.4

cd /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_atac_sig_color

#ls /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_label/20ct.*.bed.sort.bed > ideas_file_list.txt

for file in $(cat ideas_file_list.txt)
do
	filename=$(echo "$file" | awk -F '/' '{print $10}')
	ct=$(echo "$filename" | awk -F '.' '{print $2}')
	echo $ct
	time python ~/group/software/signorm/bin/generate_bws/ideas_state_color_by_atac_sig.py -a '/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_label/20ct.'$ct'.bed' -b 5 -c '/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/atac_sig/'$ct'.atac.pkn16.txt' -d 1 -e ideas_state_id_color_name_list.txt -u 16.0 -l 0.0 -s /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes -o 'ideas.v15.8.'$ct'.bed'
	sort -k1,1 -k2,2n 'ideas.v15.8.'$ct'.bed' > 'ideas.v15.8.'$ct'.sort.bed'
	./bedToBigBed 'ideas.v15.8.'$ct'.sort.bed' /storage/home/gzx103/group/projects/vision/input_norm/mm10.chrom.sizes 'ideas.v15.8.'$ct'.bb'
done
