##################################
script_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/bin/'
analysis_folder='/Volumes/MAC_Data/data/labs/hardison_lab/vision/gene_rnaseq_atac/'

input_folder='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/test_data/input_signal/'
output_folder_t_r_file='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/test_data/signorm_signal_t_r_info/'
output_folder_normed_sig_file='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/test_data/signorm_normed_signal/'
output_folder_r2='/Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/test_data/signorm_normed_signal_r2/'

##################################
	if [ -d "$output_folder_t_r_file" ]; then  
		rm -r $output_folder_t_r_file
	fi
	mkdir $output_folder_t_r_file

	if [ -d "$output_folder_normed_sig_file" ]; then  
		rm -r $output_folder_normed_sig_file
	fi
	mkdir $output_folder_normed_sig_file

	if [ -d "$output_folder_r2" ]; then  
		rm -r $output_folder_r2
	fi
	mkdir $output_folder_r2

##################################

##################################
cd $analysis_folder
### run pipeline
while read LINE
do
	sig1=$(echo "$LINE" | awk '{print $1}')
	sig2=$(echo "$LINE" | awk '{print $2}')
	echo $sig1 $sig2
	### norm sig1 and sig2 based on total reads
	#python $script_folder'total_reads_norm.py' -a $input_folder$sig1 -b $input_folder$sig2 -m mean -o $input_folder$sig1'.totalreads_normed.txt'
	### get ncis t_r matrix
	time python $script_folder'get_ncis_t_a_b.py' -i $input_folder$sig1 -j $input_folder$sig2 -o $output_folder_t_r_file$sig1'_vs_'$sig2'.txt'
	### get scale factor and normalize x-axis signal
	time Rscript $script_folder'signorm.R' $input_folder$sig1 $input_folder$sig2 $output_folder_t_r_file$sig1'_vs_'$sig2'.txt' $output_folder_normed_sig_file$sig1'.norm.txt' var PELT 5 polynorm 50 1000000 0.95 2017 0 1000 3 4 $script_folder
done < $input_folder'info_table.txt'
##################################

### get r2 and r between replicates 
# cp info_table_compare_r2_signorm.txt info_table_compare_r2_totalmean_norm.txt
# :%s/.signorm.txt/.signorm.totalmean.txt/g
# :%s/5end_histone_histone_signorm\///g
# ls 5end_histone_histone_signorm/*.bamtobed5endintersect.signal.histone_ncis_subtract_input.txt.signorm.txt > info_table_compare_r2_signorm.txt
time Rscript $script_folder'get_relictes_r2.R' 'info_table_compare_r2_signorm.txt' $output_folder_normed_sig_file $output_folder_r2 500000
time Rscript $script_folder'get_relictes_r2.R' 'info_table_compare_r2_totalmean_norm.txt' $output_folder_normed_sig_file $output_folder_r2 500000





