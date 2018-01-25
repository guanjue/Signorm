### download data
# scp -r gzx103@aci-i.aci.ics.psu.edu:/storage/home/gzx103/group/projects/vision/nbp_ideas_state17 ./
mkdir sig_bed
mkdir binary_bed
mkdir binary_allbins_bed
mkdir hist
### get threshold
for cm in $(cat cell_marker_list.txt)
do
	echo $cm
	paste ~/group/projects/vision/merged_input/200_noblack.11_22_2017.bed '/storage/home/gzx103/scratch/vision/5end/fisher_p/marknorm/'$cm'.fisher_p.txt.signorm.txt.marknorm.txt' | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,".", $4}' > $cm'.fsm.signal.bed'
	Rscript /storage/home/gzx103/group/software/CD_viewer/bin/plot_density_fdrthresh.R $cm'.fsm.signal.bed' $cm'.fsm.signal.hist.png' T 0.1 0.05
	cat $cm'.fsm.signal.bed.binary.allbins.bed' | awk -F '\t' -v OFS='\t' '{if ($5==1) print 1}'  > $cm'.fsm.signal.bed.binary.bed'
	cat $cm'.fsm.signal.bed.binary.allbins.bed' | awk -F '\t' -v OFS='\t' '{print $5}'  > $cm'.fsm.signal.bed.binary.allbins.txt'
	wc -l $cm'.fsm.signal.bed.binary.bed'
	### mv to folder
	mv $cm'.fsm.signal.hist.png' hist/
	mv $cm'.fsm.signal.bed' sig_bed/
	mv $cm'.fsm.signal.bed.binary.allbins.bed' binary_allbins_bed/
	mv $cm'.fsm.signal.bed.binary.allbins.txt' binary_allbins_bed/
	mv $cm'.fsm.signal.bed.binary.bed' binary_bed/
done


paste binary_allbins_bed/ERY_fl*.fsm.signal.bed.binary.allbins.txt | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7"_"$8}' > ERY_fl.binary_matrix.txt
#paste sig_bed/*.fsm.signal.bed | awk -F '\t' -v OFS='\t' '{print $1,$2,$3,int($4*100)/100"_"int($8*100)/100"_"int($12*100)/100"_"int($16*100)/100"_"int($20*100)/100"_"int($24*100)/100"_"int($28*100)/100"_"int($32*100)/100}' > signal_matrix.txt

### plot density of index set number
Rscript /Volumes/MAC_Data/data/labs/zhang_lab/01projects/CD_viewer/bin/plot_index_set_num.R ERY_fl.binary_matrix.txt ERY_fl.binary_matrix.png T 1

### run CD-viewer
time python /Volumes/MAC_Data/data/labs/zhang_lab/01projects/CD_viewer/bin/get_index_vision_ideas.py




for cm in $(cat cell_marker_list.txt)
do
	cat $cm'.fsm.signal.bed.binary.allbins.bed' | awk -F '\t' -v OFS='\t' '{print $5}'  > $cm'.fsm.signal.bed.binary.allbins.txt'
done

