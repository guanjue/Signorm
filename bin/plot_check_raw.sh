
for mark in $(cat mark_list.txt)
do
	Rscript /Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/bin/plot_check.R $mark'_check_table.txt' 1 1 $mark'_check_table.png'
done

