
for mark in $(cat mark_list.txt)
do
	Rscript /Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/bin/plot_check.R $mark'_check_table.txt' 4 6 $mark'_check_table_ncis_subtract.png'
done

for mark in $(cat mark_list.txt)
do
	Rscript /Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/bin/plot_check.R $mark'_check_table.txt' 2 6 $mark'_check_table_ncis_divide.png'
done

for mark in $(cat mark_list.txt)
do
	Rscript /Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/bin/plot_check.R $mark'_check_table.txt' 3 6 $mark'_check_table_totalmean_subtract.png'
done

for mark in $(cat mark_list.txt)
do
	Rscript /Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/bin/plot_check.R $mark'_check_table.txt' 1 6 $mark'_check_table_totalmean_divide.png'
done

