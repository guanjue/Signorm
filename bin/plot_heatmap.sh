for info in $(cat /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_range_color_hex.txt)
do
	i=$(echo "$info" | awk -F '-' '{print $1}')
	color=$(echo "$info" | awk -F '-' '{print $3}')
	ylim=$(echo "$info" | awk -F '-' '{print $4}')
	echo $i
	echo $color
	echo $ylim
	plotHeatmap --colorList 'white, '$color -m 'CMP.'$i'.matrix.mat.gz' -out 'CMP.'$i'.matrix.mat.png' --sortRegions no --zMax 0.2 --zMin 0 --yMin 0 --yMax $ylim
	plotHeatmap --colorList 'white, '$color -m 'ERY_fl.'$i'.matrix.mat.gz' -out 'ERY_fl.'$i'.matrix.mat.png' --sortRegions no --zMax 0.2 --zMin 0 --yMin 0 --yMax $ylim
done

for info in $(cat /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_range_color_hex.txt)
do
        i=$(echo "$info" | awk -F '-' '{print $1}')
        color=$(echo "$info" | awk -F '-' '{print $3}')
        ylim=$(echo "$info" | awk -F '-' '{print $4}')
        echo $i
        echo $color
        echo $ylim
        echo $i
        ### compute matrix
        #computeMatrix scale-regions -S 'CMP.'$i'.state_num.bw' 'ERY_fl.'$i'.state_num.bw' -R gencode.vM4.annotation.pc.cmp_ery_fc.bed --beforeRegionStartLength 50000 --regionBodyLength 50000 --afterRegionStartLength 50000 -o 'CMP_ERY_fl.'$i'.matrix.mat.gz' --binSize 1000 --numberOfProcessors max --sortRegions keep --missingDataAsZero --averageTypeBins mean
        ### sort matrix based on state21
        #time computeMatrixOperations sort -m $i'.matrix.mat.gz' -o $i'.matrix.mat.sort.gz' -R state21_sort.bed
        ### plot heatmap
        plotHeatmap --colorList 'white, '$color -m 'CMP_ERY_fl.'$i'.matrix.mat.gz' -out 'CMP_ERY_fl.'$i'.matrix.mat.png' --sortRegions no --zMin 0 --zMax 0.2 --yMin 0
done

for info in $(cat /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_range_color_hex_black.txt)
do
        i=$(echo "$info" | awk -F '-' '{print $1}')
        color=$(echo "$info" | awk -F '-' '{print $3}')
        ylim=$(echo "$info" | awk -F '-' '{print $4}')
        echo $i
        echo $color
        echo $ylim
        echo $i
        ### compute matrix
        #computeMatrix scale-regions -S 'CMP.'$i'.state_num.bw' 'ERY_fl.'$i'.state_num.bw' -R gencode.vM4.annotation.pc.cmp_ery_fc.bed --beforeRegionStartLength 50000 --regionBodyLength 50000 --afterRegionStartLength 50000 -o 'CMP_ERY_fl.'$i'.matrix.mat.gz' --binSize 1000 --numberOfProcessors max --sortRegions keep --missingDataAsZero --averageTypeBins mean
        ### sort matrix based on state21
        #time computeMatrixOperations sort -m $i'.matrix.mat.gz' -o $i'.matrix.mat.sort.gz' -R state21_sort.bed
        ### plot heatmap
        plotHeatmap --colorList 'black, '$color -m 'CMP_ERY_fl.'$i'.matrix.mat.gz' -out 'CMP_ERY_fl.'$i'.matrix.mat.png' --sortRegions no --zMin 0 --zMax 0.2 --yMin 0
done


for info in $(cat /storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_range_color_hex_black.txt)
do
        i=$(echo "$info" | awk -F '-' '{print $1}')
        color=$(echo "$info" | awk -F '-' '{print $3}')
        ylim=$(echo "$info" | awk -F '-' '{print $4}')
        echo $i
        echo $color
        echo $ylim
        echo $i
        ### compute matrix
        #computeMatrix scale-regions -S 'CMP.'$i'.state_num.bw' 'ERY_fl.'$i'.state_num.bw' -R gencode.vM4.annotation.pc.cmp_ery_fc.bed --beforeRegionStartLength 50000 --regionBodyLength 50000 --afterRegionStartLength 50000 -o 'CMP_ERY_fl.'$i'.matrix.mat.gz' --binSize 1000 --numberOfProcessors max --sortRegions keep --missingDataAsZero --averageTypeBins mean
        ### sort matrix based on state21
        #time computeMatrixOperations sort -m $i'.matrix.mat.gz' -o $i'.matrix.mat.sort.gz' -R state21_sort.bed
        ### plot heatmap
        plotHeatmap --colorList 'blue, '$color -m 'ERY_fl.'$i'.matrix.mat.gz' -out 'ERY_fl.'$i'.matrix.mat.png' --sortRegions no --zMin 0 --zMax 0.2 --yMin 0
done



