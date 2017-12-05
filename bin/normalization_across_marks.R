args = commandArgs(trailingOnly=TRUE)
mark_list_file = args[1]

mark_list = read.table(mark_list_file, header = F)

median_list = c()
i=1
for (mark in mark_list){
	d1=read.table(paste(mark, '.all_sample.txt', sep = ''), header =F)
	### get non-zero median
	median_data = median(d1[d1!=0,])
	median_list[i] = median_data
	i = i+1
	print(median_data)
}

median_sf = max(median_list) / median_list
i=1
for (mark in mark_list){
	write.table(median_sf[i], paste(mark, '.median_sf.txt', sep = ''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
	i=i+1
}