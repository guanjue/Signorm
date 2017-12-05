############
args = commandArgs(trailingOnly=TRUE)
mark_list_file = args[1]
print(mark_list_file)
mark_list = as.matrix(read.table(mark_list_file, header = F))
print(mark_list)
median_list = c()

for (i in seq(1,dim(mark_list)[1])){
	mark = mark_list[i,]
	print(mark)
	d1=read.table(paste(mark, '.all_sample.txt', sep = ''), header =F)
	### get non-zero median
	median_data = median(d1[d1!=0,])
	median_list[i] = median_data
	print(median_data)
}

median_sf = max(median_list) / median_list
for (i in seq(1,dim(mark_list)[1])){
	mark = mark_list[i,]
	write.table(median_sf[i], paste(mark, '.median_sf.txt', sep = ''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}