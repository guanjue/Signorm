############
args = commandArgs(trailingOnly=TRUE)
mark_list_file = args[1]
print(mark_list_file)
mark_list = as.matrix(read.table(mark_list_file, header = F))
print(mark_list)
mean_list = c()
sd_list = c()

for (i in seq(1,dim(mark_list)[1])){
	mark = mark_list[i,]
	print(mark)
	d1=read.table(paste(mark, '.all_sample.txt', sep = ''), header =F)
	### get non-zero mean
	mean_data = mean(d1[,])
	sd_data = sd(d1[,])
	mean_list[i] = mean_data
	sd_list[i] = sd_data
	print(mean_data)
	print(sd_data)
}

#mean_sd_sf = max(mean_list) / mean_list
for (i in seq(1,dim(mark_list)[1])){
	mark = mark_list[i,]
	write.table(c(mean_list[i], sd_list[i]), paste(mark, '.mean_sd_sf.txt', sep = ''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}