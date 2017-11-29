library(LSD)

### get parameters
args = commandArgs(trailingOnly=TRUE)
input_table = args[1]
col_num = as.numeric(args[2])
col_num_dist = as.numeric(args[3])
output_filename = args[4]

data = read.table(input_table, header = F)
used_id = seq(col_num, dim(data)[1], col_num_dist)
data=data[used_id,]

header_name = c('s1', 's2', 'p_cor_all', 's_cor_all', 'MSE', 'p_cor_non0', 's_cor_non0', 'R2', 'total_reads_d1', 'total_reads_d2')

colnames(data) = header_name

data = data[(data[,3]!=1),]

total_r = seq(dim(data)[1])#data[, 9]+data[, 10]

O = order(total_r)

data_sort = data[,]

png(output_filename)
par(mfrow=c(2,2))
plot(total_r, data_sort[,3], main=header_name[3], ylim = c(0,1))
abline(h=0.5)
#lines(total_r, data_sort[,3])
plot(total_r, data_sort[,4], main=header_name[4], ylim = c(0,1))
abline(h=0.5)
#lines(total_r, data_sort[,4])
plot(total_r, data_sort[,6], main=header_name[6], ylim = c(0,1))
abline(h=0.5)
#lines(total_r, data_sort[,6])
plot(total_r, data_sort[,7], main=header_name[7], ylim = c(0,1))
abline(h=0.5)
#lines(total_r, data_sort[,7])
dev.off()




