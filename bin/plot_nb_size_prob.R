### get parameters
args = commandArgs(trailingOnly=TRUE)

cell_marker = args[1]
tail = args[2]
input_folder = args[3]

### extract filenames of the cell marker
file_list = list.files(input_folder, pattern=paste('^', '(.*)', cell_marker, '(.*)', tail, '$', sep='') )
cell_type = apply(as.matrix(file_list), 1, function(x) unlist(strsplit(x, "\\."))[1])
### read files of the cell marker
mean_matrix = NULL
for (file in file_list){
	d = read.table(file, header = F)
	mean_matrix = cbind(mean_matrix, d[,1])
}

var_matrix = NULL
for (file in file_list){
	d = read.table(file, header = F)
	var_matrix = cbind(var_matrix, d[,2])
}

size_matrix = NULL
for (file in file_list){
	d = read.table(file, header = F)
	size_matrix = cbind(size_matrix, d[,3])
}

prob_matrix = NULL
for (file in file_list){
	d = read.table(file, header = F)
	prob_matrix = cbind(prob_matrix, d[,4])
}

index = seq(1:dim(mean_matrix)[2])
###
png(paste(cell_marker, '.nb_para.png', sep=''), width = 1000, height = 1000)
par(mfrow=c(2,2))
### means
plot(index, mean_matrix[1,], col='black', ylim=c(0, max(mean_matrix)), type='l', xaxt="n")
axis(1, at=index,labels=cell_type, col.axis="black", las=2)
points(index, mean_matrix[1,], col='black')
lines(index, mean_matrix[2,], col='blue')
points(index, mean_matrix[2,], col='blue')
lines(index, mean_matrix[3,], col='red')
points(index, mean_matrix[3,], col='red')
### var
plot(index, var_matrix[1,], col='black', ylim=c(0, max(var_matrix)), type='l', xaxt="n")
axis(1, at=index,labels=cell_type, col.axis="black", las=2)
points(index, var_matrix[1,], col='black')
lines(index, var_matrix[2,], col='blue')
points(index, var_matrix[2,], col='blue')
lines(index, var_matrix[3,], col='red')
points(index, var_matrix[3,], col='red')
### size
plot(index, size_matrix[1,], col='black', ylim=c(0, max(size_matrix)), type='l', xaxt="n")
axis(1, at=index,labels=cell_type, col.axis="black", las=2)
points(index, size_matrix[1,], col='black')
lines(index, size_matrix[2,], col='blue')
points(index, size_matrix[2,], col='blue')
lines(index, size_matrix[3,], col='red')
points(index, size_matrix[3,], col='red')
### prob
plot(index, prob_matrix[1,], col='black', ylim=c(0, max(prob_matrix)), type='l', xaxt="n")
axis(1, at=index,labels=cell_type, col.axis="black", las=2)
points(index, prob_matrix[1,], col='black')
lines(index, prob_matrix[2,], col='blue')
points(index, prob_matrix[2,], col='blue')
lines(index, prob_matrix[3,], col='red')
points(index, prob_matrix[3,], col='red')
dev.off()


