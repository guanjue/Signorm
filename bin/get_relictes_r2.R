library(LSD)
### get parameters
args = commandArgs(trailingOnly=TRUE)

input_file_list = args[1]
output_folder = args[2]
input_file_names = read.table('info_table_compare_r2_signorm.txt', header = FALSE, sep='\t')

test_r2 = function(sig1, sig2){
	r2 =  1 - sum((sig1[,i] - sig2[,i])^2)/sum((sig1[,i] - mean(sig1[,i]))^2)
	return(r2)
}

for ( i in seq(dim(input_file_names[i,1])) ){
	sig1_file=toString(input_file_names[i,1])
	sig2_file=toString(input_file_names[i,2])
	print(sig1_file)
	print(sig2_file)
	### read signal file
	sig1=read.table(sig1_file, header=FALSE)
	sig2=read.table(sig2_file, header=FALSE)
	### calculate r2 & r
	r2 = test_r2(sig1, sig2)
	print(r2)
	r = cor(sig1, sig2)
	print(r)
	### plot scatterplot
	png(paste(output_folder, sig1_file, sig2_file, '.png', sep=''))
	heatscatter(sig1, sig2, log='xy', pch = 20)
	abline(0, 1, col='red')
	dev.off()
}




