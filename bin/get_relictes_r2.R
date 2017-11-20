library(LSD)
### get parameters
args = commandArgs(trailingOnly=TRUE)

input_file_list = args[1]
input_folder = args[2]
output_folder = args[3]

ignore_sig = as.numeric(args[4])
sampling_num = as.numeric(args[5])
random_seed = as.numeric(args[6])

input_file_names = read.table(input_file_list, header = FALSE, sep='\t')

test_MSE = function(sig1, sig2){
	MSE =  sum((sig1 - sig2)^2)/length(sig1)
	#r2 =  1 - sum((sig1 - sig2)^2)/sum((sig1 - mean(sig1))^2)
	return(MSE)
}

r2_vector = c()
r_vector = c()
for ( i in seq(dim(input_file_names)[1]) ){
	sig1_file=toString(input_file_names[i,1])
	sig2_file=toString(input_file_names[i,2])
	print(paste(input_folder, sig1_file, sep=''))
	print(paste(input_folder, sig2_file, sep=''))
	### read signal file
	print('read signal files')
	sig1=read.table(paste(input_folder, sig1_file, sep=''), header=FALSE)
	sig2=read.table(paste(input_folder, sig2_file, sep=''), header=FALSE)
	### random sample
	print('random sample')
	set.seed(random_seed)
	used_id = sample(dim(sig1)[1],sampling_num)
	sig1=sig1[used_id,]
	sig2=sig2[used_id,]
	### nonzero
	used_id_nonzero = as.logical((sig1>ignore_sig) * (sig2>ignore_sig))
	sig1 = (sig1[used_id_nonzero])
	sig2 = (sig2[used_id_nonzero])

	### calculate r2 & r
	print('calculate MSE & R')
	print('MSE: ')
	r2 = test_MSE(sig1, sig2)
	print(r2)
	r2_vector[i] = r2
	print('R: ')
	r = cor(sig1, sig2)
	print(r)
	r_vector[i] = r
	### plot scatterplot
	print('plot scatterplot')
	png(paste(output_folder, sig1_file, sig2_file, '.png', sep=''))
	heatscatter(sig1, sig2, log='xy', pch = 20)
	abline(0, 1, col='red')
	dev.off()
}

### write output
print('write output')
r2r_matrix = cbind(r2_vector, r_vector)
write.table(r2r_matrix, paste(input_file_list, '.r2_r.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


