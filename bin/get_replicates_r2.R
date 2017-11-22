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
	return(MSE)
}

test_r2 = function(sig1, sig2){
	r2 =  1 - sum((sig1 - sig2)^2)/sum((sig1 - mean(sig1))^2)
	return(MSE)
}

mse_vector = c()
mse_log_vector = c()
r2_vector = c()
r2_log_vector = c()
r_vector = c()
total_mean = c()
total_mean_log2 = c()
r_vector_log2 = c()
r_name = c()

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

	### log scale
	sig1_log = log2(sig1)
	sig2_log = log2(sig2)

	### calculate r2 & r
	r_name[i] = toString(input_file_names[i,1])
	print('calculate MSE & R')

	print('MSE: ')
	mse = test_MSE(sig1, sig2)
	r2 = test_r2(sig1, sig2)
	print(r2)
	mse_vector[i] = mse
	r2_vector[i] = r2

	print('MSE (log2): ')
	mse_log_vector = test_MSE(sig1_log, sig2_log)
	r2_log2 = test_r2(sig1_log, sig2_log)
	print(r2_log2)
	mse_log_vector[i] = mse_log_vector
	r2_log_vector[i] = r2_log2

	print('R: ')
	r = cor(sig1, sig2)
	print(r)
	r_vector[i] = r

	print('R (log2): ')
	r_log2 = cor(sig1_log, sig2_log)
	print(r_log2)
	r_vector_log2[i] = r_log2

	### get total mean
	total_mean[i] = (mean(sig1) + mean(sig2))/2
	total_mean_log2[i] = (mean(sig1_log) + mean(sig2_log))/2

	### plot scatterplot
	print('plot scatterplot')
	png(paste(output_folder, sig1_file, sig2_file, '.png', sep=''))
	heatscatter(sig1, sig2, log='xy', pch = 20)
	abline(0, 1, col='red')
	dev.off()
}

### write output
print('write output')
r2r_matrix = cbind(r_name, mse_vector, r2_vector, r_vector, total_mean, mse_log_vector, r2_log_vector, r_vector_log2, total_mean_log2)
write.table(r2r_matrix, paste(input_file_list, '.r2_r.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


