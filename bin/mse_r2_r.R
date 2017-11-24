### get parameters
args = commandArgs(trailingOnly=TRUE)
signorm_sig_matrix = args[1]
totalmean_sig_matrix = args[2]
cor_list_file = args[3]
output_table_file = args[4]

##############
test_MSE = function(sig1, sig2){
	MSE =  mean((sig1 - sig2)^2)
	return(MSE)
}

test_r2 = function(sig1, sig2){
	r2 =  1 - sum((sig1 - sig2)^2)/sum((sig1 - mean(sig1))^2)
	return(r2)
}
##############
### read input
data1 = read.table(signorm_sig_matrix, header = F)
data1_sig = t(apply(data1[,c(-1:-6)], 1, as.numeric))
data1_sig = data1_sig[!is.na(rowSums(data1_sig)),]

data2 = read.table(totalmean_sig_matrix, header = F)
data2_sig = t(apply(data2[,c(-1:-6)], 1, as.numeric))
data2_sig = data2_sig[!is.na(rowSums(data2_sig)),]

### input column number
cor_list = read.table(cor_list_file, header = F)
cell_type = cor_list[,1]

### calculate MSE, R2, R
r2_vector = c()
for (i in c(1: dim(cor_list)[1])){
	print(i)
	r2 = test_r2(data1_sig[,cor_list[i,2]], data1_sig[,cor_list[i,3]])
	r2_vector[i] = r2
}

r2_vector_totalmean = c()
for (i in c(1: dim(cor_list)[1])){
	print(i)
	r2 = test_r2(data2_sig[,cor_list[i,2]], data2_sig[,cor_list[i,3]])
	r2_vector_totalmean[i] = r2
}

r_vector = c()
for (i in c(1: dim(cor_list)[1])){
	print(i)
	r = cor(data1_sig[,cor_list[i,2]], data1_sig[,cor_list[i,3]])
	r_vector[i] = r
}

r_vector_totalmean = c()
for (i in c(1: dim(cor_list)[1])){
	print(i)
	r = cor(data2_sig[,cor_list[i,2]], data2_sig[,cor_list[i,3]])
	r_vector_totalmean[i] = r
}

mse_vector = c()
for (i in c(1: dim(cor_list)[1])){
	print(i)
	mse = test_MSE(data1_sig[,cor_list[i,2]], data1_sig[,cor_list[i,3]])
	mse_vector[i] = mse
}

mse_vector_totalmean = c()
for (i in c(1: dim(cor_list)[1])){
	print(i)
	mse = test_MSE(data2_sig[,cor_list[i,2]], data2_sig[,cor_list[i,3]])
	mse_vector_totalmean[i] = mse
}

### calculate difference
r2_vector
r2_vector_totalmean
r2_dif = r2_vector - r2_vector_totalmean

r_vector
r_vector_totalmean
r_dif = r_vector - r_vector_totalmean

mse_vector
mse_vector_totalmean
mse_dif = -(mse_vector - mse_vector_totalmean)

result = data.frame(cell_type, mse_dif, r2_dif, r_dif, mse_vector, mse_vector_totalmean, r2_vector, r2_vector_totalmean, r_vector, r_vector_totalmean)

### write output
write.table(result, output_table_file, quote=F, sep='\t', row.names = FALSE, col.names = FALSE)



