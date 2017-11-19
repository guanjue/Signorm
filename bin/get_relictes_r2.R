### get parameters
args = commandArgs(trailingOnly=TRUE)

input_file_list = args[1]

input_file_names = read.table(input_file_list, header = FALSE)

test_r2 = function(sig1, sig2){
	r2 =  1 - sum((sig1[,i] - sig2[,i])^2)/sum((sig1[,i] - mean(sig1[,i]))^2)
	return(r2)
}

for (files in input_file_names){
	print(files[1])
	print(files[2])
	### read signal file
	sig1=read.table(files[1], header=FALSE)
	sig2=read.table(files[2], header=FALSE)
	### calculate r2 & r
	r2 = test_r2(sig1, sig2)
	print(r2)
	r = cor(sig1, sig2)
	print(r)
}




