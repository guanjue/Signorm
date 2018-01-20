### get parameters
args = commandArgs(trailingOnly=TRUE)
input_matrix = args[1]
output_filename = args[2]

data = read.table(input_matrix, header = F)

### get row 1 / variance-to-mean ratio (1/VMR) 
data_vmr = apply(data, 1, function(x)  (abs(mean(x))+0.001) / (var(x)+0.001) )

### get output
write.table(data_vmr, output_filename, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


