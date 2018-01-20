### get parameters
args = commandArgs(trailingOnly=TRUE)
input_matrix = args[1]
output_filename = args[2]

data = read.table(input_matrix, header = F)

### get row 1 / Coefficient of Variation (1/CV)
data_sv = apply(data, 1, function(x)  mean(x) / sd(x) )

### get output
write.table(data_sv, output_filename, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


