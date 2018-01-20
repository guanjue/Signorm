### get parameters
args = commandArgs(trailingOnly=TRUE)
input_matrix = args[1]
output_filename = args[2]

data = read.table(input_matrix, header = F)

### get row Coefficient of Variation (CV)
data_sv = apply(data, 1, function(x) sd(x) / mean(x))

### get output
write.table(data_sv, output_filename, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


