### get parameters
args = commandArgs(trailingOnly=TRUE)
input_filename = args[1]
output_filename = args[2]
round_digit = as.numeric(args[3])

### read input signal
signal = read.table(input_filename, header = FALSE)

### round signal
signal_round = round(a, digit = round_digit)

### write signal vector
write.table(signal_round, output_filename, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
