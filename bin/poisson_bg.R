### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_filename = args[1]
background_filename = args[2]
output_filename = args[3]

### read input signal
signal = read.table(signal_filename, header = FALSE)[,1]+1
bg = read.table(background_filename, header = FALSE)[,1]+1

### get negative log10 p-value 
neg_log10_p = -log10(dpois(signal, lambda=bg))

### write signal vector
write.table(neg_log10_p, output_filename, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


