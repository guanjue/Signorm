### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_filename = args[1]
output_filename = args[2]
threshold = as.numeric(args[3])
### read input signal
signal = read.table(signal_filename, header = FALSE)[,1]+1
lamda = mean(signal)
### get negative log10 p-value 
p1 = ppois(signal, lambda=lamda, lower=FALSE)

### get background signals
signal_bg = signal[p1>=threshold]
lamda_bg = mean(signal_bg)


### get negative log10 p-value 
neg_log10_p = -log10(ppois(signal, lambda=lamda_bg, lower=FALSE)+0.01)

### write signal vector
write.table(neg_log10_p, output_filename, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


