library(hash)
library(LSD)
library(changepoint)
### get parameters
args = commandArgs(trailingOnly=TRUE)

signal_track_file = args[1]
input_track_file = args[2]
bg_bins_file = args[3]
output_name = args[4]

### read data
sig = read.table(signal_track_file, header = F)
input = read.table(signal_track_file, header = F)
bg_bins = read.table(bg_bins_file, header = F)

### get sig bg regions
sig_bg = sig[bg_bins[,1]==1,]
sig_bg_mean = mean(sig_bg)
sig_bg_var = var(sig_bg)
print(paste('check signal track overdispersion in background regions, var/mean=', toString(sig_bg_var/sig_bg_mean, digits=3))))

### get negative binomial parameters from signal track bg regions
sig_bg_prob = sig_bg_mean / sig_bg_var
sig_bg_size = sig_bg_mean * sig_bg_prob / (1-sig_bg_prob)

### get input bg regions
input_bg = input[bg_bins[,1]==1,]
input_bg_mean = mean(input_bg+1)
inpy_bg_var = var(input_bg)
print(paste('check input track overdispersion in background regions, var/mean=', toString(inpy_bg_var/input_bg_mean, digits=3))))

### get negative binomial p-value
sig_input = cbind(sig, input)
nb_pval = apply(sig_input, MARGIN=1, function(x) 1-pnbinom(x[1], sig_bg_size * (x[2]+1)/input_bg_mean, sig_bg_prob) )
### get -log10(p-value)
neglog10_nb_pval = -log10(nb_pval)

### write output
write.table(neglog10_nb_pval, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

