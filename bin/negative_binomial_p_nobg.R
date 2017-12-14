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
input = read.table(input_track_file, header = F)
bg_bins = read.table(bg_bins_file, header = F)

### get sig bg regions
#sig_bg = sig[bg_bins[,1]==1,]
sig_bg = sig[,1]
sig_bg_non0 = sig_bg[sig_bg>0]
sig_bg_mean = mean(sig_bg_non0)
sig_bg_var = var(sig_bg_non0)
print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_bg_var/sig_bg_mean, digits=3)) ))
print(sig_bg_mean)
print(sig_bg_var)
print(length(sig_bg_non0))

### get negative binomial parameters from signal track bg regions
sig_bg_prob = sig_bg_mean / sig_bg_var
if (sig_bg_prob<0.1){
	sig_bg_prob = 0.1
}
sig_bg_size = sig_bg_mean * sig_bg_prob / (1-sig_bg_prob)
### get input bg regions
#input_bg = input[bg_bins[,1]==1,]
input_bg = input[,1]
input_bg_non0 = input_bg[input_bg>0]
input_bg_mean = mean(input_bg_non0)
inpy_bg_var = var(input_bg_non0)
print(paste('check input track overdispersion in background regions, var/mean=', toString(round(inpy_bg_var/input_bg_mean, digits=3)) ))
print(sig_bg_prob)
print(sig_bg_size)
print(length(input_bg_non0))

print(head(input_bg))
print(summary(input_bg))
print(input_bg_mean)
print(inpy_bg_var)

### get negative binomial p-value
sig_input = cbind(sig, input)
nb_pval = apply(sig_input, MARGIN=1, function(x) 1-pnbinom(x[1], sig_bg_size, sig_bg_prob) )
### get -log10(p-value)
nb_pval[nb_pval==0] = 0.1^16

############### second round
### get sig bg regions
#sig_bg = sig[bg_bins[,1]==1,]
sig_bg = sig[nb_pval>=0.001,]
sig_bg_non0 = sig_bg[sig_bg>0]
sig_bg_mean = mean(sig_bg_non0)
sig_bg_var = var(sig_bg_non0)
print(paste('check signal track overdispersion in background regions, var/mean=', toString(round(sig_bg_var/sig_bg_mean, digits=3)) ))
print(sig_bg_mean)
print(sig_bg_var)
print(length(sig_bg_non0))

### get negative binomial parameters from signal track bg regions
sig_bg_prob = sig_bg_mean / sig_bg_var
if (sig_bg_prob<0.1){
	sig_bg_prob = 0.1
}
sig_bg_size = sig_bg_mean * sig_bg_prob / (1-sig_bg_prob)
### get input bg regions
input_bg = input[bg_bins[,1]==1,]
input_bg_non0 = input_bg[input_bg>0]
input_bg_mean = mean(input_bg_non0)
inpy_bg_var = var(input_bg_non0)
print(paste('check input track overdispersion in background regions, var/mean=', toString(round(inpy_bg_var/input_bg_mean, digits=3)) ))
print(sig_bg_prob)
print(sig_bg_size)
print(length(input_bg_non0))

print(head(input_bg))
print(summary(input_bg))
print(input_bg_mean)
print(inpy_bg_var)

### get negative binomial p-value
sig_input = cbind(sig, input)
nb_pval = apply(sig_input, MARGIN=1, function(x) 1-pnbinom(x[1], sig_bg_size * (x[2]+1)/(input_bg_mean+1), sig_bg_prob) )
### get -log10(p-value)
nb_pval[nb_pval==0] = 0.1^16
neglog10_nb_pval = -log10(nb_pval)

### write output
write.table(neglog10_nb_pval, output_name, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

