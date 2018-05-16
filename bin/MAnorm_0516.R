#library(robustreg)
library(MASS)
library(affy)
#library(R.basic) # source("http://www.braju.com/R/hbLite.R") 
		 # hbLite("R.basic")
cat('test\n')
#common_peak_count_read1=read.table(input_sig1,header=FALSE)
#common_peak_count_read2=read.table(input_sig2,header=FALSE)
### get parameters
args = commandArgs(trailingOnly=TRUE)
input_sig1 = args[1]
input_sig2 = args[2]
input_sig3 = args[3]

input_sig1_p = args[4]
input_sig2_p = args[5]

output = args[6]

small_num = as.numeric(args[7])
random_sample_num = args[8]


ct = unlist(strsplit(input_sig3, "[.]"))[1]

sig1 = scan(input_sig1)
sig2 = scan(input_sig2)
sig3 = scan(input_sig3)
totalmean_sf = sum(sig1) / sum(sig2)
#sig4 = (sig2+small_num) * totalmean_sf #- small_num
sig4 = (sig2) * totalmean_sf
#sig4[sig4>upperlim] = upperlim
#sig4[sig4<lowerlim] = lowerlim

sig1_p = scan(input_sig1_p)
sig2_p = scan(input_sig2_p)
sig1_binary = p.adjust(sig1_p, method='fdr') <= 0.05
sig2_binary = p.adjust(sig2_p, method='fdr') <= 0.05

peak_binary_pk = as.logical(sig1_binary * sig2_binary) 
peak_binary = peak_binary_pk 

bg_binary_bg = as.logical((sig1_binary + sig2_binary)==0)
bg_binary = bg_binary_bg 

#plot_binary = (sig1 != sig1[1]) & (sig2 != sig2[1]) 

common_peak_count_read1 = sig1[peak_binary]
common_peak_count_read2 = sig2[peak_binary]


M=log2((common_peak_count_read2+small_num)/(common_peak_count_read1+small_num))
A=0.5*log2((common_peak_count_read2+small_num)*(common_peak_count_read1+small_num))
M = as.matrix(M)
A = as.matrix(A)

linear=lm(M~A)$coefficients
#b=lm(M~A)$coefficients
#b=robustRegBS(M,A,beta=linear)
b=rlm(M~A)$coefficients

cat("M = b[1] + b[2] * A\n")
log2_peak_count_read1 = log2(common_peak_count_read1 + small_num)
log2_peak_count_read2 = log2(common_peak_count_read2 + small_num)
log2_peak_count_read2_rescaled = (2-b[2])*log2_peak_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);
M_rescaled = (log2_peak_count_read2_rescaled - log2_peak_count_read1);
A_rescaled = (log2_peak_count_read2_rescaled + log2_peak_count_read1)/2;

ylim = max(c(abs(min(M)), abs(max(M)), abs(min(M_rescaled)), abs(max(M_rescaled))))

png(paste(output,".MAplot_before_rescaling.png", sep=''), width = 8, height = 8, units = 'in', res = 300)
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
plot(A,M,main="MA plot before rescaling (common peaks)", ylim=c(-ylim, ylim), pch=16, cex=1)
abline(h=0,col="red",lwd=3)
abline(b,col="dodgerblue",lwd=3, lty=2)
dev.off()


png(paste(output,".MAplot_after_rescaling.png", sep=''), width = 8, height = 8, units = 'in', res = 300)
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
plot(as.matrix(A_rescaled),as.matrix(M_rescaled),main=" MA plot after rescaling (all peaks)", ylim=c(-ylim, ylim), pch=16, cex=1)
abline(h=0,col="red",lwd=3)
abline(h=0,col="dodgerblue",lwd=3, lty=2)
dev.off()


log2_allregion_count_read1 = log2(sig1 + small_num)
log2_allregion_count_read2 = log2(sig2 + small_num)
log2_allregion_count_read2_rescaled = (2-b[2])*log2_allregion_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);
sig2_rescaled = 2^log2_allregion_count_read2_rescaled 
#sig2_rescaled[sig2_rescaled > upperlim] = upperlim
#sig2_rescaled[sig2_rescaled < lowerlim] = lowerlim

log2_allregion_count_read2_rescaled = log2(sig2_rescaled )


lims_max = max(c(max(log2_allregion_count_read1), max(log2_allregion_count_read2), max(log2_allregion_count_read2_rescaled)))
lims_min =max(c(min(log2_allregion_count_read1), min(log2_allregion_count_read2), min(log2_allregion_count_read2_rescaled)))

set.seed(2018)
sample_id = sample(length(sig2_rescaled), random_sample_num)
png(paste(output,".scatterplot_before_rescaling.png", sep=''), width = 8, height = 8, units = 'in', res = 300)
pk_points_read1 = log2_allregion_count_read1[sample_id][peak_binary[sample_id]]
pk_points_read2 = log2_allregion_count_read2[sample_id][peak_binary[sample_id]]
bg_points_read1 = log2_allregion_count_read1[sample_id][bg_binary[sample_id]]
bg_points_read2 = log2_allregion_count_read2[sample_id][bg_binary[sample_id]]

pk_mean_read1 = log2(mean(sig1[peak_binary]))
pk_mean_read2 = log2(mean(sig2[peak_binary]))
bg_mean_read1 = log2(mean(sig1[bg_binary]))
bg_mean_read2 = log2(mean(sig2[bg_binary]))
total_mean_read1 = log2(mean(sig1))
total_mean_read2 = log2(mean(sig2))

plot(log2_allregion_count_read2[sample_id], log2_allregion_count_read1[sample_id], col = 'dodgerblue', pch=16, xlim=c(lims_min, lims_max), ylim=c(lims_min, lims_max), cex=1, main='Raw_Signal', ylab='reference', xlab=ct)
points(pk_points_read2, pk_points_read1, col='darkorange1', pch=16, cex=1)
points(bg_points_read2, bg_points_read1, col='gray28', pch=16, cex=1)
points(pk_mean_read2, pk_mean_read1, col='black', pch=16, cex=2)
points(bg_mean_read2, bg_mean_read1, col='black', pch=16, cex=2)
points(total_mean_read2, total_mean_read1, col='red', pch=16, cex=2)
abline(0,1,lwd=3,col='black')
lines(c(bg_mean_read2, pk_mean_read2), c(bg_mean_read1, pk_mean_read1), col='green3', lty=2, lwd=3)
dev.off()

png(paste(output,".scatterplot_after_MAnorm.png", sep=''), width = 8, height = 8, units = 'in', res = 300)
pk_points_read1 = log2_allregion_count_read1[sample_id][peak_binary[sample_id]]
pk_points_read2_rescaled = log2_allregion_count_read2_rescaled[sample_id][peak_binary[sample_id]]
bg_points_read1 = log2_allregion_count_read1[sample_id][bg_binary[sample_id]]
bg_points_read2_rescaled = log2_allregion_count_read2_rescaled[sample_id][bg_binary[sample_id]]

pk_mean_read1 = log2(mean(sig1[peak_binary]))
pk_mean_read2_rescaled = log2(mean(sig2_rescaled[peak_binary]))
bg_mean_read1 = log2(mean(sig1[bg_binary]))
bg_mean_read2_rescaled = log2(mean(sig2_rescaled[bg_binary]))
total_mean_read1 = log2(mean(sig1))
total_mean_read2_rescaled = log2(mean(sig2_rescaled))

plot(log2_allregion_count_read2_rescaled[sample_id], log2_allregion_count_read1[sample_id], col = 'dodgerblue', pch=16, xlim=c(lims_min, lims_max), ylim=c(lims_min, lims_max), cex=1, main='MAnorm', ylab='reference', xlab=ct)
points(pk_points_read2_rescaled, pk_points_read1, col='darkorange1', pch=16, cex=1)
points(bg_points_read2_rescaled, bg_points_read1, col='gray28', pch=16, cex=1)
points(pk_mean_read2_rescaled, pk_mean_read1, col='black', pch=16, cex=2)
points(bg_mean_read2_rescaled, bg_mean_read1, col='black', pch=16, cex=2)
points(total_mean_read2_rescaled, total_mean_read1, col='red', pch=16, cex=2)
abline(0,1,lwd=3,col='black')
lines(c(bg_mean_read2_rescaled, pk_mean_read2_rescaled), c(bg_mean_read1, pk_mean_read1), col='green3', lty=2, lwd=3)
dev.off()

png(paste(output,".scatterplot_after_PKnorm.png", sep=''), width = 8, height = 8, units = 'in', res = 300)

log2_allregion_count_read1 = log2(sig1 + small_num)
log2_allregion_count_read3 = log2(sig3 + small_num)

pk_points_read1 = log2_allregion_count_read1[sample_id][peak_binary[sample_id]]
pk_points_read3 = log2_allregion_count_read3[sample_id][peak_binary[sample_id]]
bg_points_read1 = log2_allregion_count_read1[sample_id][bg_binary[sample_id]]
bg_points_read3 = log2_allregion_count_read3[sample_id][bg_binary[sample_id]]

pk_mean_read1 = log2(mean(sig1[peak_binary]))
pk_mean_read3 = log2(mean(sig3[peak_binary]))
bg_mean_read1 = log2(mean(sig1[bg_binary]))
bg_mean_read3 = log2(mean(sig3[bg_binary]))
total_mean_read1 = log2(mean(sig1))
total_mean_read3 = log2(mean(sig3))

plot(log2_allregion_count_read3[sample_id], log2_allregion_count_read1[sample_id], col = 'dodgerblue', pch=16, xlim=c(lims_min, lims_max), ylim=c(lims_min, lims_max), cex=1, main='PKnorm', ylab='reference', xlab=ct)
points(pk_points_read3, pk_points_read1, col='darkorange1', pch=16, cex=1)
points(bg_points_read3, bg_points_read1, col='gray28', pch=16, cex=1)
points(pk_mean_read3, pk_mean_read1, col='black', pch=16, cex=2)
points(bg_mean_read3, bg_mean_read1, col='black', pch=16, cex=2)
points(total_mean_read3, total_mean_read1, col='red', pch=16, cex=2)
abline(0,1,lwd=3,col='black')
lines(c(bg_mean_read3, pk_mean_read3), c(bg_mean_read1, pk_mean_read1), col='green3', lty=2, lwd=3)
dev.off()

png(paste(output,".scatterplot_after_totalmeannorm.png", sep=''), width = 8, height = 8, units = 'in', res = 300)

log2_allregion_count_read1 = log2(sig1 + small_num)
log2_allregion_count_read4 = log2(sig4 + small_num)

pk_points_read1 = log2_allregion_count_read1[sample_id][peak_binary[sample_id]]
pk_points_read4 = log2_allregion_count_read4[sample_id][peak_binary[sample_id]]
bg_points_read1 = log2_allregion_count_read1[sample_id][bg_binary[sample_id]]
bg_points_read4 = log2_allregion_count_read4[sample_id][bg_binary[sample_id]]

pk_mean_read1 = log2(mean(sig1[peak_binary]))
pk_mean_read4 = log2(mean(sig4[peak_binary]))
bg_mean_read1 = log2(mean(sig1[bg_binary]))
bg_mean_read4 = log2(mean(sig4[bg_binary]))
total_mean_read1 = log2(mean(sig1))
total_mean_read4 = log2(mean(sig4))

plot(log2_allregion_count_read4[sample_id], log2_allregion_count_read1[sample_id], col = 'dodgerblue', pch=16, xlim=c(lims_min, lims_max), ylim=c(lims_min, lims_max), cex=1, main='Total_Signal', ylab='reference', xlab=ct)
points(pk_points_read4, pk_points_read1, col='darkorange1', pch=16, cex=1)
points(bg_points_read4, bg_points_read1, col='gray28', pch=16, cex=1)
points(pk_mean_read4, pk_mean_read1, col='black', pch=16, cex=2)
points(bg_mean_read4, bg_mean_read1, col='black', pch=16, cex=2)
points(total_mean_read4, total_mean_read1, col='red', pch=16, cex=2)
abline(0,1,lwd=3,col='black')
lines(c(bg_mean_read4, pk_mean_read4), c(bg_mean_read1, pk_mean_read1), col='green3', lty=2, lwd=3)
dev.off()



sig2_rescaled_FRiP = sum(sig2_rescaled[sig2_binary]) / sum(sig2_rescaled)
sig1_FRiP = sum(sig1[sig1_binary]) / sum(sig1)
sig2_FRiP = sum(sig2[sig2_binary]) / sum(sig2)
sig3_FRiP = sum(sig3[sig2_binary]) / sum(sig3)
sig4_FRiP = sum(sig4[sig2_binary]) / sum(sig4)

info = rbind(c(sum(sig2), sum(sig1), b[1], b[2], totalmean_sf), c(sig1_FRiP, sig2_FRiP, sig2_rescaled_FRiP, sig3_FRiP, sig4_FRiP))

write.table(info, paste(output, '.MA.norm.info.txt', sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(sig2_rescaled, paste(output,".MAnorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(sig4, paste(output,".totalsig_norm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

mvfile2folder <- function(from, to) {
    todir <- dirname(to)
    if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
    file.rename(from = from,  to = to)
}

dir.create(paste(output, '_totalmean_MAnorm', sep=''), showWarnings = FALSE)
mvfile2folder(from=paste(output, '.MA.norm.info.txt', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.MA.norm.info.txt', sep=''))
mvfile2folder(from=paste(output,".MAnorm.txt", sep=''), to=paste(output, '_totalmean_MAnorm/', output,".MAnorm.txt", sep=''))
mvfile2folder(from=paste(output, '.totalsig_norm.txt', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.totalsig_norm.txt', sep=''))
mvfile2folder(from=paste(output, '.scatterplot_after_totalmeannorm.png', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.scatterplot_after_totalmeannorm.png', sep=''))
mvfile2folder(from=paste(output, '.scatterplot_after_PKnorm.png', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.scatterplot_after_PKnorm.png', sep=''))
mvfile2folder(from=paste(output, '.scatterplot_after_MAnorm.png', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.scatterplot_after_MAnorm.png', sep=''))
mvfile2folder(from=paste(output, '.scatterplot_before_rescaling.png', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.scatterplot_before_rescaling.png', sep=''))
mvfile2folder(from=paste(output, '.MAplot_after_rescaling.png', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.MAplot_after_rescaling.png', sep=''))
mvfile2folder(from=paste(output, '.MAplot_before_rescaling.png', sep=''), to=paste(output, '_totalmean_MAnorm/', output, '.MAplot_before_rescaling.png', sep=''))

