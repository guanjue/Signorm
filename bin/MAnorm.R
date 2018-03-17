#library(robustreg)
library(MASS)
library(affy)
#library(R.basic) # source("http://www.braju.com/R/hbLite.R") 
		 # hbLite("R.basic")
cat('test\n')
#common_peak_count_read1<-read.table(input_sig1,header=FALSE)
#common_peak_count_read2<-read.table(input_sig2,header=FALSE)
### get parameters
args = commandArgs(trailingOnly=TRUE)
input_sig1 = args[1]
input_sig2 = args[2]
output = args[3]

small_num = 0.01
random_sample_num = 50000

sig1 = scan(input_sig1)
sig2 = scan(input_sig2)

sig1_binary = 10^(-sig1) <= 0.01
sig2_binary = 10^(-sig2) <= 0.01

peak_binary_pk = as.logical(sig1_binary * sig2_binary) 
peak_binary = peak_binary_pk & (sig1 != sig1[1]) & (sig2 != sig2[1])

bg_binary_bg = as.logical((sig1_binary + sig2_binary)==0)
bg_binary = bg_binary_bg & (sig1 != sig1[1]) & (sig2 != sig2[1]) 

plot_binary = (sig1 != sig1[1]) & (sig2 != sig2[1]) 

common_peak_count_read1 = sig1[peak_binary]+small_num
common_peak_count_read2 = sig2[peak_binary]+small_num


M<-log2((common_peak_count_read2+small_num)/(common_peak_count_read1+small_num))
A<-0.5*log2((common_peak_count_read2+small_num)*(common_peak_count_read1+small_num))
M <- as.matrix(M)
A <- as.matrix(A)

linear<-lm(M~A)$coefficients
#b<-lm(M~A)$coefficients
#b<-robustRegBS(M,A,beta=linear)
b<-rlm(M~A)$coefficients

cat("M = b[1] + b[2] * A\n")
log2_peak_count_read1 <- log2(common_peak_count_read1 + small_num)
log2_peak_count_read2 <- log2(common_peak_count_read2 + small_num)
log2_peak_count_read2_rescaled <- (2-b[2])*log2_peak_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);
M_rescaled <- (log2_peak_count_read2_rescaled - log2_peak_count_read1);
A_rescaled <- (log2_peak_count_read2_rescaled + log2_peak_count_read1)/2;

ylim = max(c(abs(min(M)), abs(max(M)), abs(min(M_rescaled)), abs(max(M_rescaled))))

png(paste(output,".MAplot_before_rescaling.png", sep=''))
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
plot(A,M,cex=1,main="MA plot before rescaling (common peaks)", ylim=c(-ylim, ylim))
abline(b,col="dodgerblue",lwd=3, lty=2)
abline(h=0,col="red",lwd=3)
dev.off()


png(paste(output,".MAplot_after_rescaling.png", sep=''))
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
plot(as.matrix(A_rescaled),as.matrix(M_rescaled),cex=1,main=" MA plot after rescaling (all peaks)", ylim=c(-ylim, ylim))
abline(h=0,col="dodgerblue",lwd=3, lty=2)
abline(h=0,col="red",lwd=3)
dev.off()


log2_allregion_count_read1 <- log2(sig1 + small_num)
log2_allregion_count_read2 <- log2(sig2 + small_num)
log2_allregion_count_read2_rescaled <- (2-b[2])*log2_allregion_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);

sig2_rescale = 2^log2_allregion_count_read2_rescaled - small_num

lims_max = max(c(max(log2_allregion_count_read1), max(log2_allregion_count_read2), max(log2_allregion_count_read2_rescaled)))
lims_min =max(c(min(log2_allregion_count_read1), min(log2_allregion_count_read2), min(log2_allregion_count_read2_rescaled)))

set.seed(2018)
sample_id = sample(length(sig2_rescale[plot_binary]), random_sample_num)
png(paste(output,".scatterplot_before_rescaling.png", sep=''))
pk_points_read1 = log2_allregion_count_read1[plot_binary][sample_id][peak_binary[sample_id]]
pk_points_read2 = log2_allregion_count_read2[plot_binary][sample_id][peak_binary[sample_id]]
bg_points_read1 = log2_allregion_count_read1[plot_binary][sample_id][bg_binary[sample_id]]
bg_points_read2 = log2_allregion_count_read2[plot_binary][sample_id][bg_binary[sample_id]]
plot(log2_allregion_count_read1[plot_binary][sample_id], log2_allregion_count_read2[plot_binary][sample_id], col = 'dodgerblue', pch='.', xlim=c(lims_min, lims_max), ylim=c(lims_min, lims_max), cex=3)
points(pk_points_read1, pk_points_read2, col='darkorange1', pch='.', cex=3)
points(bg_points_read1, bg_points_read2, col='gray', pch='.', cex=3)
points(mean(pk_points_read1), mean(pk_points_read2), col='black', pch='.', cex=5)
points(mean(bg_points_read1), mean(bg_points_read2), col='black', pch='.', cex=5)
lines(c(mean(bg_points_read1), mean(pk_points_read1)), c(mean(bg_points_read2), mean(pk_points_read2)), col='gray', lty=2, lwd=3)
abline(0,1,lwd=3,col='black')
dev.off()

png(paste(output,".scatterplot_after_rescaling.png", sep=''))
pk_points_read1 = log2_allregion_count_read1[plot_binary][sample_id][peak_binary[sample_id]]
pk_points_read2 = log2_allregion_count_read2_rescaled[plot_binary][sample_id][peak_binary[sample_id]]
bg_points_read1 = log2_allregion_count_read2_rescaled[plot_binary][sample_id][bg_binary[sample_id]]
bg_points_read2 = log2_allregion_count_read2_rescaled[plot_binary][sample_id][bg_binary[sample_id]]
plot(log2_allregion_count_read1[plot_binary][sample_id], log2_allregion_count_read2_rescaled[plot_binary][sample_id], col = 'dodgerblue', pch='.', xlim=c(lims_min, lims_max), ylim=c(lims_min, lims_max), cex=3)
points(pk_points_read1, pk_points_read2, col='darkorange1', pch='.', cex=3)
points(bg_points_read1, bg_points_read2, col='gray', pch='.', cex=3)
points(mean(pk_points_read1), mean(pk_points_read2), col='black', pch='.', cex=5)
points(mean(bg_points_read1), mean(bg_points_read2), col='black', pch='.', cex=5)
lines(c(mean(bg_points_read1), mean(pk_points_read1)), c(mean(bg_points_read2), mean(pk_points_read2)), col='gray', lty=2, lwd=3)
abline(0,1,lwd=3,col='black')
dev.off()


sig2_rescale_FRiP = sum(sig2_rescale[peak_binary]) / sum(sig2_rescale)
sig1_FRiP = sum(sig1[peak_binary]) / sum(sig1)
sig2_FRiP = sum(sig2[peak_binary]) / sum(sig2)

info = rbind(c(sum(sig2)/sum(sig1), b[1], b[2]), c(sig1_FRiP, sig2_rescale_FRiP, sig2_FRiP))

write.table(info, paste(output, '.MA.norm.info.txt', sep=''))
write.table(sig2_rescale, paste(output,".MAnorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

