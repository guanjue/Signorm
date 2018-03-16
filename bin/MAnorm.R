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
png(paste(output,".MAplot_before_rescaling.png", sep=''))
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
ma.plot(A,M,cex=1,main="MA plot before rescaling (common peaks)")
abline(b,col="green")
dev.off()

cat("M = b[1] + b[2] * A\n")
log2_peak_count_read1 <- log2(common_peak_count_read1 + small_num)
log2_peak_count_read2 <- log2(common_peak_count_read2 + small_num)
log2_peak_count_read2_rescaled <- (2-b[2])*log2_peak_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);
M_rescaled <- (log2_peak_count_read2_rescaled - log2_peak_count_read1);
A_rescaled <- (log2_peak_count_read2_rescaled + log2_peak_count_read1)/2;

png(paste(output,".MAplot_after_rescaling.png", sep=''))
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
ma.plot(as.matrix(A_rescaled),as.matrix(M_rescaled),cex=1,main=" MA plot after rescaling (all peaks)")
dev.off()


log2_allregion_count_read1 <- log2(sig1 + small_num)
log2_allregion_count_read2 <- log2(sig2 + small_num)
log2_allregion_count_read2_rescaled <- (2-b[2])*log2_allregion_count_read2/(2+b[2]) - 2*b[1]/(2+b[2]);

sig2_rescale = 2^log2_allregion_count_read2_rescaled - small_num

set.seed(2018)
sample_id = sample(nrow(sig2_rescale), random_sample_num)
png(paste(output,".scatterplot_before_rescaling.png", sep=''))
plot(log2_allregion_count_read1[sample_id], log2_allregion_count_read2[sample_id], col = 'blue')
points(log2_allregion_count_read1[sample_id][peak_binary[sample_id]], log2_allregion_count_read2[sample_id][peak_binary[sample_id]], col='orange')
points(log2_allregion_count_read1[sample_id][bg_binary[sample_id]], log2_allregion_count_read2[sample_id][bg_binary[sample_id]], col='gray')
dev.off()

png(paste(output,".scatterplot_after_rescaling.png", sep=''))
plot(log2_allregion_count_read1[sample_id], log2_allregion_count_read2_rescaled[sample_id], col = 'blue')
points(log2_allregion_count_read1[sample_id][peak_binary[sample_id]], log2_allregion_count_read2_rescaled[sample_id][peak_binary[sample_id]], col='orange')
points(log2_allregion_count_read1[sample_id][bg_binary[sample_id]], log2_allregion_count_read2_rescaled[sample_id][bg_binary[sample_id]], col='gray')
dev.off()


sig2_rescale_FRiP = sum(sig2_rescale[peak_binary]) / sum(sig2_rescale)
sig1_FRiP = sum(sig1[peak_binary]) / sum(sig1)
sig2_FRiP = sum(sig2[peak_binary]) / sum(sig2)

info = rbind(c(sum(sig2)/sum(sig1), b[1], b[2]), c(sig1_FRiP, sig2_rescale_FRiP, sig2_FRiP))

write.table(info, paste(output, '.MA.norm.info.txt', sep=''))
write.table(sig1_rescale, paste(output,".MAnorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

