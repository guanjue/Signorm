#library(robustreg)
library(MASS)
library(affy)
#library(R.basic) # source("http://www.braju.com/R/hbLite.R") 
		 # hbLite("R.basic")
cat('test\n')
#common_peak_count_read1<-read.table(input_sig1,header=FALSE)
#common_peak_count_read2<-read.table(input_sig2,header=FALSE)

sig1 = read.table(input_sig1,header=FALSE)
sig2 = read.table(input_sig2,header=FALSE)

sig1_binary = 10**(-sig1) <= 0.01
sig2_binary = 10**(-sig2) <= 0.01

peak_binary_pk = as.logical(sig1_binary[:,0] * sig2_binary[:,0]) 
peak_binary = peak_binary_pk & (sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]) 

bg_binary_bg = as.logical((sig1_binary[:,0] + sig2_binary[:,0])==0)
bg_binary = bg_binary_bg & (sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]) 

common_peak_count_read1 = sig1[peak_binary,0]+small_num
common_peak_count_read2 = sig2[peak_binary,0]+small_num


M<-log2((common_peak_count_read1+1)/(common_peak_count_read2+1))
A<-0.5*log2((common_peak_count_read1+1)*(common_peak_count_read2+1))
M <- as.matrix(M)
A <- as.matrix(A)

linear<-lm(M~A)$coefficients
#b<-lm(M~A)$coefficients
#b<-robustRegBS(M,A,beta=linear)
b<-rlm(M~A)$coefficients
png('MAplot_before_rescaling.png')
#ma.plot(A,M,cex=1,main=paste(dataname," MA plot before rescaling (common peaks)",sep=""))
ma.plot(A,M,cex=1,main="MA plot before rescaling (common peaks)")
abline(b,col="green")
dev.off()

cat("M = b[1] + b[2] * A\n")
log2_peak_count_read1 <- log2(common_peak_count_read1 + 1)
log2_peak_count_read2 <- log2(common_peak_count_read2 + 1)
log2_peak_count_read1_rescaled <- (2-b[2])*log2_peak_count_read1/(2+b[2]) - 2*b[1]/(2+b[2]);
M_rescaled <- (log2_peak_count_read1_rescaled - log2_peak_count_read2);
A_rescaled <- (log2_peak_count_read1_rescaled + log2_peak_count_read2)/2;

png('MAplot_after_rescaling.png')
#ma.plot(A_rescaled,M_rescaled,cex=1,main=paste(dataname," MA plot after rescaling (all peaks)",sep=""))
ma.plot(as.matrix(A_rescaled),as.matrix(M_rescaled),cex=1,main=" MA plot after rescaling (all peaks)")
dev.off ()


log2_allregion_count_read1 <- log2(sig1 + 1)
log2_allregion_count_read2 <- log2(sig2 + 1)
log2_allregion_count_read1_rescaled <- (2-b[2])*log2_allregion_count_read1/(2+b[2]) - 2*b[1]/(2+b[2]);

sig1_rescale = 2^log2_allregion_count_read1_rescaled - 1

write.table(sig1_rescale, paste(output,".MAnorm.txt", sep=''),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

