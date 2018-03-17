### get parameters
args = commandArgs(trailingOnly=TRUE)
info_list = args[1]
outputname = args[2]

info_list_files = read.table(info_list, header = F)

### sf vectors
pk_sf = c()
bg_sf = c()
totalreads_sf = c()

### FSiP vectors
sig1_ref_fsip = c()
sig2_raw_fsip = c()
sig2_manorm_fsip = c()
sig2_pknorm_fsip = c()
sig2_totalsignorm_fsip = c()

### cell type
cell_type = c()

for ( i in c(1:dim(info_list_files)[1])){
	print(info_list_files[i,1])
	info = read.table(toString(info_list_files[i,1]), header = F)
	### sf vector
	#totalreads_sf[i] = info[1,1]
	#pk_sf[i] = info[1,2]
	#bg_sf[i] = info[1,3]
	### FSiP vector
	sig1_ref_fsip[i] = info[2,1]
	sig2_raw_fsip[i] = info[2,2]
	sig2_manorm_fsip[i] = info[2,3]
	sig2_pknorm_fsip[i] = info[2,4]
	sig2_totalsignorm_fsip[i] = info[2,5]
	### cell type
	cell_type[i] = unlist(strsplit(toString(info_list_files[i,1]), "[rep]"))[1]
}


ymin_FSiP = min(c(sig1_ref_fsip, sig2_raw_fsip, sig2_manorm_fsip, sig2_pknorm_fsip, sig2_totalsignorm_fsip))-0.2
ymax_FSiP = max(c(sig1_ref_fsip, sig2_raw_fsip, sig2_manorm_fsip, sig2_pknorm_fsip, sig2_totalsignorm_fsip))+0.2

pdf(paste(outputname, '.FSiP.pdf', sep=''))
plot(sig2_raw_fsip, pch=20, col='black', ylim=c(ymin_FSiP, ymax_FSiP), axes=FALSE, xlab='')
axis(2)
axis(1, at=c(1: length(cell_type)),labels=cell_type, las=2)
points(sig2_manorm_fsip, pch=20, col='orange')
points(sig2_pknorm_fsip, pch=20, col='green')
points(sig2_totalsignorm_fsip, pch=20, col='red')
lines(sig1_ref_fsip, pch=20, col='black')
for (i in c(1:length(sig2_raw_fsip))){
	x_tmp = i
	ymin = min(c(sig2_raw_fsip[i], sig2_manorm_fsip[i], sig2_pknorm_fsip[i], sig2_totalsignorm_fsip[i]))
	ymax = max(c(sig2_raw_fsip[i], sig2_manorm_fsip[i], sig2_pknorm_fsip[i], sig2_totalsignorm_fsip[i]))
	segments(x_tmp, ymin, x_tmp, ymax, col='gray', lwd=1.5)
}
dev.off()


######### plot FSiP factor plot
FSiP = cbind(sig2_raw_fsip, sig2_totalsignorm_fsip, sig2_manorm_fsip, sig2_pknorm_fsip)
colnames(FSiP) = c('raw_signal', 'totalsignorm', 'MAnorm', 'PKnorm')
pdf(paste(outputname, '.FSiP.box.pdf', sep=''))
boxplot(FSiP, ylim=c(ymin_FSiP, ymax_FSiP))
abline(h=sig1_ref_fsip, pch=20, col='blue', lwd=1.5, lty=2)
dev.off()








