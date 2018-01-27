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
sig1_fsip = c()
sig2_norm_fsip = c()
sig2_fsip = c()

### cell type
cell_type = c()

for (i in c(1:dim(info_list_files)[1])){
	print(info_list_files[i,1])
	info = read.table(toString(info_list_files[i,1]), header = F)
	### sf vector
	totalreads_sf[i] = info[1,1]
	pk_sf[i] = info[1,2]
	bg_sf[i] = info[1,3]
	### FSiP vector
	sig1_fsip[i] = info[2,1]
	sig2_norm_fsip[i] = info[2,2]
	sig2_fsip[i] = info[2,3]
	### cell type
	cell_type[i] = unlist(strsplit(toString(info_list_files[i,1]), "[.]"))[1]
}


######### plot scale factor plot
ymin_sf = min(c(totalreads_sf, pk_sf, bg_sf))
ymax_sf = max(c(totalreads_sf, pk_sf, bg_sf))

pdf(paste(outputname, '.sf.pdf', sep=''))
plot(totalreads_sf, pch=20, col='black', ylim=c(ymin_sf, ymax_sf), axes=FALSE, xlab='')
axis(2)
axis(1, at=c(1: length(cell_type)),labels=cell_type, las=2)
points(pk_sf, pch=20, col='orange')
points(bg_sf, pch=20, col='blue')
for (i in c(1:length(totalreads_sf))){
	x_tmp = i
	ymin = min(c(totalreads_sf[i], pk_sf[i], bg_sf[i]))
	ymax = max(c(totalreads_sf[i], pk_sf[i], bg_sf[i]))
	segments(x_tmp, ymin, x_tmp, ymax, col='gray', lwd=1.5)
}
dev.off()

ymin_FSiP = min(c(sig2_fsip, sig2_norm_fsip, sig1_fsip))-0.1
ymax_FSiP = max(c(sig2_fsip, sig2_norm_fsip, sig1_fsip))+0.1

pdf(paste(outputname, '.FSiP.pdf', sep=''))
plot(sig2_fsip, pch=20, col='black', ylim=c(ymin_FSiP, ymax_FSiP), axes=FALSE, xlab='')
axis(2)
axis(1, at=c(1: length(cell_type)),labels=cell_type, las=2)
points(sig2_norm_fsip, pch=20, col='orange')
lines(sig1_fsip, pch=20, col='blue')
for (i in c(1:length(sig2_norm_fsip))){
	x_tmp = i
	ymin = min(c(sig2_fsip[i], sig2_norm_fsip[i], sig1_fsip[i]))
	ymax = max(c(sig2_fsip[i], sig2_norm_fsip[i], sig1_fsip[i]))
	segments(x_tmp, ymin, x_tmp, ymax, col='gray', lwd=1.5)
}
dev.off()


######### plot FSiP factor plot
FSiP = cbind(sig2_fsip, sig2_norm_fsip)
pdf(paste(outputname, '.FSiP.box.pdf', sep=''))
boxplot(FSiP)#, ylim=c(0,0.8))
dev.off()








