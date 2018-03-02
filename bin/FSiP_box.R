### get parameters
args = commandArgs(trailingOnly=TRUE)
info_list1 = args[1]
info_list2 = args[2]
info_list3 = args[3]
info_list4 = args[4]
outputname = args[5]
ref_ct = args[6]

info_list_files1 = read.table(info_list1, header = F)
info_list_files2 = read.table(info_list2, header = F)
info_list_files3 = read.table(info_list3, header = F)
info_list_files4 = read.table(info_list4, header = F)

### FSiP vectors
sig1_fsip = c()
sig2_fsip = c()
sig3_fsip = c()
sig4_fsip = c()

### cell type
filename = c()
cell_type = c()

for (i in c(1:dim(info_list_files1)[1])){
	print(info_list_files1[i,1])
	info = read.table(toString(info_list_files1[i,1]), header = F)
	### FSiP vector
	sig1_fsip[i] = info[1,1]
	### cell type
	filename[i] = unlist(strsplit(toString(info_list_files1[i,1]), "[/]"))[2]
	cell_type[i] = unlist(strsplit(toString(filename[i]), "[.]"))[1]
	if (cell_type[i] == ref_ct){
		ref_fsip = sig1_fsip[i]
	}
}

for (i in c(1:dim(info_list_files2)[1])){
	print(info_list_files2[i,1])
	info = read.table(toString(info_list_files2[i,1]), header = F)
	### FSiP vector
	sig2_fsip[i] = info[1,1]
	### cell type
	filename[i] = unlist(strsplit(toString(info_list_files1[i,1]), "[/]"))[2]
	cell_type[i] = unlist(strsplit(toString(filename[i]), "[.]"))[1]
	if (cell_type[i] == ref_ct){
		ref_fsip = sig1_fsip[i]
	}
}

for (i in c(1:dim(info_list_files3)[1])){
	print(info_list_files3[i,1])
	info = read.table(toString(info_list_files3[i,1]), header = F)
	### FSiP vector
	sig3_fsip[i] = info[1,1]
}

for (i in c(1:dim(info_list_files4)[1])){
	print(info_list_files4[i,1])
	info = read.table(toString(info_list_files4[i,1]), header = F)
	### FSiP vector
	sig4_fsip[i] = info[1,1]
}


######### plot boxplot for FSiP
ymin_FSiP = min(c(sig1_fsip, sig2_fsip, sig3_fsip, sig4_fsip))-0.1
ymax_FSiP = max(c(sig1_fsip, sig2_fsip, sig3_fsip, sig4_fsip))+0.1

pdf(paste(outputname, '.FSiP.pdf', sep=''))
plot(sig1_fsip, pch=20, col='black', ylim=c(ymin_FSiP, ymax_FSiP), axes=FALSE, xlab='')
axis(2)
axis(1, at=c(1: length(cell_type)),labels=cell_type, las=2)
points(sig2_fsip, pch=20, col='blue')
points(sig3_fsip, pch=20, col='green')
points(sig4_fsip, pch=20, col='orange')
abline(h=ref_fsip, col='blue', lwd=1.5)

for (i in c(1:length(sig1_fsip))){
	x_tmp = i
	ymin = min(c(sig1_fsip[i], sig2_fsip[i], sig3_fsip[i], sig4_fsip[i]))
	ymax = max(c(sig1_fsip[i], sig2_fsip[i], sig3_fsip[i], sig4_fsip[i]))
	segments(x_tmp, ymin, x_tmp, ymax, col='gray', lwd=1.5)
}
dev.off()


######### plot FSiP factor plot
FSiP = cbind(sig1_fsip, sig2_fsip, sig3_fsip, sig4_fsip)
pdf(paste(outputname, '.FSiP.box.pdf', sep=''))
boxplot(FSiP, ylim=c(0,1), range = 1.5)
abline(h=ref_fsip, col='blue', lwd=1.5)
dev.off()








