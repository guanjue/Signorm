d=read.table('rnaHtseqCountsall.cmp_ery.bedmatched.txt', header= F)

d_s = d[order(-d[,8]),]
d_s = d_s[d_s[,1]!='chrM',c(1:6)]
write.table(d_s, 'gencode.vM4.annotation.pc.cmp_ery_fc.bed', quote=F, col.names=F, row.names=F, sep='\t')


