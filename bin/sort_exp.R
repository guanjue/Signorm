d=read.table('rnaHtseqCountsall.cmp_ery.bedmatched.txt', header= F)


rpkm_cmp = d[,7]*1000000/51105650/(d[,3]-d[,2])+10
rpkm_ery = d[,8]*1000000/51105650/(d[,3]-d[,2])+10

d_s = d[order(-(rpkm_cmp/rpkm_ery),]
d_s = d_s[d_s[,1]!='chrM',c(1:6)]
write.table(d_s, 'gencode.vM4.annotation.pc.cmp_ery_fc.bed', quote=F, col.names=F, row.names=F, sep='\t')


