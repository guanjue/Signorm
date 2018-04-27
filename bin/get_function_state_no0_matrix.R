d_fuc = read.table('atac_20cell.function.matrix.txt', header=F)
d_sig = read.table('atac_20cell.signal.matrix.txt', header=F)
d_index = read.table('atac_20cell.index.matrix.txt', header=F)

### get all function
all_fun = d_fuc[,-c(1,2,3,4)]
### get all state sum
all_state_sum = rowSums(all_fun)

### remove all 0 state
d_fuc_no0 = d_fuc[all_state_sum!=0,]
d_sig_no0 = d_sig[all_state_sum!=0,]
d_index_no0 = d_index[all_state_sum!=0,]

write.table(d_fuc_no0, 'atac_20cell.function.matrix.no0.txt', quote=F, row.names=F, col.names=F, sep='\t')
write.table(d_sig_no0, 'atac_20cell.signal.matrix.no0.txt', quote=F, row.names=F, col.names=F, sep='\t')
write.table(d_index_no0, 'atac_20cell.index.matrix.no0.txt', quote=F, row.names=F, col.names=F, sep='\t')


