library(pheatmap)

file_list = read.table('ideas_list.txt', header = F)
names = read.table('ideas_list.txt', header = F, sep='.')[,2]

count_vector = c()

for (i in seq(1,dim(file_list)[1])){
	print(i)
	ideas_label_tmp = read.table(toString(file_list[i,1]), header=F)[,5]
	count_vector_tmp = table(ideas_label_tmp)
	count_vector = cbind(count_vector, count_vector_tmp)
}

row_sum = apply(count_vector,1, sum)

count_vector_p = count_vector / row_sum
colnames(count_vector_p) = file_list[,2]
colnames(count_vector) = file_list[,2]
count_vector_p_add100 = (count_vector+100) / (row_sum+100*20)

breaksList = seq(0, 1, by = 0.01)
my_colorbar=colorRampPalette(c('white', 'blue'))(n = length(breaksList))

colnames(count_vector_p) = names

pdf('ideas_state_count_p.pdf')
pheatmap(count_vector_p, color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

pdf('ideas_state_count_p_cluster.pdf')
pheatmap(count_vector_p, color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

pdf('ideas_state_count_p_add100.pdf')
pheatmap(count_vector_p_add100, color=my_colorbar, breaks = breaksList, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()

write.table(count_vector_p, 'count_vector_p.txt', sep='\t', quote=FALSE)
write.table(count_vector_p_add100, 'count_vector_p_add100.txt', sep='\t', quote=FALSE)


