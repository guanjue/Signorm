library(metap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

cell_marker = args[1]
tail = args[2]
input_folder = args[3]
siglim = as.numeric(args[4])
### extract filenames of the cell marker
file_list = list.files(input_folder, pattern=paste('^', cell_marker, '(.*)', tail, '$', sep='') )
print(file_list)
### read files of the cell marker
data_matrix = NULL
for (file in file_list){
	d = read.table(paste(input_folder, file, sep=''), header = F)
	data_matrix = cbind(data_matrix, d[,])
}
### get fisher method combined p-value
get_fisher_p = function(x){
	if (length(x)!=1){
		x_p = 10^(-x)
		fp = sumlog(x_p)$p
		if (fp<=0.1^siglim){
			fp = 0.1^siglim
		}		
	} else{
		fp = 10^(-x)
		if (fp<=0.1^siglim){
			fp = 0.1^siglim
		}		
	}

	fp_neglog10 = -log10(fp)
	return(fp_neglog10)
}

fisher_p = apply(data_matrix, 1, function(x) get_fisher_p(x))

### write output
write.table(fisher_p, paste(cell_marker, '.fisher_p.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

