import os
import numpy as np
from scipy import stats

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used, sep):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(sep)]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
def quantile_normalization(data_matrix, signal_col, outputname):
	###### read inputs
	### read data matrix
	data_matrix = read2d_array(data_matrix, str, '\t')

	### split header; info matrix; signal matrix
	header = [data_matrix[0,:]]
	data_matrix_info = data_matrix[1:, 0:signal_col-1]
	data_matrix_sig = np.array(data_matrix[1:, signal_col-1:], dtype=float)

	print('data_matrix_sig.shape:')
	print(data_matrix_sig.shape)
	######
	### put all sig with the same rank into one vector in a dict
	rank_sig_dict = {}
	data_rank_matrix = []
	for i in range(0, data_matrix_sig.shape[1]):
		### get sample data
		sample_data = data_matrix_sig[:,i]
		### get sample data rank
		sample_data_rank=stats.rankdata(sample_data, method='ordinal')
		### generate data rank matrix
		data_rank_matrix.append(sample_data_rank)

		### put sig into rank_sig_dict
		for sig, rank in zip(sample_data, sample_data_rank):
			if not (rank in rank_sig_dict):
				rank_sig_dict[rank] = [ sig ]
			else:
				rank_sig_dict[rank].append(sig)

	### transpose data_rank_matrix to original shape
	data_rank_matrix = np.transpose( np.array(data_rank_matrix) )
	print('data_rank_matrix.shape:')
	print(np.array(data_rank_matrix).shape)
	
	###### quantile normalization
	quantile_norm_sig={}
	for i in rank_sig_dict:
		### calculate mean of signal with the same rank
		qn=np.mean( rank_sig_dict[i] )
		### put rank mean sig to quantile_norm_sig dict
		quantile_norm_sig[i]=qn

	###### get quantile normalized matrix
	data_qn_matrix = []
	for rs in data_rank_matrix:
		qn_tmp = []
		for r in rs:
			### replace data_rank_matrix by rank_mean_sig
			qn_tmp.append(quantile_norm_sig[r])
		data_qn_matrix.append(qn_tmp)

	### add back data info matrix
	data_qn_matrix = np.concatenate((data_matrix_info, data_qn_matrix),  axis = 1)
	### add back header
	data_qn_matrix = np.concatenate((header, data_qn_matrix),  axis = 0)
	
	###### write quantile normed data
	write2d_array(data_qn_matrix, outputname)


############################################################################
############################################################################
#time python quantile_normalization.py -i all_rna_sampe_target_class.txt -n 6 -o all_rna_sampe_target_class_qn.txt

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:n:o:")
	except getopt.GetoptError:
		print 'time python quantile_normalization.py -i data_matrix -n signal_col -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python quantile_normalization.py -i data_matrix -n signal_col -o outputname'
			sys.exit()
		elif opt=="-i":
			data_matrix=str(arg.strip())
		elif opt=="-n":
			signal_col=int(arg.strip())	
		elif opt=="-o":
			outputname=str(arg.strip())		

	quantile_normalization(data_matrix, signal_col, outputname)

if __name__=="__main__":
	main(sys.argv[1:])