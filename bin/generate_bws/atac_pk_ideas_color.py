#module load python/2.7
import os
import numpy as np
from subprocess import call

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used, split='\t'):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(split)]
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


def atac_pk_ideas_color(atac_binary_matrix_file, atac_start_col, atac_pk_list, ideas_state_matrix_file, ideas_start_col, ideas_state_list, ideas_state_id_color_name_list, outputname):
	### read atac pk
	atac_binary_matrix = read2d_array(atac_binary_matrix_file, str)
	atac_pk_bed = atac_binary_matrix[:, 0:3]
	atac_pk_binary = atac_binary_matrix[:, (atac_start_col-1):]

	atac_pk_list = read2d_array(atac_pk_list, str)
	### matrix to dict
	atac_pk_dict = {}
	for i in range(0, atac_pk_binary.shape[1]):
		ct = atac_pk_list[i][1]
		atac_pk_dict[ct] = atac_pk_binary[:,i]


	### read atac pk ideas state
	ideas_state_matrix = read2d_array(ideas_state_matrix_file, str)
	ideas_state_matrix_id = ideas_state_matrix[:, (ideas_start_col-1):]

	ideas_state_list = read2d_array(ideas_state_list, str)
	### matrix to dict
	ideas_state_dict = {}
	for i in range(0, ideas_state_matrix_id.shape[1]):
		ct = ideas_state_list[i][1]
		ideas_state_dict[ct] = ideas_state_matrix_id[:,i]


	### read ideas_state_id_color_name_list
	ideas_state_info = read2d_array(ideas_state_id_color_name_list, str)
	ideas_state_info_dict = {}
	for infos in ideas_state_info:
		ideas_state_info_dict[infos[0]] = [ infos[1], infos[2] ]


	### binary pk colored by ideas state
	for ct in atac_pk_dict:
		print(ct)
		if ct in ideas_state_dict:
			ct_atac_pk_ideas_color_bed = []
			atac_pk_vec = atac_pk_dict[ct].reshape(atac_pk_dict[ct].shape[0],1)
			ideas_state_vec = ideas_state_dict[ct].reshape(ideas_state_dict[ct].shape[0],1)
			### concatenate bed, atac_pk, ideas_state
			bed_atac_pk_ideas_state = np.concatenate((atac_pk_bed, atac_pk_vec, ideas_state_vec), axis = 1)

			for info in bed_atac_pk_ideas_state:
				### check if atac_pk != 0
				if info[3] != '0':
					name = ideas_state_info_dict[info[4]][1] + ';' + info[4]
					color = ideas_state_info_dict[info[4]][0]
					ct_atac_pk_ideas_color_bed.append([ info[0],info[1],info[2], name, '1000', '.', info[1],info[2], color ])

			ct_atac_pk_ideas_color_bed = np.array(ct_atac_pk_ideas_color_bed)
			write2d_array(ct_atac_pk_ideas_color_bed, outputname+'.'+ct+'.bed')


###########################################
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ha:b:c:d:e:f:g:o:")
	except getopt.GetoptError:
		print 'python atac_pk_ideas_color.py -a atac_binary_matrix_file -b atac_start_col -c atac_pk_list -d ideas_state_matrix_file -e ideas_start_col -f ideas_state_list -g ideas_state_id_color_name_list -o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python atac_pk_ideas_color.py -a atac_binary_matrix_file -b atac_start_col -c atac_pk_list -d ideas_state_matrix_file -e ideas_start_col -f ideas_state_list -g ideas_state_id_color_name_list -o outputname'
			sys.exit()
		elif opt=="-a":
			atac_binary_matrix_file=str(arg.strip())
		elif opt=="-b":
			atac_start_col=int(arg.strip())
		elif opt=="-c":
			atac_pk_list=str(arg.strip())
		elif opt=="-d":
			ideas_state_matrix_file=str(arg.strip())
		elif opt=="-e":
			ideas_start_col=int(arg.strip())
		elif opt=="-f":
			ideas_state_list=str(arg.strip())
		elif opt=="-g":
			ideas_state_id_color_name_list=str(arg.strip())

		elif opt=="-o":
			outputname=str(arg.strip())


	atac_pk_ideas_color(atac_binary_matrix_file, atac_start_col, atac_pk_list, ideas_state_matrix_file, ideas_start_col, ideas_state_list, ideas_state_id_color_name_list, outputname)

if __name__=="__main__":
	main(sys.argv[1:])








