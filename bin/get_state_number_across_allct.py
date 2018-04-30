#module load python/2.7
import os
import numpy as np

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data.readline()
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(' ')]
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
### get state number
def get_state_number_across_allct(ideas_state_matrix_file, target_state_label):
	target_state = target_state_label
	ideas_state_matrix = read2d_array(ideas_state_matrix_file, str)

	ideas_state_wig = []
	i=0
	for records in ideas_state_matrix:
		if i%100000==0:
			print(i)
		i = i+1
		state = records[4:len(records)-1]
		state_num = np.sum(state==target_state)
		ideas_state_wig.append([records[1],records[2],records[3],str(state_num)])

	ideas_state_wig = np.array(ideas_state_wig)
	write2d_array(ideas_state_wig, target_state+'.state_num.bedgraph')



############################################################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:s:")
	except getopt.GetoptError:
		print 'time python get_state_number_across_allct.py -i ideas_state_matrix -s target_state_label'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python get_state_number_across_allct.py -i ideas_state_matrix -s target_state_label'
			sys.exit()
		elif opt=="-i":
			ideas_state_matrix_file=str(arg.strip())
		elif opt=="-s":
			target_state_label=str(arg.strip())

	get_state_number_across_allct(ideas_state_matrix_file, target_state_label)

if __name__=="__main__":
	main(sys.argv[1:])

