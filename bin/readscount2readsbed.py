#module load python/2.7
import os
import numpy as np

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
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

def readscount2readsbed(input_reads_count_file, output_reads_file):
	rc = read2d_array(input_reads_count_file, str)

	reads = []
	i = 0
	for rc_i in rc:
		if i%100000 ==0:
			print(i)
		i = i+1
		num = int(rc_i[3])
		for j in range(0, j):
			reads.append(rc_i[0:3])

	reads = np.array(reads)
	write2d_array(reads, output_reads_file)



############################################################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:")
	except getopt.GetoptError:
		print 'time python readscount2readsbed.py -i input_reads_count_file -o output_reads_file'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python readscount2readsbed.py -i input_reads_count_file -o output_reads_file'
			sys.exit()
		elif opt=="-i":
			input_reads_count_file=str(arg.strip())
		elif opt=="-o":
			output_reads_file=str(arg.strip())

	readscount2readsbed(input_reads_count_file, output_reads_file)

if __name__=="__main__":
	main(sys.argv[1:])



