import numpy as np

def total_reads_norm(sig1, sig2, method, sig1_normed_output):
	data1 = open(sig1,'r')
	data2 = open(sig2,'r')

	### skip header
	#data1.readline()
	#data2.readline()

	sig1 = []
	sig2 = []
	for rec1,rec2 in zip(data1, data2): ### loop matrix at x axis
		d1 = [ x.strip() for x in rec1.split('\t') ]
		d2 = [ x.strip() for x in rec2.split('\t') ]

		### add to vector
		sig1.append(d1[0])
		sig2.append(d2[0])		
	data1.close()
	data2.close()
	### vector to float np.array
	sig1 = np.array(sig1, dtype=float)
	sig2 = np.array(sig2, dtype=float)

	### get total reads scale factor
	if method == 'mean':
		scale_factor = np.mean(sig2) / np.mean(sig1)
	elif method != 'mean':
		scale_factor = np.percentile(sig2, float(method)) / np.percentile(sig1, float(method))

	### get normed signal1
	sig1_normed = np.around(sig1 * scale_factor, decimals = 0)

	### write output
	result = open(sig1_normed_output,'w')
	for records in sig1_normed:
		result.write(str(records)+'\n')
	result.close()


###########################################
# python total_reads_norm.py -a signal_1_file -b signal_2_file -m mean -o normed_sig1_output_file
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ha:b:m:o:")
	except getopt.GetoptError:
		print 'python total_reads_norm.py -a signal_1_file -b signal_2_file -m method(mean or percentile #) -o normed_sig1_output_file'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python total_reads_norm.py -a signal_1_file -b signal_2_file -m method(mean or percentile #) -o normed_sig1_output_file'
			sys.exit()
		elif opt=="-a":
			sig1=str(arg.strip())
		elif opt=="-b":
			sig2=str(arg.strip())
		elif opt=="-m":
			method=str(arg.strip())
		elif opt=="-o":
			sig1_normed_output=str(arg.strip())
	total_reads_norm(sig1, sig2, method, sig1_normed_output)

if __name__=="__main__":
	main(sys.argv[1:])

