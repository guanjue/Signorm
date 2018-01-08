import os
import numpy as np
from subprocess import call
from scipy import stats
from scipy.stats import nbinom

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

### get bakcground adjusted -log10(p-value) with negative binomial distribution
def nbp_bg_adj(sample_sig_file, input_sig_file, outputname):
	### read input files
	sample = read2d_array(sample_sig_file, float)
	background = read2d_array(input_sig_file, float)

	### set threshold to ignore
	thesh = 0

	### get sample mean & var 
	sample_non0 = sample[sample>thesh]
	sample_mean = np.mean(sample_non0)
	sample_var = np.var(sample_non0)

	### get negative binomial parameters from sample track regions
	sample_prob = sample_mean / sample_var
	if sample_prob<0.1:
		sample_prob = 0.1
	if sample_prob>=0.9:
		sample_prob = 0.9

	### get size parameter for negative binomial distribution p-value (1st round)
	sample_size = sample_mean * sample_prob / (1-sample_prob)

	### get background mean & var
	background_non0 = background[background>thesh]
	bg_mean = np.mean(background_non0)
	bg_var = np.var(background_non0)

	print('check input track overdispersion in background regions, var/mean=' + str(round(bg_var/bg_mean, 3)) )
	print(sample_prob)
	print(sample_size)
	print(len(background_non0))

	print(bg_mean)
	print(bg_var)

	### 1st round negative binomial p-value
	i=0
	nb_pval_list = np.empty((0,), float)
	print(nb_pval_list.shape)
	for sig in sample:
		if i%10000 == 0:
			print(i)
		i=i+1
		nb_pval_tmp = np.array( 1 - nbinom.cdf(sig, sample_size, sample_prob, loc=0) )
		print(nb_pval_tmp.shape)
		nb_pval_list = np.concatenate((nb_pval_list, nb_pval_tmp))

	############### second round
	### get sample bg regions
	sample_bg = sample[nb_pval_list>=0.001,]
	sample_bg_non0 = sample_bg[sample_bg>thesh]
	sample_bg_mean = np.mean(sample_bg_non0)
	sample_bg_var = np.var(sample_bg_non0)

	print('check signal track overdispersion in background regions, var/mean=' + str(round(sample_bg_var/sample_bg_mean, 3)) )
	print(sample_bg_mean)
	print(sample_bg_var)
	print(len(sample_bg_non0))

	### get negative binomial parameters from signal track bg regions
	sample_bg_prob = sample_bg_mean / sample_bg_var
	if sample_bg_prob<0.1:
		sample_bg_prob = 0.1

	if sample_bg_prob>=0.9:
		sample_bg_prob = 0.9


	### get size parameter for negative binomial distribution p-value (2nd round)
	sample_bg_size = sample_bg_mean * sample_bg_prob / (1-sample_bg_prob)

	### get background bg regions
	background_bg = background[nb_pval_list>=0.001,]
	background_bg_non0 = background_bg[background_bg>thesh]
	background_bg_mean = np.mean(background_bg_non0)
	background_bg_var = np.var(background_bg_non0)

	print('check input track overdispersion in background regions, var/mean=' + str(round(background_bg_var/background_bg_mean, 3)) )
	print(sig_bg_prob)
	print(sig_bg_size)
	print(len(input_bg_non0))

	print(input_bg_mean)
	print(inpy_bg_var)

	### 2nd round negative binomial p-value
	i=0
	nb_pval_list = np.empty((0,), float)
	for sig, bg in zip(sample, background):
		if i%10000 == 0:
			print(i)
		i = i+1
		nb_pval_tmp = 1 - nbinom.cdf(sig, sample_bg_size * (bg+1)/(background_bg_mean+1), sample_bg_prob, loc=0)
		nb_pval_list = np.concatenate((nb_pval_list,np.array([nb_pval_tmp])))
	### convert to np array
	nb_pval_list = -np.log10(nb_pval_list)
	nb_pval_list = nb_pval_list.reshape(nb_pval_list.shape[0], 1)

	### write output
	write2d_array(nb_pval_list, outputname+'.nbp_2r.txt')

	### info vector
	mvsp = open(outputname+'.mvsp.txt', 'w')
	mvsp.write(str(sample_bg_mean) + '\t')
	mvsp.write(str(sample_bg_var) + '\t')
	mvsp.write(str(sample_bg_size) + '\t')
	mvsp.write(str(sample_bg_prob) + '\n')
	mvsp.close()

###########################################
# time python negative_binomial_p_2r_bgadj.py -i <sample_track_signal> -t <background_track_signal> -o tmp.txt
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hs:b:o:")
	except getopt.GetoptError:
		print 'python negative_binomial_p_2r_bgadj.py -i <sample_track_signal> -t <background_track_signal> -o <output filename>'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python negative_binomial_p_2r_bgadj.py -i <sample_track_signal> -t <background_track_signal> -o <output filename>'
			sys.exit()
		elif opt=="-s":
			sample_sig_file=str(arg.strip())
		elif opt=="-b":
			input_sig_file=str(arg.strip())
		elif opt=="-o":
			outputname=str(arg.strip())
	nbp_bg_adj(sample_sig_file, input_sig_file, outputname)

if __name__=="__main__":
	main(sys.argv[1:])



