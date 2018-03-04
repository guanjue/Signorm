#module load python/2.7
import os
import numpy as np
from subprocess import call

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

################################################################################################
### FSiP
def FSiP(wg_bed, peak_bed, sig2_col_list, sig2_wg_raw):
	sig2_col_list = sig2_col_list.split(',')

	sig2_output_name = sig2_wg_raw.split('.')[0] + '.' + sig2_wg_raw.split('.')[1]

	######
	### get sig2 columns
	sig2_col_id_plus = ''
	for i in range(0,len(sig2_col_list)-1):
		sig2_col_id_plus = sig2_col_id_plus + '$' + str(sig2_col_list[i]) + '+'
	sig2_col_id_plus = sig2_col_id_plus + '$' + str(sig2_col_list[len(sig2_col_list)-1])

	### get sig2 column
	print('get sig2 column...')
	call('tail -n+2 ' + peak_bed + ' | awk -F \'\t\' -v OFS=\'\t\' \'{' + 'if (' + sig2_col_id_plus + ' > 0)' + 'print $1, $2, $3, $4, ' + sig2_col_id_plus + ' }\' > ' + sig2_output_name + '.pk.bed', shell=True)
	### intersect with wg bed
	call('bedtools map' + ' -a ' + sig2_output_name + '.pk.bed' + ' -b ' + wg_bed + ' -c 5 -o sum -null 0' + ' > ' + sig2_output_name + '.pk.sig.bed', shell=True)

	### read table
	wg_table = read2d_array(wg_bed, str)
	wg_sig = wg_table[:,4].astype(float)

	pk_table = read2d_array(sig2_output_name + '.pk.sig.bed', str)
	pk_sig = pk_table[:,5].astype(float)

	pk_region_len = np.sum(pk_table[:,2].astype(float) - pk_table[:,1].astype(float))
	pk_num = pk_sig.shape[0]
	pk_sig_sum = np.sum(pk_sig)
	wg_sig_sum = np.sum(wg_sig)
	FSiP = pk_sig_sum / wg_sig_sum

	print('pk length: ' + str(pk_region_len))
	print('pk num: ' + str(pk_num))
	print('pk signal: ' + str(pk_sig_sum))
	print('wg signal: ' + str(wg_sig_sum))
	print('FSiP: ' + str(FSiP))

	write2d_array(np.array([[FSiP, pk_sig_sum, wg_sig_sum, pk_region_len, pk_num]]), sig2_output_name+'.info.txt')
	write2d_array(np.array(np.concatenate((wg_table[:,0:3], np.array((wg_sig+0.1)/(wg_sig_sum)*8315353-0.1, dtype=str).reshape((wg_sig.shape[0], 1))), axis=1)), sig2_output_name+'.totalmean.bedgraph')

	call('rm ' + sig2_output_name + '.pk.bed', shell=True)
	call('rm ' + sig2_output_name + '.pk.sig.bed', shell=True)
############################################################################
#time python /Volumes/MAC_Data/data/labs/zhang_lab/01projects/signorm/bin/peaknorm_order.py -w 200_noblack.11_22_2017.bed -p atacTable180124peaksFiltered.txt -n 100000 -a $sig1_col -b $sig1'.upperlim.txt' -c $sig2_col -d $sig2'.upperlim.txt' -u 32

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hw:p:a:b:")
	except getopt.GetoptError:
		print 'time python index_label2meansig.py -i input_file_list -o outputname -l log2 -s small_num -d lowerlim -u upperlim'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python index_label2meansig.py -i input_file_list -o outputname -l log2 -s small_num -d lowerlim -u upperlim'
			sys.exit()
		elif opt=="-w":
			wg_bed=str(arg.strip())
		elif opt=="-p":
			peak_bed=str(arg.strip())
		elif opt=="-a":
			sig2_col_list=str(arg.strip())
		elif opt=="-b":
			sig2_wg_raw=str(arg.strip())

	FSiP(wg_bed, peak_bed, sig2_col_list, sig2_wg_raw)

if __name__=="__main__":
	main(sys.argv[1:])



