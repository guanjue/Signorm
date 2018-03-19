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



def ideas_state_color_by_atac_sig(ideas_state_bed, ideas_state_col, atac_sig_bed, atac_sig_col, ideas_state_id_color_name_list, signal_upperlim, signal_lowerlim, outputname):
	### read ideas state
	ideas_state_all = read2d_array(ideas_state_bed, str)

	### read atac signal
	atac_sig_all = read2d_array(atac_sig_bed, str)

	### read ideas_state_id_color_name_list
	ideas_state_info = read2d_array(ideas_state_id_color_name_list, str)
	ideas_state_info_dict = {}
	for infos in ideas_state_info:
		ideas_state_info_dict[infos[0]] = [ infos[1], infos[2] ]

	ideas_state_bigbed = []
	for ideas_peak, atac_peak in zip(ideas_state_all, atac_sig_all):
		ideas_pk_state_chr = ideas_peak[0]
		ideas_pk_state_start = ideas_peak[1]
		ideas_pk_state_end = ideas_peak[2]
		ideas_pk_state_id = ideas_peak[ideas_state_col-1]
		ideas_pk_state_color = ideas_state_info_dict[ideas_pk_state_id][0]
		ideas_pk_state_name = ideas_state_info_dict[ideas_pk_state_id][1]
		atac_pk_sig = atac_peak[atac_sig_col-1]

		### change rgb color based on atac-seq signal
		signal_range = signal_upperlim - signal_lowerlim
		sigdif = (signal_upperlim - float(atac_pk_sig)) / signal_range
		ideas_pk_state_atac_sig_color_r = int(float(ideas_pk_state_color[0]) * (1-sigdif) + 255.0 * sigdif)
		ideas_pk_state_atac_sig_color_g = int(float(ideas_pk_state_color[1]) * (1-sigdif) + 255.0 * sigdif)
		ideas_pk_state_atac_sig_color_b = int(float(ideas_pk_state_color[2]) * (1-sigdif) + 255.0 * sigdif)
		ideas_pk_state_atac_sig_color = str(ideas_pk_state_atac_sig_color_r)+','+str(ideas_pk_state_atac_sig_color_g)+','+str(ideas_pk_state_atac_sig_color_b)

		### merge all information
		ideas_state_bigbed.append([ ideas_pk_state_chr, ideas_pk_state_start, ideas_pk_state_end, ideas_pk_state_name+';'+ideas_pk_state_id, '1000', '.', ideas_pk_state_start, ideas_pk_state_end, ideas_pk_state_atac_sig_color ])


	ideas_state_bigbed = np.array(ideas_state_bigbed)

	write2d_array(ideas_state_bigbed, outputname)



###########################################
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"ha:b:c:d:e:u:l:o:")
	except getopt.GetoptError:
		print 'python ideas_state_color_by_atac_sig.py -a ideas_state_bed -b ideas_state_col -c atac_sig_bed -d atac_sig_col -e ideas_state_id_color_name_list -u signal_upperlim -l signal_lowerlim o outputname'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python ideas_state_color_by_atac_sig.py -a ideas_state_bed -b ideas_state_col -c atac_sig_bed -d atac_sig_col -e ideas_state_id_color_name_list -u signal_upperlim -l signal_lowerlim o outputname'
			sys.exit()
		elif opt=="-a":
			ideas_state_bed=str(arg.strip())
		elif opt=="-b":
			ideas_state_col=int(arg.strip())
		elif opt=="-c":
			atac_sig_bed=str(arg.strip())
		elif opt=="-d":
			atac_sig_col=int(arg.strip())
		elif opt=="-e":
			ideas_state_id_color_name_list=str(arg.strip())
		elif opt=="-u":
			signal_upperlim=float(arg.strip())
		elif opt=="-l":
			signal_lowerlim=float(arg.strip())

		elif opt=="-o":
			outputname=str(arg.strip())


	ideas_state_color_by_atac_sig(ideas_state_bed, ideas_state_col, atac_sig_bed, atac_sig_col, ideas_state_id_color_name_list, signal_upperlim, signal_lowerlim, outputname)

if __name__=="__main__":
	main(sys.argv[1:])






