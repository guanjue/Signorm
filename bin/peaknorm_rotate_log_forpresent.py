#module load python/2.7
import os
import numpy as np
from subprocess import call
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
###
def pknorm(wg_bed, peak_bed, sample_num, sig1_col_list, sig1_wg_raw, sig2_col_list, sig2_wg_raw, upperlim, lowerlim):
	sig1_col_list = sig1_col_list.split(',')
	sig2_col_list = sig2_col_list.split(',')

	sig1_output_name = sig1_wg_raw.split('.')[0]
	sig2_output_name = sig2_wg_raw.split('.')[0]
	######
	### get sig1 columns
	sig1_col_id_plus = ''
	for i in range(0,len(sig1_col_list)-1):
		sig1_col_id_plus = sig1_col_id_plus + '$' + str(sig1_col_list[i]) + '+'
	sig1_col_id_plus = sig1_col_id_plus + '$' + str(sig1_col_list[len(sig1_col_list)-1])
	### get sig1 column
	call('tail -n+2 ' + peak_bed + ' | awk -F \'\t\' -v OFS=\'\t\' \'{' + 'if (' + sig1_col_id_plus + ' > 0)' + 'print $1, $2, $3, $4, ' + sig1_col_id_plus + ' }\' > ' + sig1_output_name + '.bed', shell=True)
	### intersect with wg bed
	call('bedtools intersect' + ' -a ' + wg_bed + ' -b ' + sig1_output_name + '.bed' + ' -c' + ' > ' + sig1_output_name + '.wg.bed', shell=True)

	call('cut -f4 ' + sig1_output_name + '.wg.bed' + ' > ' + sig1_output_name + '.wg.txt', shell=True)

	######
	### get sig2 columns
	sig2_col_id_plus = ''
	for i in range(0,len(sig2_col_list)-1):
		sig2_col_id_plus = sig2_col_id_plus + '$' + str(sig2_col_list[i]) + '+'
	sig2_col_id_plus = sig2_col_id_plus + '$' + str(sig2_col_list[len(sig2_col_list)-1])

	### get sig2 column
	call('tail -n+2 ' + peak_bed + ' | awk -F \'\t\' -v OFS=\'\t\' \'{' + 'if (' + sig2_col_id_plus + ' > 0)' + 'print $1, $2, $3, $4, ' + sig2_col_id_plus + ' }\' > ' + sig2_output_name + '.bed', shell=True)
	### intersect with wg bed
	call('bedtools intersect' + ' -a ' + wg_bed + ' -b ' + sig2_output_name + '.bed' + ' -c' + ' > ' + sig2_output_name + '.wg.bed', shell=True)
	### extract binary column
	call('cut -f4 ' + sig2_output_name + '.wg.bed' + ' > ' + sig2_output_name + '.wg.txt', shell=True)



	### read whole genome signals
	sig1 = read2d_array(sig1_wg_raw, float)
	sig2 = read2d_array(sig2_wg_raw, float)

	### read whole genome binary label
	sig1_binary = read2d_array(sig1_output_name + '.wg.txt', int)
	sig2_binary = read2d_array(sig2_output_name + '.wg.txt', int)

	### peak region (both != 0 in sig1 & sig2)
	peak_binary = (sig1_binary[:,0] * sig2_binary[:,0]) != 0
	peak_binary = peak_binary & (sig1_binary[:,0] < upperlim) & (sig2_binary[:,0] < upperlim)
	### background region (both == 0 in sig1 & sig2)
	bg_binary = (sig1_binary[:,0] + sig2_binary[:,0]) == 0
	bg_binary = bg_binary & (sig1_binary[:,0] < upperlim) & (sig2_binary[:,0] < upperlim)

	### get transformation factor
	B = ( np.log2(np.mean(sig1[peak_binary,0]+1)) - np.log2(np.mean(sig1[bg_binary,0]+1)) ) / ( np.log2(np.mean(sig2[peak_binary,0]+1)) - np.log2(np.mean(sig2[bg_binary,0]+1)) )
	A = np.mean(sig1[peak_binary,0]+1) / (np.mean(sig2[peak_binary,0]+1) ** B)

	print('transformation: '+'B: '+str(B)+'; A: '+str(A))
	### transformation
	sig2_norm = []
	for s in sig2[:,0]:
		if (s > lowerlim) and (s < upperlim):
			s_norm = A * (s ** B)
			if s_norm >= upperlim:
				s_norm = upperlim
			elif s_norm <= lowerlim:
				s_norm = lowerlim
		elif (s >= upperlim) or (s <= lowerlim):
			s_norm = s
		sig2_norm.append(s_norm)

	sig2_norm = np.array(sig2_norm, float)
	
	print(sig2[0:10])
	print(sig2_norm[0:10])
	### total reads sf (for compare)
	total_mean_sf = np.sum(sig1) / np.sum(sig2)

	### convert to float np.array
	sig2_norm = np.array(sig2_norm, float)
	### reshape for writing oputput
	sig2_norm = np.reshape(sig2_norm, (sig2_norm.shape[0],1))

	###FRiP score
	sig2_norm_FRiP = np.sum(sig2_norm[(sig2_binary[:,0]!=0),0]) / np.sum(sig2_norm)
	sig2_FRiP = np.sum(sig2[(sig2_binary[:,0]!=0),0]) / np.sum(sig2)
	sig1_FRiP = np.sum(sig1[(sig1_binary[:,0]!=0),0]) / np.sum(sig1)

	### write output: normalized signal
	write2d_array(sig2_norm, sig2_output_name + '.pknorm.txt')

	### write output: sf & FRiP
	info = np.array([[total_mean_sf, B, A], [sig1_FRiP, sig2_norm_FRiP, sig2_FRiP]])
	write2d_array(info, sig2_output_name + '.info.txt')


	### plot scatter plot
	np.random.seed(2018)
	idx = np.random.randint(sig2_norm.shape[0], size=sample_num)
	peak_binary_sample = peak_binary[idx]
	bg_binary_sample = bg_binary[idx]
	plot_x = np.log2(sig2_norm[idx,0]+1)
	plot_y = np.log2(sig1[idx,0]+1)
	lims_max = np.max(np.concatenate((plot_x, plot_y)))
	lims_min = np.min(np.concatenate((plot_x, plot_y)))

	plt.figure()
	plt.scatter(plot_x, plot_y, marker='.')
	plt.scatter(plot_x[peak_binary_sample], plot_y[peak_binary_sample], marker='.', color='coral')
	plt.scatter(np.mean(plot_x[peak_binary_sample]), np.mean(plot_y[peak_binary_sample]), marker='.', color='r')
	plt.scatter(np.mean(plot_x[bg_binary_sample]), np.mean(plot_y[bg_binary_sample]), marker='.', color='r')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'r')
	plt.xlabel(sig2_output_name + '.pknorm')
	plt.ylabel(sig1_output_name + '.pknorm')
	plt.xlim(lims_min, lims_max)
	plt.ylim(lims_min, lims_max)
	plt.savefig(sig2_output_name + '.pknorm.scatterplot.png')


	plot_x = np.log2(sig2[idx,0]+1)
	plot_y = np.log2(sig1[idx,0]+1)
	lims_max = np.max(np.concatenate((plot_x, plot_y)))
	lims_min = np.min(np.concatenate((plot_x, plot_y)))

	plt.figure()
	plt.scatter(plot_x, plot_y, marker='.')
	plt.scatter(plot_x[peak_binary_sample], plot_y[peak_binary_sample], marker='.', color='coral')
	plt.scatter(np.mean(plot_x[peak_binary_sample]), np.mean(plot_y[peak_binary_sample]), marker='.', color='r')
	plt.scatter(np.mean(plot_x[bg_binary_sample]), np.mean(plot_y[bg_binary_sample]), marker='.', color='r')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'r')
	plt.xlabel(sig2_output_name + '.pknorm')
	plt.ylabel(sig1_output_name + '.pknorm')
	plt.xlim(lims_min, lims_max)
	plt.ylim(lims_min, lims_max)
	plt.savefig(sig2_output_name + '.scatterplot.png')

############################################################################
#time python plot_violin.py -i atac_list.txt -o atac_list -l T -s 2 -l 4 -u 100

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hw:p:n:a:b:c:d:u:l:")
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
		elif opt=="-n":
			sample_num=int(arg.strip())
		elif opt=="-a":
			sig1_col_list=str(arg.strip())			
		elif opt=="-b":
			sig1_wg_raw=str(arg.strip())		
		elif opt=="-c":
			sig2_col_list=str(arg.strip())			
		elif opt=="-d":
			sig2_wg_raw=str(arg.strip())		
		elif opt=="-u":
			upperlim=float(arg.strip())
		elif opt=="-l":
			lowerlim=float(arg.strip())

	pknorm(wg_bed, peak_bed, sample_num, sig1_col_list, sig1_wg_raw, sig2_col_list, sig2_wg_raw, upperlim, lowerlim)

if __name__=="__main__":
	main(sys.argv[1:])



