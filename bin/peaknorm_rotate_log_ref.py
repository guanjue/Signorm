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
### p-value adjust (fdr & bonferroni)
def p_adjust(pvalue, method):
	p = pvalue
	n = len(p)
	p0 = np.copy(p, order='K')
	nna = np.isnan(p)
	p = p[~nna]
	lp = len(p)
	if method == "bonferroni":
		p0[~nna] = np.fmin(1, lp * p)
	elif method == "fdr":
		i = np.arange(lp, 0, -1)
		o = (np.argsort(p))[::-1]
		ro = np.argsort(o)
		p0[~nna] = np.fmin(1, np.minimum.accumulate((p[o]/i*lp)))[ro]
	else:
		print "Method is not implemented"
		p0 = None
	return p0

################################################################################################
### gradient descent
def gradientDescent(sig1_pk,sig1_bg, sig2_pk,sig2_bg, A, B, alpha, beta, numIterations):
	best_loss0 = 1e+10
	p = 0
	for i in range(0, numIterations):
		h_sig2_pk0 = A*(sig2_pk**(1*B))
		h_sig2_bg0 = A*(sig2_bg**(1*B))
		h_sig2_pk0_mean = np.mean(h_sig2_pk0)
		h_sig2_bg0_mean = np.mean(h_sig2_bg0)
		sig1_pk_mean = np.mean(sig1_pk**1)
		sig1_bg_mean = np.mean(sig1_bg**1)

		#loss0 = abs(h_sig2_pk0_mean - sig1_pk_mean) + abs(h_sig2_bg0_mean - sig1_bg_mean)
		loss0 =abs( (h_sig2_pk0_mean / h_sig2_bg0_mean) - (sig1_pk_mean / sig1_bg_mean) )
		#loss0 = abs(np.sqrt(np.mean(h_sig2_pk0**2)) - np.sqrt(np.mean(sig1_pk)**2+np.var(sig1_pk))) + abs(np.sqrt(np.mean(h_sig2_bg0**2)) - np.sqrt(np.mean(sig1_bg)**2+np.var(sig1_bg)))

		### next step
		h_sig2_pk_B = A*(sig2_pk**(1*(B+beta)))
		h_sig2_bg_B = A*(sig2_bg**(1*(B+beta)))
		h_sig2_pk0_mean_B = np.mean(h_sig2_pk_B)
		h_sig2_bg0_mean_B = np.mean(h_sig2_bg_B)
		loss_B = abs( (h_sig2_pk0_mean_B / h_sig2_bg0_mean_B) - (sig1_pk_mean / sig1_bg_mean) )
		#loss_B = abs(np.sqrt(np.mean(h_sig2_pk_B**2)) - np.sqrt(np.mean(sig1_pk)**2+np.var(sig1_pk))) + abs(np.sqrt(np.mean(h_sig2_bg_B**2)) - np.sqrt(np.mean(sig1_bg)**2+np.var(sig1_bg)))

		print(loss0)
		if loss0 < best_loss0:
			best_AB = [A, B]
			best_loss0 = loss0
			p=0
		else:
			p = p+1
			if p > 20:
				print('opt')
				break
		# avg cost per example (the 2 in 2*m doesn't really matter here.
		# But to be consistent with the gradient, I include it)
		print("Iteration %d | Cost: %f" % (i, loss0))
		if loss0 <= 0.0000000001:
			print('converged!')
			break
		# avg gradient per example
		gradientB = (-loss0+loss_B)
		print(gradientB)
		# update
		B = B - 1e-3 * gradientB / abs(gradientB) * abs(loss0)
		A = sig1_bg_mean / h_sig2_bg0_mean_B

		print([A,B])
	print('best')
	print(best_AB)
	print(best_loss0)
	return np.array(best_AB)

################################################################################################
###
def pknorm(wg_bed, peak_bed, sample_num, sig1_col_list, sig1_wg_raw, sig2_col_list, sig2_wg_raw, upperlim, lowerlim):
	sig1_col_list = sig1_col_list.split(',')
	sig2_col_list = sig2_col_list.split(',')

	sig1_output_name = sig1_wg_raw.split('.')[0]+'_'+sig1_wg_raw.split('.')[1]
	sig2_output_name = sig2_wg_raw.split('.')[0]+'_'+sig2_wg_raw.split('.')[1]

	### read whole genome signals
	sig1 = read2d_array(sig1_wg_raw, float)
	sig2 = read2d_array(sig2_wg_raw, float)

	### add small_number
	small_num = 1e-1

	### total reads norm
	print('ref sum')
	print(np.sum(sig1))
	'''
	sig1 = sig1 / np.sum(sig1) * 5000000
	print(np.sum(sig1))
	#sig1[sig1 > upperlim] = upperlim
	'''
	if sig1_output_name == sig2_output_name:
		sig2 = sig1
	


	### read whole genome binary label
	#sig1_binary = p_adjust(10**(-sig1), 'fdr') < 0.05
	#bg1_binary = p_adjust(10**(-sig1), 'fdr') >= 0.05
	sig1_binary = 10**(-sig1) <= 0.001
	bg1_binary = 10**(-sig1) > 0.001
	print(sum(sig1_binary))
	#sig2_binary = p_adjust(10**(-sig2), 'fdr') < 0.05
	#bg2_binary = p_adjust(10**(-sig2), 'fdr') >= 0.05
	sig2_binary = 10**(-sig2) <= 0.001
	bg2_binary = 10**(-sig2) > 0.001
	print(sum(sig2_binary))

	### peak region (both != 0 in sig1 & sig2)
	peak_binary_pk = (sig1_binary[:,0] * sig2_binary[:,0]) != 0
	print(sum(peak_binary_pk))
	peak_binary = peak_binary_pk & (sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]) #& (sig1[:,0] < upperlim) & (sig2[:,0] < upperlim)
	print(sum(peak_binary))
	### background region (both == 0 in sig1 & sig2)
	bg_binary_bg = (sig1_binary[:,0] + sig2_binary[:,0]) == 0
	print(sum(bg_binary_bg))
	bg_binary = bg_binary_bg & (sig1[:,0] != sig1[0,0]) & (sig2[:,0] != sig2[0,0]) #& (sig1[:,0] < upperlim) & (sig2[:,0] < upperlim)
	print(sum(bg_binary))

	### get transformation factor
	AB = gradientDescent(sig1[sig1_binary[:,0],0]+small_num,sig1[bg_binary,0]+small_num, sig2[sig2_binary[:,0],0]+small_num,sig2[bg_binary,0]+small_num, 1.0, 1.0, 0.0001, 0.0001, 200)
	A=AB[0]
	B=AB[1]

	print('transformation: '+'B: '+str(B)+'; A: '+str(A))
	### transformation
	sig2_norm = []
	for s in sig2[:,0]:
		s = s
		s_norm = (A*(s+small_num)**B)-small_num
		if s_norm >= upperlim:
			s_norm = upperlim
		elif s_norm <= lowerlim:
			s_norm = lowerlim
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

	### rotated means for sig2 for plotting
	sig1_1log_pk_m_od = np.log2(np.mean(sig1[sig1_binary[:,0],0]+small_num))
	sig2_1log_pk_m_od = np.log2(np.mean(sig2[sig2_binary[:,0],0]+small_num))

	sig1_1log_bg_m_od = np.log2(np.mean(sig1[bg1_binary[:,0],0]+small_num))
	sig2_1log_bg_m_od = np.log2(np.mean(sig2[bg2_binary[:,0],0]+small_num))

	sig2_1log_pk_m_pkn = np.log2(np.mean(sig2_norm[sig2_binary[:,0],0]+small_num))
	sig2_1log_bg_m_pkn = np.log2(np.mean(sig2_norm[bg2_binary[:,0],0]+small_num))

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
	peak_binary_sample = peak_binary_pk[idx]
	bg_binary_sample = bg_binary_bg[idx]
	plot_x = np.log2(sig2_norm[idx,0]+small_num)
	plot_y = np.log2(sig1[idx,0]+small_num)
	lims_max = np.max(np.concatenate((plot_x, plot_y)))
	lims_min = np.min(np.concatenate((plot_x, plot_y)))

	plt.figure()
	plt.scatter(plot_x, plot_y, marker='.', color='dodgerblue')
	plt.scatter(plot_x[bg_binary_sample], plot_y[bg_binary_sample], marker='.', color='gray')
	plt.scatter(plot_x[peak_binary_sample], plot_y[peak_binary_sample], marker='.', color='coral')
	plt.scatter(sig2_1log_pk_m_pkn, sig1_1log_pk_m_od, marker='.', color='k')
	plt.scatter(sig2_1log_bg_m_pkn, sig1_1log_bg_m_od, marker='.', color='k')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'k')
	plt.plot([sig2_1log_bg_m_pkn, sig2_1log_pk_m_pkn], [sig1_1log_bg_m_od, sig1_1log_pk_m_od])
	#plt.scatter(np.mean(plot_x[peak_binary_sample]), np.mean(plot_y[peak_binary_sample]), marker='.', color='k')
	#plt.scatter(np.mean(plot_x[bg_binary_sample]), np.mean(plot_y[bg_binary_sample]), marker='.', color='k')
	#plt.plot([np.mean(plot_x[bg_binary_sample]), np.mean(plot_x[peak_binary_sample])], [np.mean(plot_y[bg_binary_sample]), np.mean(plot_y[peak_binary_sample])])
	plt.xlabel(sig2_output_name + '.pknorm')
	plt.ylabel(sig1_output_name + '.pknorm')
	plt.xlim(lims_min, lims_max)
	plt.ylim(lims_min, lims_max)
	plt.savefig(sig2_output_name + '.pknorm.scatterplot.png')


	plot_x = np.log2(sig2[idx,0]+small_num)
	plot_y = np.log2(sig1[idx,0]+small_num)
	lims_max = np.max(np.concatenate((plot_x, plot_y)))
	lims_min = np.min(np.concatenate((plot_x, plot_y)))

	plt.figure()
	plt.scatter(plot_x, plot_y, marker='.', color='dodgerblue')
	plt.scatter(plot_x[bg_binary_sample], plot_y[bg_binary_sample], marker='.', color='gray')
	plt.scatter(plot_x[peak_binary_sample], plot_y[peak_binary_sample], marker='.', color='coral')
	plt.scatter(sig2_1log_pk_m_od, sig1_1log_pk_m_od, marker='.', color='k')
	plt.scatter(sig2_1log_bg_m_od, sig1_1log_bg_m_od, marker='.', color='k')
	plt.plot([lims_min, lims_max], [lims_min, lims_max], 'k', color = 'k')
	plt.plot([sig2_1log_bg_m_od, sig2_1log_pk_m_od], [sig1_1log_bg_m_od, sig1_1log_pk_m_od])
	#plt.scatter(np.mean(plot_x[peak_binary_sample]), np.mean(plot_y[peak_binary_sample]), marker='.', color='k')
	#plt.scatter(np.mean(plot_x[bg_binary_sample]), np.mean(plot_y[bg_binary_sample]), marker='.', color='k')
	#plt.plot([np.mean(plot_x[bg_binary_sample]), np.mean(plot_x[peak_binary_sample])], [np.mean(plot_y[bg_binary_sample]), np.mean(plot_y[peak_binary_sample])])
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



