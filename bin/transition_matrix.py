#module load python/2.7
import os
import numpy as np
from subprocess import call
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm
import seaborn

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used, split='\t'):
	import numpy as np
	data=open(filename,'r')
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



state_input = 'test_20cell.location_pknorm_16lim_refmean_samplemean_nosub.state'

state_all = read2d_array(state_input, str, ' ') 
state = state_all[:,4:state_all.shape[1]-1].astype(int)

state_freq = np.bincount(state.flatten())
state_freq = np.unique(state, return_counts=True)


kmer_count = {}
for i in range(0, state.shape[1]):
	wg = state[:,i]
	for j in range(0, wg.shape[0]-1):
		kmer_ij = str(wg[j]) + str(wg[j+1])
		if kmer_ij in kmer_count:
			kmer_count[kmer_ij] = kmer_count[kmer_ij] + 1
		else:
			kmer_count[kmer_ij] = 1


state_transition_matrix = []
allnum = sum(state_freq[1])
for i in range(0,state_freq[0].shape[0]):
	state_transition_i = []
	for j in range(0,state_freq[0].shape[0]):
		if str(i)+str(j) in kmer_count:
			counts = kmer_count[str(i)+str(j)]
		else:
			counts = 0
		exp_i = state_freq[1][i]
		exp_j = state_freq[1][j]
		exp = exp_i * exp_j / allnum
		if i != j:
			fc = (counts+100) / (exp+100)
		else:
			fc = 1
		state_transition_i.append(fc)
	state_transition_matrix.append(state_transition_i) 

state_transition_matrix = np.array(state_transition_matrix, dtype=float)

write2d_array(state_transition_matrix, 'state_transition_matrix.txt')

fig, ax = plt.subplots(figsize=(10,10))  
seaborn.clustermap(state_transition_matrix, annot=True, linewidths=.5, ax=ax)
fig.savefig('state_transition_matrix.png')






