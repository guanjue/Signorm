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


cd_tree_file = 'cd_tree.txt'
cd_tree = read2d_array(cd_tree_file, str, ',') 


state_input = 'test_20cell.location_pknorm_16lim_refmean_samplemean_nosub.state'

state_input = 'test.state'

state_all = read2d_array(state_input, str, ' ')
state_all_ct = state_all[0]

kmer_count = {}
for ct in cd_tree:
	s = ct[0]
	t = ct[1]
	s_col_id = np.where(state_all_ct==s)[0][0]
	t_col_id = np.where(state_all_ct==t)[0][0]
	s_state_array = state_all[1:,s_col_id]
	t_state_array = state_all[1:,t_col_id]
	for s_state,t_state in zip(s_state_array, t_state_array):
		kmer_ij = str(s_state) + str(t_state)
		if kmer_ij in kmer_count:
			kmer_count[kmer_ij] = kmer_count[kmer_ij] + 1
		else:
			kmer_count[kmer_ij] = 1


state = state_all[1:,4:state_all.shape[1]-1].astype(int)

state_freq = np.bincount(state.flatten())
state_freq = np.unique(state, return_counts=True)


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

write2d_array(state_transition_matrix, 'state_transition_matrix_ct.txt')

fig, ax = plt.subplots(figsize=(10,10))  
seaborn.heatmap(state_transition_matrix, annot=True, linewidths=.5)
fig.savefig('state_transition_matrix_ct.png')


fig, ax = plt.subplots(figsize=(10,10))  
fig = seaborn.clustermap(state_transition_matrix)
fig.savefig('state_transition_matrix_ct.png')




