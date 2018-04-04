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


list1 = read2d_array('/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/peak_list.txt', str)

ccRE_type = read2d_array('/storage/home/gzx103/scratch/vision/5end/index_set_20cell_pknorm/ideas_atac_sig_color/ideas_state_id_color_name_list.txt', str)[:,2]

ccRE_type_matrix = []

for name in list1:
	filename = 'atac_pk.' + name + '.bed'
	d_dict = {}
	d = read2d_array(filename, str)
	for info in d:
		if info[3] in d_dict:
			d_dict[info[3]] = d_dict[info[3]] +1
		else:
			d_dict[info[3]] = 0
	### get to 1d dict
	count_array = []
	for ccRE in ccRE_type:
		if ccRE in d_dict:
			count_array.append(d_dict[ccRE])
		else:
			count_array.append(0)
		ccRE_type_matrix.append(count_array)

ccRE_type_matrix = np.array(ccRE_type_matrix)

output = open('count_table.txt', 'w')
output.write('celltyp'+'\t')
i = 0
for records in ccRE_type_matrix:
	output.write(list1[i]+'\t')
	for count in records:
		output.write(str(count)+'\t')
	output.write('\n')
	i = i+1

output.close()







