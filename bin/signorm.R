### get parameters
args = commandArgs(trailingOnly=TRUE)

xais_variable_file = args[1]
yais_variable_file = args[2]
input_file_t_r_matrix = args[3]

t_r_change_point_plot_file_name = paste(input_file_t_r_matrix, '.cp.png', sep='')
scatterplot_MAplot_output_file_name = paste(input_file_t_r_matrix, '.scatter_MA.png', sep='')

data_x_sig_norm_output_file = args[4]
signal_scale_factor_vector_output_file = paste(data_x_sig_norm_output_file, '.sf_vec.txt', sep='')

mean_or_var = args[5] ### mean or var
changepoint_method = args[6]
max_cp_num = as.numeric(args[7])
fit_polynorm = args[8]
polynomial_degree = as.numeric(args[9])

sampling_num = as.numeric(args[10])
quantile_lim = as.numeric(args[11])
seed = as.numeric(args[12])
ignore_t_lim_lower = as.numeric(args[13])
ignore_t_lim_upper = as.numeric(args[14])
raw_plot_lim = as.numeric(args[15])

scale_factor_type = as.numeric(args[16]) ### 1: total mean; 2:total median; 3: low Poisson mean; 4: high Poisson mean

source_code_folder = args[17]

### signorm functions
source(paste(source_code_folder, 'signorm_functions.R', sep = ''))

### read r vs t table of r vs t: r=sum(x1)/sum(x2); t=x1+x2
t_r_matrix = read.table(input_file_t_r_matrix,header = F)
print(dim(t_r_matrix))
### get t threshold
t_threshold = t_r_curve_change_point(t_r_matrix, changepoint_method, max_cp_num, t_r_change_point_plot_file_name, ignore_t_lim_lower, ignore_t_lim_upper, raw_plot_lim, mean_or_var, fit_polynorm, polynomial_degree)
print((t_threshold))
### read input reads table
data_x_od = read.table(xais_variable_file, header = FALSE)
data_x_sig = as.matrix(data_x_od[,1]) 

data_y_od = read.table(yais_variable_file, header = FALSE)
data_y_sig = as.matrix(data_y_od[,1]) 

### get scale factor based on signal part
signal_scale_factor_vector_t_threshold_modified = calculate_scale_factor_with_t_thresh(data_x_sig, data_y_sig, sampling_num, seed, t_threshold, ignore_t_lim_lower, quantile_lim, scatterplot_MAplot_output_file_name)
signal_scale_factor_vector = signal_scale_factor_vector_t_threshold_modified$sf_vector
t_threshold_modified = signal_scale_factor_vector_t_threshold_modified$t_threshold

print('signal_scale_factor_vector')
print(signal_scale_factor_vector)
print('t_threshold_modified')
print(t_threshold_modified)
is.element(scale_factor_type, c(1,2,3,4))
### norm the x-axis signal by the scale factor
if (is.element(scale_factor_type, c(1,2,3,4))){
	print('1: total mean; 2:total median; 3: low Poisson mean; 4: high Poisson mean')
	print(scale_factor_type)
	data_x_sig_norm = data_x_sig / signal_scale_factor_vector[scale_factor_type] * signal_scale_factor_vector[scale_factor_type+4]
} else if (scale_factor_type=5) {
	### 2 scale factor for 2 parts
	print('2 scale factor for 2 parts method')
	sf_low = 1/signal_scale_factor_vector[3] * signal_scale_factor_vector[3+4]
	sf_high = 1/signal_scale_factor_vector[4] * signal_scale_factor_vector[4+4]
	data_x_sig_norm = apply(cbind(data_x_sig, data_y_sig), 1, function(x) if((x[1]+x[2])<=t_threshold_modified) { x * sf_low} else { x * sf_high } )
}

### write normed signal vector
write.table(data_x_sig_norm, data_x_sig_norm_output_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

### write scale facotr vector
write.table(signal_scale_factor_vector, signal_scale_factor_vector_output_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
