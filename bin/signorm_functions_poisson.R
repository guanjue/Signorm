##############################################
t_r_curve_change_point = function(t_r_matrix, changepoint_method, max_cp_num, t_r_change_point_plot_file_name, ignore_t_lim_lower, ignore_t_lim_upper, raw_plot_lim, mean_or_var, fit_polynorm, polynomial_degree){
	library(LSD)
	library(changepoint)

	### ignore the data points with t value 1-ignore_t_lim (too noisy)
	t_r_matrix_lower = t_r_matrix[t_r_matrix[,1]>ignore_t_lim_lower,]
	t_r_matrix = t_r_matrix_lower[t_r_matrix_lower<ignore_t_lim_upper,]
	print(dim(t_r_matrix))
	### extract t value and r value
	t_od = (t_r_matrix[,1])
	r_od = (t_r_matrix[,2])/(t_r_matrix[,3]) ### log2 transform r

	### remove x1 OR x2 equals 0
	t = t_od[as.logical((!is.na(r_od)) * (is.finite(r_od)) ) ]
	r = r_od[as.logical((!is.na(r_od)) * (is.finite(r_od)) ) ]

	### polynomial regression fit the log(r) vs log(t) pattern
	print('fit polynomial regression model')
	data_for_polyfit = data.frame(x=t, y=r)

	lo = lm(y~poly(log(x), polynomial_degree, raw=TRUE), data = data_for_polyfit)
	### get polynomial regression fitted model predicted value
	lo_fit_value = predict(lo, newdata=data.frame(x=t))

	### variance change-point without polynomial regression norm
	print('find variance change-point without polynomial regression norm')
	if (mean_or_var=='var'){
			ansvar=(cpt.var(r, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=max_cp_num))
	} else if (mean_or_var=='mean'){
			ansvar=(cpt.mean(r, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=max_cp_num))
	} else if (mean_or_var=='meanvar'){
			ansvar=(cpt.meanvar(r, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=max_cp_num))
	}
	
	if (fit_polynorm=='polynorm'){
		### variance change-point with polynomial regression norm (cosider 0 for read the same sample data)
		print('find variance change-point with polynomial regression norm')
		if (max(r)!=0){
			r_straight = r-lo_fit_value
			if (mean_or_var=='var'){
				ansvar_norm=(cpt.var(r_straight, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=max_cp_num))
			} else if (mean_or_var=='mean'){
				ansvar_norm=(cpt.mean(r_straight, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=max_cp_num))
			} else if (mean_or_var=='meanvar'){
				ansvar_norm=(cpt.meanvar(r_straight, class=FALSE, method = changepoint_method, penalty = 'BIC', Q=max_cp_num))
			}
			print(ansvar_norm)
		} else{
				ansvar_norm=c(0)
		}
	} else{
		ansvar_norm=ansvar
	}
	### use the first change point as t_threshold
	if (ansvar_norm[1] != 0){ 
		t_variance_change_point = t[ansvar_norm[1]]
	} else{
		t_variance_change_point = 0
	}
	
	### plot r vs t pattern without polynomial regression norm & and variance change point
	png(paste(t_r_change_point_plot_file_name, '.raw.png', sep=''))
	par(mfrow=c(1,1))
	heatscatter(t, r, pch = 20, ylim=c(-raw_plot_lim,raw_plot_lim), xlim=c(1,10000), log='x', main=paste(toString(t[ansvar[1]]), 'VS', toString(t[ansvar_norm[1]]), sep=' ') )
	lines(t, lo_fit_value, col = 'blue', lty=2)
	abline(v = t[ansvar[1]], col = 'gray', lty=2)
	abline(v = t_variance_change_point, col = 'red', lty=2)
	dev.off()

	### plot r vs t pattern with polynomial regression norm & and variance change point
	png(paste(t_r_change_point_plot_file_name, '.polynorm.png', sep=''))
	par(mfrow=c(1,1))
	heatscatter(t, r-lo_fit_value, pch = 20, ylim=c(-raw_plot_lim,raw_plot_lim), xlim=c(1,10000), log='x', main=toString(t[ansvar_norm[1]]))
	abline(v = t_variance_change_point, col = 'red', lty=2)
	dev.off()

	output = list('tcp'=t_variance_change_point, 'pnmodel'=lo)
	return(output)
} 
##############################################

##############################################
scale_factor = function(data_x, data_y, method){
	### get mean
	if (method=='mean'){
			merge_x = mean(data_x)
			merge_y = mean(data_y)
		} else{
			merge_x = median(data_x)
			merge_y = median(data_y)
		}

	### get scale factor log scale
	sf = merge_y / merge_x
	print(sf)

	### get M & A for MA plot
	M = log(data_y/1) - log(data_x/1)
	A = 0.5*(log(data_y/1) + log(data_x/1) )

	### remove 0s for plotting
	for_plotting_id = as.logical( (!is.na(M)) * (!is.na(A)) )
	M = M[for_plotting_id]
	A = A[for_plotting_id]

	output = list('sf'=sf, 'M'=M, 'A'=A, 'merge_x'=merge_x, 'merge_y'=merge_y)
	return(output)
}
##############################################

##############################################
plotting_scatterplot_MAplot = function(data_x_high_t, data_y_high_t, data_x_low_t, data_y_low_t, data_x, data_y, log, od_A, od_M, sf_info_total_mean, sf_info_low_t, sf_info_high_t, output_file_name){
	### plot sig_vs_sig only
	print('scatter plot')
	png(paste(output_file_name, '.high_t.hs.png', sep =''))
	if (log == 'T'){
		heatscatter(data_x_high_t, data_y_high_t, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
	} else{
		heatscatter(data_x_high_t, data_y_high_t, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
	}
	abline(0,1,col = 'blue')
	dev.off()

	### plot low_vs_sig only
	print('scatter plot')
	png(paste(output_file_name, '.low_t.hs.png', sep =''))
	if (log == 'T'){
		heatscatter(data_x_low_t, data_y_low_t, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
	} else{
		heatscatter(data_x_low_t, data_y_low_t, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
	}
	abline(0,1,col = 'blue')
	dev.off()

	print('scatter plot')
	png(paste(output_file_name, '.all.log.hs.png', sep =''))
	if (log == 'T'){
		heatscatter(data_x, data_y, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
	} else{
		heatscatter(data_x, data_y, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
	}
	abline(0,1,col = 'blue')
	dev.off()

	###### MA plot for different strategy
	#png(paste(output_file_name, '.4.MA.png', sep =''))
	#par(mfrow=c(2,2))
	###
	#heatscatter(od_A, od_M, pch = 20, ylim=c(-10,10), main='original signal', xlim=c(-1,8.5))
	#abline(h=0,col = 'red')
	###
	#heatscatter(sf_info_total_mean$A, sf_info_total_mean$M, pch = 20, ylim=c(-10,10), main='original signal mean', xlim=c(-1,8.5))
	#abline(h=0,col = 'red')
	###
	#heatscatter(sf_info_low_t$A, sf_info_low_t$M, pch = 20, ylim=c(-10,10), main='low t mean', xlim=c(-1,8.5))
	#abline(h=0,col = 'red')
	###
	#heatscatter(sf_info_high_t$A, sf_info_high_t$M, pch = 20, ylim=c(-10,10), main='high t mean', xlim=c(-1,8.5))
	#abline(h=0,col = 'red')
	###
	#dev.off()	
}
##############################################

##############################################
calculate_scale_factor_with_t_thresh = function(data_x, data_y, sampling_num, seed, t_threshold, ignore_t_lim, quantile_lim, scatterplot_MAplot_output_file_name){
	### if sampling_num != 0, sampling calculate scale factor & plotting 
	if (sampling_num != 0){
		set.seed(seed)
		used_id = sample(length(data_y)[1],sampling_num)
		data_x = data_x[used_id,]
		data_y = data_y[used_id,]		
	}

	### get input t value
	data_t = data_x+data_y

	### get t threshold
	print('t threshold')
	print(sum(data_t>t_threshold))
	print(t_threshold)

	### if data_t<=t_threshold have more than X% of the bins, use X% bin as the threshold
	data_t_ignore_ts = data_t[data_t>ignore_t_lim] 
	if ( (sum(data_t_ignore_ts<=t_threshold) / length(data_t_ignore_ts)) > quantile_lim ){
		print('use X% quantile')
		t_threshold = quantile(data_t_ignore_ts, quantile_lim)
		print(t_threshold)
		write.table(t_threshold, paste(scatterplot_MAplot_output_file_name,'.threshold_qunatile.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
	}

	### get bins with t <= t_threshold
	data_x_low_t = data_x[data_t<=t_threshold]
	data_y_low_t = data_y[data_t<=t_threshold]

	### get bins with t <= t_threshold
	data_x_high_t = data_x[data_t>t_threshold]
	data_y_high_t = data_y[data_t>t_threshold]

	### remove zero for plotting
	data_x_non0 = data_x[as.logical((data_x!=0) * (data_y!=0))]
	data_y_non0 = data_y[as.logical((data_x!=0) * (data_y!=0))]

	### remove zero for plotting
	data_x_low_t_non0 = data_x_low_t[as.logical((data_x_low_t!=0) * (data_y_low_t!=0))]
	data_y_low_t_non0 = data_y_low_t[as.logical((data_x_low_t!=0) * (data_y_low_t!=0))]
	print('length(data_y_low_t_non0):')
	print(length(data_y_low_t_non0))

	### remove zero for plotting
	data_x_high_t_non0 = data_x_high_t[as.logical((data_x_high_t!=0) * (data_y_high_t!=0))]
	data_y_high_t_non0 = data_y_high_t[as.logical((data_x_high_t!=0) * (data_y_high_t!=0))]
	print('length(data_y_high_t_non0):')
	print(length(data_y_high_t_non0))

	### original signal M & A
	od_M = log(data_y) - log(data_x) 
	od_A = 0.5*(log(data_y) + log(data_x) )
	### remove 0s for plotting
	for_plotting_id = as.logical( (!is.na(od_M)) * (!is.na(od_A)) )
	od_M = od_M[for_plotting_id]
	od_A = od_A[for_plotting_id]

	### sf total mean
	sf_info_total_mean = scale_factor(data_x, data_y, 'mean')

	### sf total median
	sf_info_total_median = scale_factor(data_x, data_y, 'median')

	### scale low_t
	sf_info_low_t = scale_factor(data_x_low_t, data_y_low_t, 'mean')

	### scale high_t
	sf_info_high_t = scale_factor(data_x_high_t, data_y_high_t, 'mean')

	if (t_threshold > 0){
		### plot scatterplot & MAplot
		log = 'T' ### use log scale
		plotting_scatterplot_MAplot(data_x_high_t, data_y_high_t, data_x_low_t, data_y_low_t, data_x, data_y, log, od_A, od_M, sf_info_total_mean, sf_info_low_t, sf_info_high_t, scatterplot_MAplot_output_file_name)
		### get scale factor vector
		sf_vector = c(sf_info_total_mean$merge_x, sf_info_total_median$merge_x, sf_info_low_t$merge_x, sf_info_high_t$merge_x,   sf_info_total_mean$merge_y, sf_info_total_median$merge_y, sf_info_low_t$merge_y, sf_info_high_t$merge_y)
	} else{
		sf_vector = c(sf_info_total_mean$merge_x, sf_info_total_median$merge_x, sf_info_total_mean$merge_x, sf_info_total_mean$merge_x,   sf_info_total_mean$merge_y, sf_info_total_median$merge_y, sf_info_total_mean$merge_y, sf_info_total_mean$merge_y)
	}

	sf_vector_t_threshold = list("sf_vector" = sf_vector, "t_threshold" = t_threshold)
	
	return(sf_vector_t_threshold)
}




