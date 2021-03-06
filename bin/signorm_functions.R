##############################################
t_r_curve_change_point = function(t_r_matrix, changepoint_method, max_cp_num, t_r_change_point_plot_file_name, ignore_t_lim_lower, ignore_t_lim_upper, raw_plot_lim, mean_or_var, fit_polynorm, polynomial_degree, round_factor, round_type){
	library(LSD)
	library(changepoint)

	### ignore the data points with t value 1-ignore_t_lim (too noisy)
	t_r_matrix_lower = t_r_matrix[t_r_matrix[,1]>ignore_t_lim_lower,]
	t_r_matrix = t_r_matrix_lower[t_r_matrix_lower<ignore_t_lim_upper,]
	print(dim(t_r_matrix))
	### extract t value and r value
	### round t for robustness
	if (round_type=='log2'){
		### add 2 because r add 1 (see line 26)
		t_od_round = 2**(round( log2(t_r_matrix[,1]) / round_factor) * round_factor)
	} else{
		t_od_round = (round( (t_r_matrix[,1]) / round_factor) * round_factor) 
	}
	
	### use convert data frame
	t_r_matrix_df = as.data.frame(cbind(t_od_round, t_r_matrix[,2], t_r_matrix[,3]))
	### add read sun with the same rounded t
	t_r_matrix_round = aggregate(. ~ t_od_round, data=t_r_matrix_df, FUN=mean)

	### get robust r
	t_od = t_r_matrix_round[,1]
	r_od = log2(t_r_matrix_round[,2]+1)-log2(t_r_matrix_round[,3]+1) ### log2 transform r

	### remove x1 OR x2 equals 0
	t = t_od[as.logical((!is.na(r_od)) * (is.finite(r_od)) ) ]
	r = r_od[as.logical((!is.na(r_od)) * (is.finite(r_od)) ) ]
	print(summary(t))
	print(summary(r))
	### polynomial regression fit the log(r) vs log(t) pattern
	print('fit polynomial regression model')
	data_for_polyfit = data.frame(x=t, y=r)

	#lo = lm(y~poly(log(x), polynomial_degree, raw=TRUE), data = data_for_polyfit)
	lo = loess(y ~ log(x), span=0.5, data = data_for_polyfit, control=loess.control(surface="direct"))
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
	if (length(data_x_high_t)>2){
		png(paste(output_file_name, '.high_t.hs.png', sep =''))
		if (log == 'T'){
			heatscatter(data_x_high_t, data_y_high_t, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
		} else{
			heatscatter(data_x_high_t, data_y_high_t, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
		}
		abline(0,1,col = 'blue')
		dev.off()	
	}


	### plot low_vs_sig only
	print('scatter plot')
	if (length(data_x_low_t)>2){
		png(paste(output_file_name, '.low_t.hs.png', sep =''))
		if (log == 'T'){
			heatscatter(data_x_low_t, data_y_low_t, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
		} else{
			heatscatter(data_x_low_t, data_y_low_t, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
		}
		abline(0,1,col = 'blue')
		dev.off()
	}

	print('scatter plot')
	if (length(data_x)>2){
		png(paste(output_file_name, '.all.log.hs.png', sep =''))
		if (log == 'T'){
			heatscatter(data_x, data_y, log='xy', pch = 20, ylim=c(1,10000), xlim=c(1,10000))
		} else{
			heatscatter(data_x, data_y, pch = 20, ylim=c(-10,5), xlim=c(20,10000))
		}
		abline(0,1,col = 'blue')
		dev.off()
	}
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
calculate_scale_factor_with_t_thresh = function(data_x, data_y, sampling_num, seed, t_threshold, ignore_t_lim, quantile_lim, scatterplot_MAplot_output_file_name, round_factor, round_type){
	### if sampling_num != 0, sampling calculate scale factor & plotting 
	t_all = data_x+data_y
	if (sampling_num != 0){
		set.seed(seed)
		used_id = sample(length(data_y)[1],sampling_num)
		data_x = data_x[used_id,]
		data_y = data_y[used_id,]		
	}

	### get input t value
	data_t = data_x+data_y

	### round t for robustness
	if (round_type=='log2'){
			data_t = 2**(round(log2(data_t) / round_factor) * round_factor )
		} else{
			data_t = (round((data_t) / round_factor) * round_factor )
		}
	

	### get t threshold
	print('t threshold')
	print(sum(data_t>t_threshold))
	print(t_threshold)

	### if data_t<=t_threshold have more than X% of the bins, use X% bin as the threshold
	data_t_passlim = data_t[data_t>(ignore_t_lim)] ### add 2 because t-r fit t add 2 & rs add 1
	data_used_p = sum(data_t_passlim<=t_threshold)/length(data_t_passlim)
	if ( data_used_p >= 0 ){
		print('use X% quantile')
		print(summary(data_t_passlim))
		t_threshold = quantile(data_t_passlim, quantile_lim, type=1)
		print(t_threshold)
		write.table(t_threshold, paste(scatterplot_MAplot_output_file_name,'.threshold_qunatile.txt', sep=''), quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
	}

	### get bins with t <= t_threshold
	data_x_low_t = data_x[data_t<=t_threshold]
	data_y_low_t = data_y[data_t<=t_threshold]

	### get bins with t <= t_threshold
	data_x_high_t = data_x[data_t>t_threshold]
	data_y_high_t = data_y[data_t>t_threshold]

	### get background vs foreground bins
	bg_fg_10 = ((t_all)<=t_threshold)*1 ### t_all add 2 because t-r fit t add 2 & rs add 1 so threshod added 2
	print(head(bg_fg_10))
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

	sf_vector_t_threshold = list("sf_vector" = sf_vector, "t_threshold" = t_threshold, 'bg_fg_10' = bg_fg_10)
	
	return(sf_vector_t_threshold)
}
##############################################

##############################################
MAnorm = function(data_x_sig, data_y_sig, sampling_num, seed, MAplot_output_file_name){
	### get M & A
	small_num = 0.5
	M = log2(data_x_sig + small_num) - log2(data_y_sig + small_num)
	A = 0.5*(log2(data_x_sig + small_num) + log2(data_y_sig + small_num))
	### sort M & A for plotting
	O = order(A)
	a = A[O]
	m = M[O]
	### sampling for loess fit
	### if sampling_num != 0, sampling calculate scale factor & plotting 
	if (sampling_num != 0){
		### keep the order
		used_id = round(seq(1, length(a), len = sampling_num/20))
		a = a[used_id]
		m = m[used_id]		
	}
	### fit loess line
	fit = loess(m ~ a, span=0.25)
	### get bias line for MA plot
	bias = predict(fit, newdata = data.frame(a = a))
	### get scale factors
	tsf = 1/(2**(predict(fit, newdata = data.frame(a = A))))
	### plot orignal MA plot
	png(paste(MAplot_output_file_name, '.scatterplot.png', sep=''))
	heatscatter(log2(data_x_sig + small_num), log2(data_y_sig + small_num), pch = 20, main='original signal')
	abline(0,1,col = 'blue')
	dev.off()	
	png(paste(MAplot_output_file_name, '.od_MA.png', sep=''))
	heatscatter(a, m, pch = 20, main='original signal')
	abline(h=0,col = 'blue')
	lines(a, bias, col='red', lty=1)
	dev.off()
	### plot bias corrected MA plot
	png(paste(MAplot_output_file_name, '.bias_corrected_MA.png', sep=''))
	heatscatter(a, m-bias, pch = 20, main='bias corrected signal')
	abline(h=0,col = 'blue')	
	dev.off()
	###
	return(tsf)
}


mode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

signorm_robust = function(d1, d2, p, start_point, step, cor_lim, plot_name, sampling_num, use_log_axis, ignore_sig, siglim){
	r=r2=NULL
	ignore_sig_1 = 0#mode(d1)
	ignore_sig_2 = 0#mode(d2)
	used_range = rev( p^seq(start_point, 5, step) )
	#used_range = rev(seq(start_point, 1, step))
	d12 = d1+d2
	for (i in seq(1,length(used_range))){
		used_ida = as.logical( (d1>quantile(d1[as.logical((d1>ignore_sig_1) * (d1<siglim))], 1-used_range[i])) * (d2>quantile(d2[as.logical((d2>ignore_sig_2) * (d2<siglim))], 1-used_range[i])) * (d1<siglim) * (d2<siglim) )
		#used_ida = as.logical( (d12>quantile(d12[ as.logical((d1<siglim) * (d2<siglim)) ], 1-used_range[i])) * (d1<siglim) * (d2<siglim) )
		#r2[i] = cor( (d1[used_ida]), (d2[used_ida]), method = 'spearman')
		d12_log = (cbind(log(d1[used_ida]), log(d2[used_ida])))
		r2[i] = cor( d12_log[,1], d12_log[,2], method = 'pearson')
		r[i] = sum((d1[used_ida])) / sum((d2[used_ida]))
		print(paste(i, r2[i], r[i], sep='_'))
	}

	set.seed(2017)
	used_id = sample(dim(d1)[1], 100000)

	d1_s = d1[used_id]
	d2_s = d2[used_id]

	#ansvar_norm=(cpt.meanvar(r2, class=FALSE, method = 'BinSeg', penalty = 'BIC', Q=3))
	#used_r2 = ansvar_norm[1] #
	used_r2 = which.max(r2)

	d1_thresh = quantile(d1[as.logical((d1>ignore_sig_1) * (d1<siglim))], 1-used_range[which.max(r2)] )
	d2_thresh = quantile(d2[as.logical((d2>ignore_sig_2) * (d2<siglim))], 1-used_range[which.max(r2)] )
	print(d1_thresh)
	print(d2_thresh)

	if (r2[used_r2]>=cor_lim){
		used_idb = as.logical( (d1>quantile(d1[as.logical((d1>ignore_sig_1) * (d1<siglim))], 1-used_range[which.max(r2)])) * (d2>quantile(d2[as.logical((d2>ignore_sig_2) * (d2<siglim))], 1-used_range[which.max(r2)])) * (d1<siglim) * (d2<siglim) )
		#used_idb = as.logical( (d12>quantile(d12[ as.logical((d1<siglim) * (d2<siglim)) ], 1-used_range[which.max(r2)])) * (d1<siglim) * (d2<siglim) )
		sum(used_idb)
		#heatscatter(d1_s[used_idb], d2_s[used_idb], log='xy', pch=20)
		#abline(0,1, col='red')
		sf = sum(d2[used_idb]) / sum(d1[used_idb])
		sf_totalmean = sum(d2[]) / sum(d1[])
		sf
		sf_totalmean
	} else{
		used_idb = (d1>100000)*1
		sf = sf_totalmean = sum(d2[]) / sum(d1[])
	}

	### get upper limit siglim
	d1_sf = d1[used_id]
	d1_sf[d1_sf!=siglim] = d1_sf[d1_sf!=siglim] * sf
	d1_sf[d1_sf>=siglim] = siglim
	d2_sf = d2[used_id]

	d1_sf_totalmean = d1[used_id]
	d1_sf_totalmean[d1_sf_totalmean!=siglim] = d1_sf_totalmean[d1_sf_totalmean!=siglim] * sf_totalmean
	d1_sf_totalmean[d1_sf_totalmean>=siglim] = siglim
	d2_sf_totalmean = d2[used_id]

	png(paste(plot_name, '.png', sep=''), width = 1000, height = 1000)
	par(mfrow=c(2,2))
	heatscatter(d1[used_id], d2[used_id], log=use_log_axis, pch=20, main=paste('max r2: ', toString(round(r2[used_r2], digits=3)), '; ', 'quantile_lim: ', toString(round(1-used_range[which.max(r2)], digits=3) ), sep=''))
	abline(0,1, col='red')
	abline(h=d2_thresh, col='blue')
	abline(v=d1_thresh, col='blue')

	plot(seq(1:length(r2)), r2, ylim=c(0.0, 1.01), pch=20)
	abline(h=cor_lim, col='red')
	abline(v=used_r2, col='blue')

	heatscatter(d1_sf, d2_sf, log=use_log_axis, pch=20, main=paste('signorm_sf: ', toString(round(sf, digits=3)), sep=''))
	abline(0,1, col='red')

	heatscatter(d1_sf_totalmean, d2_sf_totalmean, log=use_log_axis, pch=20, main=paste('totalmean_sf: ', toString(round(sf_totalmean, digits=3)), sep=''))
	abline(0,1, col='red')
	dev.off()


	png(paste(plot_name, '.log.png', sep=''), width = 1000, height = 1000)
	par(mfrow=c(2,2))
	heatscatter(d1[used_id], d2[used_id], log='xy', pch=20, xlim=c(0.001, 100), ylim=c(0.001, 100), main=paste('max r2: ', toString(round(r2[used_r2], digits=3)), '; ', 'quantile_lim: ', toString(round(1-used_range[which.max(r2)], digits=3) ), sep=''))
	abline(0,1, col='red')
	abline(h=d2_thresh, col='blue')
	abline(v=d1_thresh, col='blue')

	plot(seq(1:length(r2)), r2, ylim=c(0.0, 1.01), pch=20)
	abline(h=cor_lim, col='red')
	abline(v=used_r2, col='blue')

	heatscatter(d1_sf, d2_sf, log='xy', pch=20, xlim=c(0.001, 100), ylim=c(0.001, 100), main=paste('signorm_sf: ', toString(round(sf, digits=3)), sep=''))
	abline(0,1, col='red')

	heatscatter(d1_sf_totalmean, d2_sf_totalmean, log='xy', pch=20, xlim=c(0.001, 100), ylim=c(0.001, 100), main=paste('totalmean_sf: ', toString(round(sf_totalmean, digits=3)), sep=''))
	abline(0,1, col='red')
	dev.off()

	sf_vector = list("signorm_sf" = sf, "totalmean_sf" = sf_totalmean, 'bg_fg_10' = used_idb)
	return(sf_vector)
}


