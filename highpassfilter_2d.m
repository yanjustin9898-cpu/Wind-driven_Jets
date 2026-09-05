function blh_p_meanp_highpassed = highpassfilter_2d(blh_p_mean_p, window)

blh_p_mean_lowpassed1=movmean(blh_p_mean_p,window,1,'omitnan');
blh_p_mean_lowpassed2=movmean(blh_p_mean_lowpassed1,window,2,'omitnan');clear blh_p_mean_lowpassed1
blh_p_mean_lowpassed3=movmean(blh_p_mean_p-blh_p_mean_lowpassed2,window,1,'omitnan');
blh_p_mean_lowpassed4=movmean(blh_p_mean_lowpassed3,window,2,'omitnan');clear blh_p_mean_lowpassed3
blh_p_meanp_highpassed=blh_p_mean_p-(blh_p_mean_lowpassed4+blh_p_mean_lowpassed2);clear blh_p_mean_lowpassed2 blh_p_mean_lowpassed4


end