;Note the following require to have run share_omegas!!!
;Growth rates plot
pro plot_growth, growth, _extra=_extra
;plot growth with double axes. Expects growth to be a growth struct containing axes and data
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe

  max_y=10^(ceil(alog10(max(growth.data))))
  plot, growth.axes.x/om_ce, growth.data, ymargin=[4, 4], xtitle='!4x/x!3!Dce!N', ytitle='!4c!3', /ylog, yrange=[max_y-5, max_y], _extra=_extra
  axis, xaxis=1, xtickformat='px_axis', xticks=8, xtickv=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
  xyouts, (!d.x_size /2), !d.y_size-!p.charsize*18, /device, 'p!Dx!N/(m!D0!Nv!D0!N)'
end

;Spectra plots
pro plot_spect_times, spect_arr, _extra=_extra, smth=smth, corr_fac=corr_fac
;Expects an array of spectrum objects, which contain .B.axes and .B.data

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe
  IF(N_ELEMENTS(smth) EQ 0) THEN smth=0
  IF(N_ELEMENTS(corr_fac) EQ 0) THEN corr_fac = 1
  
  max_y=(ceil(alog10(max(spect_arr.B.data)))/2)
  spect_sz=(size(spect_arr))[1]
  cols=findgen(spect_sz)/spect_sz*(255-80) + 80
 
  plot, spect_arr[0].B.axes.x/om_ce, smooth(sqrt(spect_arr[0].B.data)/corr_fac, smth, /edge_wrap), /ylog, yrange=[10l^(max_y-4), 10l^max_y], xtitle='!4x/x!3!Dce!N', ytitle='Wave power', _extra=_extra
  FOR i=0, spect_sz-1 DO oplot, spect_arr[i].B.axes.x/om_ce, smooth(spect_arr[i].B.data, smth, /edge_wrap), color=cols[i]

  tims=fltarr(spect_sz)
  tims = (spect_arr.B.time[0] + spect_arr.B.time[1])/2

  legend_strings = strarr(spect_sz)
  legend_strings = 't='+STRING(format='(F5.2)', tims, /print) + ' s'
  legend, legend_strings, textcolors=cols, /right, /top
  
end

pro df_dp_plot, dist_arr, _extra=_extra, smth=smth, yspan=yspan
;Plot deriv of distribution functions
;Expects array of data arrays each a distribution

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe

  dist_sz=(size(dist_arr))[1]
  max_y=(ceil(alog10(max(deriv(dist_arr.axes.x/m0/v0, dist_arr.data)))))

  IF(N_ELEMENTS(smth) EQ 0) THEN smth=0
  IF(N_ELEMENTS(yspan) EQ 0) THEN yspan=3

  cols=findgen(dist_sz)/dist_sz*(255-80) + 80
  plot, dist_arr[0].axes.x/m0/v0, abs(smooth(deriv(dist_arr[0].axes.x/m0/v0, dist_arr[0].data), smth)), /ylog, xrange=[0, 0.5], /xsty, yrange=[10.0^(max_y-yspan), 10.0^max_y], /ysty,ytitle='df/dp!Dx!N', xtitle='p!Dx!N'
  FOR i=0, dist_sz-1 DO oplot, dist_arr[i].axes.x/m0/v0, abs(smooth(deriv(dist_arr[i].axes.x/m0/v0, dist_arr[i].data), smth)), color=cols[i]
  
  tims=fltarr(dist_sz)
  tims = (dist_arr.time[0] + dist_arr.time[1])/2

;  legend_strings = strarr(dist_sz)
;  legend_strings = 't='+STRING(format='(F5.2)', tims, /print) + ' s'
;  legend, legend_strings, textcolors=cols, /right, /top

end

