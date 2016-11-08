pro angular_distribs, spec_in, freqs=freqs, om_ce=om_ce, _extra = extr
;Plot the relative angular distribution for given freqs. Legend shows the frequency and the relative peak power at that freq
;Extra keyword params can be given and will be passed through to plot and oplot, e.g. line=, psym=, xrange= etc 

n_freqs=(size(freqs))[1]
if(n_elements(om_ce) EQ 0) THEN om_ce = 1

max_powers = findgen(n_freqs)
for i=0, n_freqs-1 DO max_powers[i] = spec_in.B.data[freqs[i]]
max_powers = max_powers/min(max_powers)
colours=50+findgen(n_freqs)*205/n_freqs

freq_vals=spec_in.B.axes.x[freqs]/om_ce
freq_strings = string(format='(F6.2)', freq_vals, /print)
power_strings = string(format='(F6.2)', max_powers, /print)
legend_strings = freq_strings + replicate(' ', n_freqs)+power_strings
legend_strings = ["f/f_pe  P/P_0", legend_strings]

ang=atan(spec_in.ang.axes.y)*180/!pi
max_y = ceil(max(spec_in.ang.data[freqs, *])*10.0)/10.0
;Set max axis to nearest 10th

plot, ang, spec_in.ang.data[freqs[0], *], xrange=[0, 30], xtitle='Angle (deg)', ytitle='Relative power', yrange = [0, max_y], /nodata, _extra=extr
for i=0, n_freqs-1 DO oplot, ang, spec_in.ang.data[freqs[i], *], color=colours[i], _extra=extr

legend, legend_strings, textcolors=[255, colours], /right, /top, box=0, charsize=1.2

end
