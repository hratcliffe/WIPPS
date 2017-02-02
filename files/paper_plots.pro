;Scripts to plot things for paper
;Mostly not reuseable

pro plot_all_tests, png=png, eps=eps, legend=legend
;Set eps to create an eps. Else set png to make a PNG DO NOT set both Set legend to add a legend
common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
growth_dir = '/Volumes/Seagate Backup Plus Drive/DiracWhistlers/GrowthTestsPaper/'
names=['fu', 'mod1', 'mod2', 'mod3', 'mod4', 'mod5']

IF(KEYWORD_SET(eps)) THEN BEGIN
  SET_PLOT, 'PS'
  device, filename=growth_dir+'all_growth.eps', xsize=18, ysize=14, bits_per_pixel=64, /encaps, /color
endIF ELSE BEGIN
 !p.background=255
 window, /free
END

n_files =(size(names))[1]
const_1 = read_deck(dir=growth_dir, file='deck_'+names[0]+'.status')
all_consts = list(const_1)
FOR i=1, n_files -1 DO all_consts.add, read_deck(dir=growth_dir, file='deck_'+names[i]+'.status')

growth_1 = read_all_growth(growth_dir+'growth_'+names[0]+'.dat')
all_growth = list(growth_1)
FOR i=1, n_files -1 DO all_growth.add, read_all_growth(growth_dir+'growth_'+names[i]+'.dat')

;colours=[0, 0, 80, 150, 230, 200, 120, 30]  
colours =[0, 0, 40, 80, 150, 210, 230]
leg_text = strarr(n_files+1)
leg_text[0] = 'Run Wpar Wperp Hpar Hperp'

IF(KEYWORD_SET(legend)) THEN xran=[0, 1.5] else xran=[0, 1]

plot, growth_1.axes.x/const_1.wce, growth_1.data/const_1.wce, /ylog, yrange=[1e-5, 0.1], xrange=xran, xtitle='!4x/x!3!Dce!N', ytitle='!4c/x!3!Dce!N', color=colours[1]
leg_text[1] = 'Fu et al '+string(format='(F5.2)',const_1.vtherm_par/v0)+' '+string(format='(F5.2)', const_1.vtherm_perp/v0)+' '+string(format='(F5.2)',const_1.vtherm_parh/v0)+' '+string(format='(F5.2)', const_1.vtherm_perph/v0) 

for i=0, n_files -1 DO BEGIN
  const_1 = all_consts[i]
  growth_1 = all_growth[i]
  leg_text[i+1] = 'M'+string(format='(I1)', i, /print)+' '+string(format='(F5.2)',const_1.vtherm_par/v0)+' '+string(format='(F5.2)', const_1.vtherm_perp/v0)+' '+string(format='(F5.2)',const_1.vtherm_parh/v0)+' '+string(format='(F5.2)', const_1.vtherm_perph/v0) 
  oplot, growth_1.axes.x/const_1.wce, growth_1.data/const_1.wce, color=colours[i+1]
    
END
IF(KEYWORD_SET(legend)) THEN BEGIN
  old_charsize = !p.charsize
  !p.charsize=1.0
  legend, leg_text, textcolors=colours, outline_color=0, /right, /top
  !p.charsize=old_charsize
END

if(keyword_set(png) && ~KEYWORD_SET(eps)) THEN BEGIN
  im=tvrd(true=1)
  details = {inputs:names}
  write_png_l, growth_dir+'all_tests.png', im, details=details
endif else BEGIN
  device, /close
  set_plot, 'X'
  a=plot_logger(growth_dir+'plots.log', {inputs:names ,colors:colours})
END

end
pro plot_time_evol, dirs, outfile=outfile, return_data=return_data, _extra=extr

base_dir = '/Volumes/Seagate Backup Plus Drive/DiracWhistlers/'

IF(N_ELEMENTS(dirs) EQ 0) THEN return
if(N_ELEMENTS(outfile) EQ 0) THEN outfile = base_dir+'time_evol_all.png'

n_dirs=(size(dirs))[1]

times_arr=list(fltarr(get_n_files(dirs[0])))
for i=1, n_dirs-1 DO times_arr.add, fltarr(get_n_files(dirs[i]))

cols = findgen(n_dirs)/float(n_dirs)*200+54
;This is how to copy a list content, rather than get a new reference to it
ens_arr = times_arr[*]
ens_arr_max= ens_arr[*]

for i=0, n_dirs-1 DO BEGIN
 ; while(file_test(base_dir+dirs[i]+'/'+string(format='(I04)', j, /print)+'.sdf')) DO BEGIN
  FOR j=0, (size(times_arr[i]))[1]-1 DO BEGIN
; FOR j=0, 5 DO BEGIN
    if(j mod 10 EQ 0) then print, j
    q=getdata(j, wkdir=base_dir+dirs[i], /ey, /ez)
    times_arr[i, j] = q.time
    ens_arr_max[i, j] = max(gauss_smooth(abs(fft(q.ey^2+q.ez^2)), 1))
    ens_arr[i, j] = total(sqrt(q.ey^2+q.ez^2))
    ;This is IDl list of arrays syntax. first index is list index, second is the array index into that list element. 
  END
END
  !p.background=255
  window, /free
  maxe = 1d-60
  foreach en, ens_arr DO maxe = max([max(en), maxe])
;  maxe=ceil(alog10(maxe))
  plot, times_arr[0], smooth(ens_arr[0], 3), xtitle='Time [s]', ytitle ='Wave Power [arb units]', color=0, yrange=[0, maxe], _extra=extr 
  for i=1, n_dirs-1 DO oplot, times_arr[i], smooth(ens_arr[i], 3), color=cols[i]
  
  im=tvrd(true=1)
  write_png_l, outfile, im, details={smooth:3}
  return_data = {times:times_arr, ens:ens_arr, dirs:dirs}
END

pro plot_2d_spec_from_file, filenames, colours=colours, outfile=outfile, _extra=extr
;Plots 2d spectra using various assumptions to produce E in V/m See below for details

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe
  n_files = (size(filenames))[1]
  
  IF(N_ELEMENTS(outfile) EQ 0) THEN outfile = './Spectra'
  IF(N_ELEMENTS(colours) EQ 0) THEN colours = findgen(n_files-1)/(n_files-1) * 200 + 40
;  window, /free
  SET_PLOT, 'PS'
  device, filename=outfile+'.eps', xsize=16, ysize=12, bits_per_pixel=64, /encaps, /color
  FOR i=0, n_files - 1 DO BEGIN
    spec=read_spect(filenames[i])
    ;Assume 2:1 ratio of n_x to n_y for 2-D. Assume 499 bins too  SHORTCUT
    IF(~N_ELEMENTS(spec)) THEN continue
    corr_fac = spec.B.space[1]^2 / 2.0 * 499.0
    IF(i EQ 0) THEN plot, spec.B.axes.x/om_ce, sqrt(spec.B.data)/corr_fac, /xsty, /ysty, xrange=[-0.8, 0.8], yrange=[0.0001, 0.1], xtitle='!4x/x!3!Dce!N', ytitle='E [V/m]', _extra=extr, /ylog, /nodata
      oplot, spec.B.axes.x/om_ce, sqrt(spec.B.data)/corr_fac,color=colours[i]
  END 
;  im= tvrd(true=1)
;  write_png_l, outfile, im, details={files:filenames, colors:colours}
  device, /close
  set_plot, 'X'
  out_dir = (stregex(outfile, '(.*/+)[^\/]+$', /extract, /subexpr))[1]
  a=plot_logger(out_dir+'plots.log', {inputs:filenames ,colors:colours})
end

pro plot_growth_with_zoom, filenames, colours=colours, outfile=outfile, _extra=extr
;Plot growth rates with and without warm component on linear scale
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe
  tmp = size(filenames)
  if((size(tmp))[1] EQ 4) THEN n_files = (size(filenames))[1] ELSE n_files=1
  
  IF(N_ELEMENTS(outfile) EQ 0) THEN outfile = './GrowAll'
  IF(N_ELEMENTS(colours) EQ 0) THEN colours = findgen(n_files-1)/(n_files-1) * 200 + 40
;  window, /free
  all_grow_b=[]
  all_grow_w=[]
  FOR i=0, n_files - 1 DO BEGIN
    tmp=read_all_growth(filenames[i]+'/grow_plot_an_both.dat')
    if(i EQ 0 && N_ELEMENTS(tmp) GT 0) THEN all_grow_b = [tmp] ELSE all_grow_b = [all_grow_b, [tmp]]  
    tmp=read_all_growth(filenames[i]+'/grow_plot_an_h.dat')
    if(i EQ 0 && N_ELEMENTS(tmp) GT 0) THEN all_grow_w = [tmp] ELSE all_grow_w = [all_grow_w, [tmp]]
  END
 
  SET_PLOT, 'PS'
  device, filename=outfile+'.eps', xsize=16, ysize=12, bits_per_pixel=64, /encaps, /color

  plot, all_grow_b[0].axes.x/om_ce, all_grow_b[0].data/om_ce, /nodata, xrange=[0.1, 1.0], /xsty, yrange=[-0.02, 0.06], /ysty,xtitle='!4x/x!3!Dce!N', ytitle='!4c/x!3!Dce!N'
  FOR i=0, n_files -1 DO oplot, all_grow_b[i].axes.x/om_ce, all_grow_b[i].data/om_ce, color=colours[i]

  device, /close
 
  device, filename=outfile+'zoom.eps', xsize=16, ysize=12, bits_per_pixel=64, /encaps, /color
  plot, all_grow_w[0].axes.x/om_ce, all_grow_w[0].data, /nodata, xrange=[0.6, 1.0], /xsty, yrange=[-0.01, 0.01], /ysty,xtitle='!4x/x!3!Dce!N', ytitle='!4c/x!3!Dce!N'
  FOR i=0, n_files -1 DO oplot, all_grow_b[i].axes.x/om_ce, all_grow_b[i].data/om_ce, color=colours[i]
  FOR i=0, n_files -1 DO oplot, all_grow_w[i].axes.x/om_ce, all_grow_w[i].data/om_ce, color=colours[i], line=2

  device, /close
  set_plot, 'X'
  out_dir = (stregex(outfile, '(.*/+)[^\/]+$', /extract, /subexpr))[1]
  a=plot_logger(out_dir+'plots.log', {inputs:filenames ,colors:colours})
END

pro plot_growth_lin_scale, growths, colours=colours, outfile=outfile, _extra=extr
;Plot some specific growth rates
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe
  tmp = size(growths)
  if((size(tmp))[1] EQ 4) THEN n_lines = (size(growths))[1] ELSE n_lines=1
  
  IF(N_ELEMENTS(outfile) EQ 0) THEN outfile = './GrowLin'
  IF(N_ELEMENTS(colours) EQ 0) THEN colours = findgen(n_lines-1)/(n_lines-1) * 200 + 40
 
  SET_PLOT, 'PS'
  device, filename=outfile+'.eps', xsize=16, ysize=12, bits_per_pixel=64, /encaps, /color

  plot, growths[0].axes.x/om_ce, growths[0].data/om_ce, /nodata, xrange=[0.1, 1.0], /xsty, yrange=[-0.02, 0.06], /ysty,xtitle='!4x/x!3!Dce!N', ytitle='!4c/x!3!Dce!N'
  FOR i=0, n_lines -1 DO oplot, growths[i].axes.x/om_ce, growths[i].data/om_ce, color=colours[i]

  device, /close
 
  set_plot, 'X'
  out_dir = (stregex(outfile, '(.*/+)[^\/]+$', /extract, /subexpr))[1]
  ;This should allow list or array to work...
  inputs = strarr(n_lines)
  FOR i=0, n_lines -1 DO inputs[i] = growths[i].filename
  a=plot_logger(out_dir+'plots.log', {inputs:inputs ,colours:colours})
END

pro plot_pancake_distribs, distribs, nrm=nrm, line_vals=line_vals, line_cols=line_cols, outfile=outfile, _extra=extr
  ;Plot 4 time distributions with overlaid pancake lines
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe
  ;These make sure everything is compiled because IDL is IDL
  RESOLVE_ROUTINE, 'pancake_functions', /IS_FUNCTION, /COMPILE_FULL_FILE
  FORWARD_FUNCTION pancake_curves
  IF(N_ELEMENTS(outfile) EQ 0) THEN outfile = './Disitrbs'
  n_lines = 1
  IF(N_ELEMENTS(line_vals) GT 1) THEN n_lines = (size(line_vals))[1]
  IF(N_ELEMENTS(line_vals) EQ 0) THEN add_lines = 0 ELSE add_lines = 1
  IF(N_ELEMENTS(line_cols) EQ 0) THEN line_cols = findgen(n_lines-1)/(n_lines-1) * 200 + 40
  IF((size(size(distribs)))[1] GT 3) THEN n_distribs = (size(distribs))[1] else n_distribs = 1
  loadct,41,file='colors1.tbl' ;Rainbow from white
  SET_PLOT, 'PS'
  !p.charsize=1.2
  device, filename=outfile+'.eps', xsize=20, ysize=20, bits_per_pixel=64, /encaps, /color
  !p.multi=[0, 2, 2]
  print, n_distribs
  FOR j =0, n_distribs-1 DO BEGIN
    print, j
    contour, alog10(distribs[j].data/nrm), distribs[j].axes.x/m0/v0, distribs[j].axes.y/m0/v0, /iso, nlev=40, /fi, xrange=[-0.5, 0.5], yrange=[-0.6, 0.6], /xsty, /ysty, xtitle='p!Dx!N/(mc)', ytitle='p!Dy!N/(mc)', title='t = '+string(format='(f6.3)', distribs[j].time[0], /print)+' s', color=0, zrange=[-8, -4], /zsty
    IF(add_lines && j GT 0) THEN BEGIN
      FOR i=0, n_lines -1 DO BEGIN
        line = pancake_curves(om_ce/om_pe, line_vals[i])
        oplot, line.v_par/v0, line.v_perp/v0, color=line_cols[i]
        oplot, -line.v_par/v0, line.v_perp/v0, color=line_cols[i]
        oplot, -line.v_par/v0, -line.v_perp/v0, color=line_cols[i]
        oplot, line.v_par/v0, -line.v_perp/v0, color=line_cols[i]
      END 
    END
  END

  device, /close
  !p.charsize=1.5
  set_plot, 'X'
  out_dir = (stregex(outfile, '(.*/+)[^\/]+$', /extract, /subexpr))[1]
  ;This should allow list or array to work...
  a=plot_logger(out_dir+'plots.log', {inputs:distribs.filename ,pancake_lines:line_vals,colours:line_cols})

END


