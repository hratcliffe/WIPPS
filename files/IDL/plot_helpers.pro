pro define_consts
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  
  q0 = 1.602176565d-19 ; elementary charge [C]
  m0 = 9.10938291d-31  ; electron mass [kg]
  v0 = 2.99792458d8    ; speed of light [m/s]
  kb = 1.3806488d-23   ; Boltzmann`s constant [J/K]
  mu0 = 4.0d-7 * !dpi  ; Magnetic constant [N/A^2]
  epsilon0 = 8.8541878176203899d-12 ; Electric constant [F/m]
  h_planck = 6.62606957d-34 ; Planck constant [J s]

end

pro share_omegas, om_ce, om_pe
  common omega_share, om_ce_sh, om_pe_sh
  om_ce_sh = om_ce
  om_pe_sh = om_pe
end

function plot_logger, filename, details, stack_levels=stack_levels
;Stringify the details into a blob and write date/time, call stack and details to plots.log in directory of filename
;Generally either supply a struct in which case the tag names are included in strings, or pre-prepared string array
;E.g. details = {source:'data/run1/time0.dat', smooth: 1}
;set stack_levels to determine max depth of calls. Default 3 levels. Number does NOT include this function
  dir = strjoin(((strsplit(filename, '/', /extract))[0:-2]), '/')+'/'
  if(strmid(filename, 0, 1) EQ '/') THEN dir = '/'+dir
  filename_part = (strsplit(filename, '/', /extract))[-1]
 ;Check this directory exists
  if(~file_test(dir)) THEN BEGIN
    print, 'plot_logger: no such directory '+dir
    RETURN, 1
  endif
  ;grab time for log
  time=systime()
  if(~isa(details, 'struct') && (~isa(details[0], 'string'))) THEN BEGIN
    PRINT, 'Invalid details, supply struct, or [array of] strings'
    RETURN, 1
  ENDIF
  details_print = time+' '+filename_part+' '
  if(N_ELEMENTS(stack_levels) EQ 0) THEN stack_levels = 3
  call_stack = ((scope_traceback(/structure)).routine)
  if((size(call_stack))[1] GE stack_levels+1) THEN call_stack = call_stack[-stack_levels-1: -1]
  ;Drop this level
  call_stack = call_stack[0:-2]
  call_stack = strjoin(call_stack, '->')
  details_print += call_stack+' '
  if(isa(details, 'struct')) THEN BEGIN
    tag_nams = tag_names(details)
    for i=0, (size(tag_nams))[1]-1 DO details_print += tag_nams[i]+' '+strjoin(strtrim(string(details.(i), /print), 2), ', ')+'; ' 
  endif else if(isa(details, /arr) && isa(details[0], 'string')) THEN BEGIN
      details_print = strjoin(details, ' ')
    endif
 ; end
  if(file_test(dir+'plots.log')) THEN BEGIN 
    openu, filenum, dir+'plots.log', /get_lun, /append
  endif else BEGIN
    openw, filenum, dir+'plots.log', /get_lun
  end
  printf, filenum, details_print
  free_lun, filenum
  RETURN, 0
end

pro write_png_l, filename, im, details=details
;Write image to png 

write_png, filename, im
if(n_elements(details) NE 0) THEN p=plot_logger(filename, details)

end

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

function average_pitch_angle, T_par=T_par, A=A, energy=energy, distrib=distrib, n_angs=n_angs
;Calculate the average pitch angle of particles at a given energy assuming parallel temperature T_par and Anisotropy of A. Distrib can be (M)axwellian or (P)ancake
;Return value should always be between 0 and 90 (degrees)
;Energy should be in keV, as should T_par
;n_angs is the number of angles to use in calculation, 500 by default

common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

IF(N_ELEMENTS(T_par) EQ 0 OR N_ELEMENTS(A) EQ 0) THEN BEGIN
  PRINT, "Distribution underspecified, require T_par and A"
  RETURN, -1
END

IF(N_ELEMENTS(n_angs) EQ 0) THEN n_angs = 500

;Get first character of distrib string
distrib = strmid(strupcase(distrib), 0, 1)

;Get energy in joules and then total momentum from energy input
en_Joules = q0*energy*1e3
mod_p = sqrt(2.0*m0*en_Joules)

;Get constants for distribution function
if(distrib EQ 'M') THEN BEGIN
  ;one /m0 from defintion of v_te and a square to convert v to p
  a_par_sq = m0*T_par*1e3*q0*2.0/sqrt(1.0 - (T_par*q0/m0)/v0/v0)
  a_perp_sq = m0*(A+1.0)*T_par*1e3 *q0*2.0/sqrt(1.0 - ((A+1.0)*T_par*q0/m0)/v0/v0)
ENDIF ELSE BEGIN
  if(distrib EQ 'P') THEN BEGIN
    PRINT, "TBC!!!"
    RETURN, -1
  
  ENDIF
ENDELSE

tot = 0.0
tot_f = 0.0
d_theta = 90.0/n_angs

FOR i=0, n_angs-1 DO BEGIN
  ang = i*d_theta*!pi/180.0
  p_x = mod_p*cos(ang)
  p_y = mod_p*sin(ang)
  ;print, p_x, p_y,  a_par_sq, a_perp_sq, exp(-p_x^2/a_par_sq)*exp(-p_y^2/a_perp_sq)
  if(distrib EQ 'M') THEN BEGIN
    tot = tot + ang* exp(-p_x^2/a_par_sq)*exp(-p_y^2/a_perp_sq)
    tot_f = tot_f + exp(-p_x^2/a_par_sq)*exp(-p_y^2/a_perp_sq)
  ENDIF ELSE BEGIN 
    IF(distrib EQ 'P') THEN BEGIN
  
    ENDIF ELSE PRINT, "Unknown distribution requested"
  
  ENDELSE
END

;Convert to degrees
tot = tot/tot_f*180.0/!pi

return, tot

end

function px_axis, axis, index, value
;Return string formatted resonant_px(omega)
;Run share_omegas first to populate om_ce, om_pe
;print, axis, index, value
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe
  omega=value*om_ce
  k = get_dispersion(omega, om_ce=om_ce, om_pe=om_pe, /k, /appr)
  p_x = (omega - om_ce)/k
  str = string(FORMAT='(F5.2)',p_x/v0)
  if(p_x NE p_x) THEN str=""
  if(p_x EQ 'Inf' OR p_x EQ '-Inf') THEN str=""
  ;Blank out an Inf or Nan string
  return, str

end


function get_n_files, dir, lead=lead, ext=ext
;Count files in directory
  if(N_ELEMENTS(lead) EQ 0) THEN lead = '/Volumes/Seagate Backup Plus Drive/DiracWhistlers/'
  if(N_ELEMENTS(ext) EQ 0) THEN ext='*.sdf'
  a=file_search(lead+dir+'/'+ext, count=cnt)
  return, cnt
end

function binary_invert, val, func_name, precision=precision, _extra=ext
;Locate val in func_name. Assumes monotonic and non-constant on args from 0 to 1

;Check increasing or decreasing
  grad = 1
  IF( call_function(func_name,0.2, _extra=ext) LT call_function(func_name,0.1, _extra=ext)) THEN grad = -1

  ;Fractional precision to terminate
  IF(N_ELEMENTS(precision) EQ 0) THEN precision = 0.00001d0
  selection = 0.5d0
  tmp = call_function(func_name,selection, _extra=ext)
  current_increment = 0.25d0
  counter = 0
  max_it = 50 ; 0.5^50 = 8e-16
  while(abs((tmp - val)) GT precision AND counter LT max_it) DO BEGIN
   tmp = call_function(func_name,selection, _extra=ext)
    if(tmp GT val) THEN selection = selection - grad*current_increment
    if(tmp LT val) THEN selection = selection + grad*current_increment
    if(tmp EQ val) THEN break
    current_increment = current_increment * 0.5
    counter = counter + 1
  END

  IF(counter GT max_it -1) THEN return, -1

  return, selection

end

function get_p_res, om_r, om_ce=om_ce, om_pe=om_pe, n=n
;Convert from omega_res to p_res using gamma v_par = (n*Om_ce - om_r)/k_res, par, assuming parallel propagation and n=-1
;NOTE om_r is om/om_ce and v is returned as v/v0
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce_imp, om_pe_imp

  ;If om_ce passed use it, else if the common block is defined use those, else use defaults
  IF(N_ELEMENTS(om_ce) EQ 0) THEN IF(N_ELEMENTS(om_ce_imp) EQ 0) THEN om_ce=2.0*!pi*10000.0 else om_ce = om_ce_imp
  IF(N_ELEMENTS(om_pe) EQ 0) THEN IF(N_ELEMENTS(om_pe_imp) EQ 0) THEN om_pe = 3.0*om_ce else om_pe=om_pe_imp
  IF(N_ELEMENTS(n) EQ 0) THEN n= -1
  k_tmp = get_dispersion(om_r*om_ce, om_ce=om_ce, om_pe=om_pe, /k)
  gv_old=abs(om_ce*(1.0*n + om_r)/k_tmp/v0)
  ;Calc gamma from this v and do a second iteration
  gamm_prev = 1.0
  precision = 0.0001
  MAX_IT = 20
  gamm = sqrt(1.0+gv_old^2)
  i=0
  gv_new = gv_old
  WHILE((abs(gamm_prev/gamm - 1.0) GT precision) && i LT MAX_IT) DO BEGIN
 ;   gamm = sqrt(1.0+gv_old^2)
    gv_new = abs(om_ce*(1.0*n + om_r*gamm)/k_tmp/v0)
    ;print, gv_1, gamm, gv_2
    gamm_prev = gamm
    gamm = sqrt(1.0+gv_new^2)
    gv_old=gv_new
   i=i+1
  END
  IF(i EQ MAX_IT) THEN print, "Exceeded max iterations!!!!"
  return, gv_new
end

function get_om_res, v_r, om_ce =om_ce, om_pe=om_pe, n=n
;Convert from v_res to om_res using gamma v_par = (n*Om_ce - om_r)/k_res, par, assuming parallel propagation and n=-1
;NOTE v_r (p) is expected as v/v0 and om is returned as om/om_ce
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce_imp, om_pe_imp

  ;If om_ce passed use it, else if the common block is defined use those, else use defaults
  IF(N_ELEMENTS(om_ce) EQ 0) THEN IF(N_ELEMENTS(om_ce_imp) EQ 0) THEN om_ce=2.0*!pi*10000.0 else om_ce = om_ce_imp
  IF(N_ELEMENTS(om_pe) EQ 0) THEN IF(N_ELEMENTS(om_pe_imp) EQ 0) THEN om_pe = 3.0*om_ce else om_pe=om_pe_imp
  IF(N_ELEMENTS(n) EQ 0) THEN n= -1
  
  ;Use basic binary interative inversion. Works as long as omega(k) strictly monotonic, with number of iterations depending on the slope
  ;I.e. we hunt for the v_res input to get out the wanted value
  om=binary_invert(v_r, 'get_p_res', om_ce=om_ce, om_pe=om_pe, n=n)
  return, om
end

function om_axis, axis, index, value
;Return string formatted resonant_px(omega)
;Run share_omegas first to populate om_ce, om_pe
;print, axis, index, value
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  common omega_share, om_ce, om_pe
  v_r=value
  om=get_om_res(v_r, om_ce=om_ce, om_pe=om_pe)
  str = string(FORMAT='(F5.2)',om)
  if(om NE om) THEN str=""
  if(om EQ 'Inf' OR om EQ '-Inf') THEN str=""
  ;Blank out an Inf or Nan string
  return, str

end

function trace_energy, v_ax, ang_ax, energy
;Trace the index pairs into v and ang axes which map to a total electron energy of energy. For each v, the corresponding alpha index is found

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

  ;Work off v axis
  sz = (size(v_ax))[1]
  index_arr = intarr(sz)
  ;Map energy to v
  gamma = energy/m0/v0^2 + 1 ; 1 is for rest-mass energy
  total_v = v0 * sqrt(gamma^2 - 1.0)/gamma
  for i = 0, sz-1 DO BEGIN
    angle = acos(abs(v_ax[i])/total_v)
    index_arr[i] = min(where(ang_ax GT angle))
  end

  return, index_arr
end

function gamma, v
;Takes v/c
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  return, 1.0/sqrt(1.0 - v^2)
end

function gamma_from_p, p
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  return, sqrt(1.0 + (p/v0/m0)^2)
end

function gamma_from_energy, energy
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

  return, energy/m0/v0^2 + 1 ; 1 is for rest-mass energy
end

function keV_to_v, energy
  ;takes energy in keV and returns velocity in v0/c
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

  gam = gamma_from_energy(energy*1e3*q0)
  return, sqrt(gam^2 -1)/gam
end

function v_to_keV, v
  ;takes velocity in v0/c and returns energy in keV
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

  gam = gamma(v)
  return, (gam - 1)*m0*v0^2/q0/1e3
end


function get_D_at_energy, D, energy, velocity=velocity, spread=spread
;Lineout D to get rms(D) as function of angle (or velocity if keyword velocity set) at given energy
;Spread is percent energy band to use, default ±10% is spread=10. Max is 99%

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

  ;Set default spread and re-range
  if(N_ELEMENTS(spread) EQ 0) THEN spread = 10
  if(spread LT 0) THEN spread = abs(spread)
  if(spread GT 99) THEN spread = 99

  n_els = (size(D.data))[2]
  other_els = (size(D.data))[1]
  if(KEYWORD_SET(velocity)) THEN BEGIN
    n_els = (size(D.data))[1]
    other_els = (size(D.data))[2]
  END

  inds_low = trace_energy(D.axes.x, D.axes.y, energy*(1.0-spread/100.0))
  inds_high = trace_energy(D.axes.x, D.axes.y, energy*(1.0+spread/100.0))
  inds_high[where(inds_high EQ -1)] = other_els - 1
  inds_low[where(inds_low EQ -1)] = other_els - 1
  inds_low[where(inds_low GT inds_high)] = 0

  D_masked = fltarr((size(D.data))[1:2])
  for i = 0, n_els - 1 DO D_masked[i, inds_low[i]:inds_high[i]] = D.data[i, inds_low[i]:inds_high[i]]/(inds_high[i] - inds_low[i] + 1)

  if(KEYWORD_SET(velocity)) THEN D_line = total(abs(D_masked^2), 2)/total(D_masked^2/(D_masked^2+1e-30), 2) else D_line = total(abs(D_masked^2), 1)/total(D_masked^2/(D_masked^2+1e-30), 1)

  return, sqrt(D_line)

end

