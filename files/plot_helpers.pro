pro define_consts
;  common 
consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  
  q0 = 1.602176565d-19 ; elementary charge [C]
  m0 = 9.10938291d-31  ; electron mass [kg]
  v0 = 2.99792458d8    ; speed of light [m/s]
  kb = 1.3806488d-23   ; Boltzmann's constant [J/K]
  mu0 = 4.0d-7 * !dpi  ; Magnetic constant [N/A^2]
  epsilon0 = 8.8541878176203899d-12 ; Electric constant [F/m]
  h_planck = 6.62606957d-34 ; Planck constant [J s]

end

pro share_omegas, om_ce, om_pe
  common omega_share, om_ce_sh, om_pe_sh
  om_ce_sh = om_ce
  om_pe_sh = om_pe
end

function plot_logger, filename, details
;Stringify the details into a blob and write date/time and details to plots.log in directory of filename
;Generally either supply a struct in which case the tag names are included in strings, or pre-prepared string array
;E.g. details = {routine:'plot_fft', source:'data/run1/time0.dat', smooth: 1}
  dir = strjoin(((strsplit(filename, '/', /extract))[0:-2]), '/')+'/'
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
  details_print = time+' '
  if(isa(details, 'struct')) THEN BEGIN
    tag_nams = tag_names(details)
    for i=0, (size(tag_nams))[1]-1 DO details_print += tag_nams[i]+' '+strtrim(string(details.(i), /print), 2)+' '
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

  
