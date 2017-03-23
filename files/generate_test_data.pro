
pro generate_fftd
;Creates mock FFT data for use in testing
;Produces an appropriately sized array containing mock data along whistler mode dispersion curve which is Gaussian in k. Spread in omega is controlled by om_fuzz parameter.
;Also produces the derived spectrum from windowed integration along the dispersion curve using calculate_energy_dens routine

;Changeable params
om_fuzz = 5 ;In percent
output_file = "FFT_data.dat"
spectrum_file = "spectrum.dat"

x_ran = [0, 1e6]
t_ran = [0, 5e-2]
n_pts = [1024, 500]

om_ce = 17588.200878
om_pe = 35176.401757
k_peak = n_pts[0]/ 12.0 ;Location of peaks in k relative to k=0 at n_pts[0]/2.0
norm_B_sq = (1e-10)^2 ;100pT field

;Fixed constants
q0 = 1.602176565d-19 ; elementary charge [C]
m0 = 9.10938291d-31  ; electron mass [kg]
v0 = 2.99792458d8    ; speed of light [m/s]
kb = 1.3806488d-23   ; Boltzmann's constant [J/K]
mu0 = 4.0d-7 * !dpi  ; Magnetic constant [N/A^2]
epsilon0 = 8.8541878176203899d-12 ; Electric constant [F/m]
h_planck = 6.62606957d-34 ; Planck constant [J s]


;Code begins
FFT_data = fltarr(n_pts)
sz=size(FFT_data)

;Real units axis resolutions
res = [ abs(x_ran[1]-x_ran[0])/n_pts[0], abs(t_ran[1]-t_ran[0])/n_pts[1] ]

;Make the axes
tmp_ax = make_axis(sz[1], res=res[0])
axes = LIST(tmp_ax, /no_copy)
FOR i=2, sz[0] DO BEGIN
  tmp_ax = make_axis(sz[i], res=res[i-1])
  axes.add, tmp_ax, /no_copy
END

;Calculate whistler dispersion, omega as function of k
dispersion = fltarr(sz[1])
dispersion = axes[0]*axes[0]*float(v0*v0)*om_ce/(axes[0]*axes[0]*float(v0*v0) + om_pe*om_pe)

;Fill the data. Based on calculate_energy_density but in reverse

om_min_ind=intarr(sz[1])
om_max_ind=intarr(sz[1])

FOR i=0, sz[1]-1 DO BEGIN
  om_min=(1.0 - om_fuzz/100.0)*dispersion[i]
  om_max=(1.0 + om_fuzz/100.0)*dispersion[i]
  om_min_ind(i)=where(abs(axes[1]-om_min) eq min(abs(axes[1]-om_min)))
  om_max_ind(i)=where(abs(axes[1]-om_max) eq min(abs(axes[1]-om_max)))
  ;calculate frequency band and corresponding indices
END

k_distrib = exp( - (axes[0] - (axes[0])[n_pts[0]/2.0 + k_peak])^2 / 0.01/max(axes[0])^2 )
k_distrib[n_pts[0]/2.0:n_pts[0]/2.0+10] = 0.0
k_distrib[0:n_pts[0]/2.0 -1] = reverse(k_distrib[n_pts[0]/2.0:n_pts[0]-1])

FOR i=0, sz[1]-1 DO BEGIN
  cells = abs(om_min_ind[i]-om_max_ind[i]-1)
  IF(cells GT 0) THEN FFT_data[i, om_min_ind[i]:om_max_ind[i]] = norm_B_sq * k_distrib[i]/cells
END

;Mirror data
FFT_data[*, 0:n_pts[1]/2-1] = (FFT_data[*, n_pts[1]-1:n_pts[1]/2:-1]) 

;Save the data
err = write_data(output_file, FFT_data, axes, id='FFTd', TIME=[0, 1], SPACE=[0, 1], B_REF = 1.0)

PRINT, "Making spectrum"

;Make spectrum back from data
calculate_energy_density, FFT_data, axes, dispersion, output=spectrum, margin=om_fuzz*2
;Make it fuzzier...
;Save the spectrum
err = write_data(spectrum_file, spectrum, list(dispersion), id="spect", TIME=[0, 1], SPACE=[0, 1], B_REF = 1.0)

end

pro generate_lookup
;Generate the lookup data we're expecting in main code
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  v_therm_par = v0*0.075
  v_therm_perp = v0*0.15
  dens = 0.17
  FORWARD_FUNCTION create_lookup_testdata
  data = create_lookup_testdata(v_therm_par, v_therm_perp, dens, v_lim=0.5, n_els=500)
  err=write_data("./lookupanalytic_bimax.dat",data.data, list(data.ax_x, data.ax_y), id="lookup")
end

function create_lookup_testdata, v_therm_par, v_therm_perp, dens,v_lim=v_lim, n_els=n_els
;Doing this in IDl so I can plot it and not go mad
;Expects v_thermal and produces a bimaxwellian with axes 
;Axes run [- v_lim * c, v_lim*c], and default v_lim = 0.5
;Or supply 2 element v_lim for x and y axes separately
;n_els is number of els in either direction

;TODO add separable version

  COMPILE_OPT IDL2
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  ;Set min and max from supplied limits
  IF(N_ELEMENTS(v_lim) eq 0) THEN v_lim = 0.5
  IF(N_ELEMENTS(v_lim) eq 1) THEN BEGIN
    p_min = - v_lim * v0
    p_max = v_lim * v0
    p_min_y = p_min
    p_max_y = p_max
  ENDIF ELSE BEGIN
    p_min = -v_lim[0]*v0
    p_max = v_lim[0]*v0
    p_min_y= -v_lim[1]*v0
    p_max_y = v_lim[1]*v0
  ENDELSE

  IF(N_ELEMENTS(n_els) EQ 0) THEN n_els = 500
  IF(N_ELEMENTS(n_els) GT 1) THEN BEGIN
    print, 'Excess elements supplied, using n_els='+n_els[0]
    n_els = n_els[0]
  ENDIF

  ;Now create some axes
  ax_x = findgen(n_els)/n_els * float(abs(p_min-p_max)) + float(p_min)
  ax_y = findgen(n_els)/n_els * float(abs(p_min_y-p_max_y)) + float(p_min_y)

  ;Calculate the normalising factors

;  a_par = sqrt(2.0)*v_therm_par/sqrt(1.0 - std::pow(parameters[vpar]/v0, 2))
  a_par = sqrt(2.0)*v_therm_par/sqrt(1.0 - (v_therm_par/v0)^2)

;  a_perp = sqrt( 2.0*std::pow(parameters[vperp], 2)/(1.0 - 2.0*std::pow(parameters[vperp]/v0, 2)));
  a_perp = v_therm_perp*sqrt(2.0)/sqrt(1.0 - 2.0*(v_therm_perp/v0)^2)
  C_norm = 1.0/(!pi*sqrt(!pi)*a_par*a_perp^2)*dens
print, a_par, a_perp, C_norm

  ;Create and populate data array
  data = fltarr(n_els, n_els)
  FOR i=0, n_els-1 DO data[*, i] = exp( - ax_x^2/a_par^2)*exp(-ax_y[i]^2/a_perp^2)
  data[*,*]=C_norm*data[*,*]

  ;Pack and return data & axes
  distrib = {ax_x:ax_x, ax_y: ax_y, data: data}
  return, distrib

end

pro generate_all

  generate_fftd
  generate_lookup

end
