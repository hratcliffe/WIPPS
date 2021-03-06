
pro generate_fftd
;Creates mock FFT data for use in testing
;Produces an appropriately sized array containing mock data along whistler mode dispersion curve which is Gaussian in k. Spread in omega is controlled by om_fuzz parameter.
;Also produces the derived spectrum from windowed integration along the dispersion curve using calculate_energy_dens routine

common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
common extra_consts, global_file_dir

;Changeable params
om_fuzz = 5 ;In percent
output_file = global_file_dir + "FFT_data.dat"
spectrum_file = global_file_dir + "spectrum.dat"

x_ran = [0, 1e6]
t_ran = [0, 5e-2]
n_pts = [1024, 500]

om_ce = 17588.200878
om_pe = 35176.401757
k_peak = n_pts[0]/ 12.0 ;Location of peaks in k relative to k=0 at n_pts[0]/2.0
norm_B_sq = (1e-10)^2/1e4 ;100pT field over 1e4 "frequency width"


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

pro generate_fitted_spect
;Creates a spectrum fitted to 2DDualMatch data at first time

om_info = {mid:0.53, delta:0.06, min:0.2, max:0.8, n_pts:1024, norm:1.5e-19}
x_info = {mid:0, delta:tan(!pi/6), min:0, max:1, n_pts:100}

inputdir = '/Volumes/Seagate Backup Plus Drive/DiracWhistlers/2DDualMatch/'
const = read_deck(dir = inputdir)

generate_mock_spect, om_info, x_info, const, outfile=inputdir+"match_spectrum.dat"

end

pro generate_intermediate_spect
;Creates a spectrum intermediate to the Albert and fitted

om_info = {mid:0.45, delta:0.1, min:0.2, max:0.7, n_pts:1024, norm:1.5e-19}
x_info = {mid:0, delta:tan(!pi/6), min:0, max:1, n_pts:100}

inputdir = '/Volumes/Seagate Backup Plus Drive/DiracWhistlers/2DDualMatch/'
const = read_deck(dir = inputdir)

generate_mock_spect, om_info, x_info, const, outfile=inputdir+"inter_spectrum.dat"

end


pro generate_Albert_spect
;Creates mock spectrum as used in Albert for example D calcs
;Truncated Gaussians in B and angle

;Changeable params
om_info = {mid:0.35, delta:0.15, min:0.05, max:0.65, n_pts:1024, norm:(1e-10)^2}
;100pT field, so if int B d omega is normed to 1 this will be correct

x_info = {mid:0, delta:tan(!pi/6), min:0, max:1, n_pts:100}
;om_mid = 0.35
;delta_om = 0.15
;om_min = 0.05
;om_max = 0.65
;delta_x = tan(!pi/6)
;x_min = 0
;x_mid = 0
;x_max = 1

;n_pts_om = 1024
;n_pts_x = 100

deck_file_prefix = "d_test"

const = read_deck(dir = global_file_dir, pref = deck_file_prefix)

generate_mock_spect, om_info, x_info, const, outfile="d_testspectrum.dat"

end

pro generate_mock_spect, om, x, const, outfile=outfile
;Creates mock spectrum of truncated Gaussians in B and angle
;Takes two structs each containing mid, delta, min, max, n_pts and (for omega only) norm one for omega, one for x, and a deck-constants structure
;Optionally also takes an output file path

common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
common extra_consts, global_file_dir

IF(N_ELEMENTS(outfile) GT 0) THEN spectrum_file = outfile ELSE spectrum_file = "mock_spectrum.dat"

om_ce = const.wce
om_pe = const.wpe
desired_norm_B_sq = om.norm
norm_B = 1

;Make the axes

om_ax = findgen(om.n_pts)*om.max*1.2/om.n_pts
x_ax = findgen(x.n_pts)*x.max*1.2/x.n_pts

B_arr = exp(- (om_ax - om.mid)^2/om.delta^2)
B_arr(where(om_ax LT om.min)) = 0.0
B_arr(where(om_ax GT om.max)) = 0.0
norm_B = float(total(B_arr)*abs(om_ax[1]-om_ax[0])*om_ce)
B_arr = desired_norm_B_sq*B_arr/norm_B


ang_arr = transpose(exp(-(x_ax - x.mid)^2/x.delta))
ang_arr(where(x_ax LT x.min)) = 0.0
ang_arr(where(x_ax GT x.max)) = 0.0

B = {data:B_arr, axes:{x:om_ax*float(const.wce)}, space:[0, 1], time:[0.0, 1.0], B_ref: 1.0}
ang = {data:ang_arr, axes:{x:[1.0], y:x_ax}, space:[0, 1], time:[0.0, 1.0], B_ref: 1.0}
spectrum = {B: B, ang: ang}

err = write_spect(spectrum_file, spectrum)

end

pro generate_lookup
;Generate the lookup data we're expecting in main code

common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
common extra_consts, global_file_dir

  v_therm_par = v0*0.075
  v_therm_perp = v0*0.15
  dens = 0.17
  FORWARD_FUNCTION create_lookup_testdata
  data = create_lookup_testdata(v_therm_par, v_therm_perp, dens, v_lim=0.5, n_els=500)
  err=write_data(global_file_dir + "lookupanalytic_bimax.dat",data.data, list(data.ax_x, data.ax_y), id="lookup")
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
  common extra_consts, global_file_dir

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

pro reprocess_padie_data
;Convert the Padie test data into my-style data array files

;read_padie_data, dir, file_start = file_start, n = n
;function write_data, filelabel, data, axes, usenum=usenum, _extra=extr

common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
common extra_consts, global_file_dir

  padie_dir = global_file_dir + "PadieTestData/"
  data = read_padie_data(padie_dir)
  outfile_name = "Data.dat"
  err = write_data(padie_dir + outfile_name, float(data.daa), list(float(data.axis)))
  outfile_name = "Data_p.dat"
  err = write_data(padie_dir + outfile_name, float(data.dpp), list(float(data.axis)))

  FOR i = -2, 2 DO BEGIN
    outfile_name = "Data_"
    if i LT 0 THEN outfile_name = outfile_name + "-"
    outfile_name = outfile_name + string(format="(I1)", abs(i), /print)
    data = read_padie_data(padie_dir, n=i)
    print, "Writing " + outfile_name
    err = write_data(padie_dir + outfile_name + ".dat", float(data.daa), list(float(data.axis)))

    err = write_data(padie_dir + outfile_name + "_p.dat", float(data.dpp), list(float(data.axis)))


  END

end

pro generate_all

  generate_fftd
  generate_lookup

end
