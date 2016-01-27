pro generate_fftd
;Creates mock FFT data for use in testing
;Produces an appropriately sized array containing mock data along whistler mode dispersion curve which is Gaussian in k. Spread in omega is controlled by om_fuzz parameter.
;Also produces the derived spectrum from windowed integration along the dispersion curve using calculate_energy_dens routine

;Changeable params
om_fuzz = 5 ;In percent
output_file = "FFT_data.dat"
spectrum_file = "spectrum.dat"

x_ran = [0, 1e7]
t_ran = [0, 5e-2]
n_pts = [1024, 500]

om_ce = 17588.200878
om_pe = 35176.401757
k_peak = n_pts[0]/ 8.0 ;Location of peaks in k relative to k=0 at n_pts[0]/2.0

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
dispersion = axes[0]*axes[0]*v0*v0*om_ce/(axes[0]*axes[0]*v0*v0 + om_pe*om_pe);

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
k_distrib[0:n_pts[0]/2.0 -1] = reverse(k_distrib[n_pts[0]/2.0:n_pts[0]-1])

FOR i=0, sz[1]-1 DO BEGIN
  cells = abs(om_min_ind[i]-om_max_ind[i]-1)
  IF(cells GT 0) THEN FFT_data[i, om_min_ind[i]:om_max_ind[i]] = k_distrib[i]/cells
END


;Save the data
err = write_data(output_file, FFT_data, axes, id='FFTd')

;Make spectrum back from data
calculate_energy_density, FFT_data, axes, dispersion, output=spectrum, margin=om_fuzz

;Save the spectrum
err = write_data(spectrum_file, spectrum, axes, id="spect")

end