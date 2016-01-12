

new_dat = fltarr(dims[0], dims[1])
for i=0, dims[0]-1 DO BEGIN
  new_dat[i,0] = sin(float(i)/25.0)
end
for i =0, dims[1]-1 DO new_dat[*,i] = new_dat[*,0]
;for i =0, dims[1]-1 DO new_dat[*,i] = shift(new_dat[*,0], 10*i)

;generate sine test data

;dat2=complexarr(dims[0], dims[1])
;plan = idlfftw_plan(new_dat, dat2, /VERBOSE)

;idlfftw, plan, new_dat, dat2
;do FFT using fftw dat2 is complex result.

;native_fft = shift(fft(new_dat), [dims[0]/2, dims[1]/2])

;plot, abs(dat2[*,0])^2/(dims[1]*dims[0])^2, xrange=[0, 500]
;plot, abs(native_fft[*,0])^2, xrange=[0, 500]

dat2=complexarr(dims[0])

;dat_in = new_dat[*,0]
dat_in = fltarr(dims[0])

plan = idlfftw_plan(dat_in, dat2, /VERBOSE)

;dat_in = new_dat[*,0]

FOR i=0, dims[0]-1 DO BEGIN
;  IF(dat_in[i] NE new_dat[i,0]) THEN dat_in[i] = new_dat[i,0]
  dat_in[i] = new_dat[i,0]

END

idlfftw, plan, dat_in, dat2
;do FFT using fftw dat2 is complex result.

;native_fft = shift(fft(new_dat[*,0]), dims[0]/2)
native_fft = fft(new_dat[*,0])



;test data checks
;print, data.axes.k[0], -0.00307184272
;print, data.axes.k[10],-0.00305684353
;print, data.axes.k[4095], 0.00307034282
;print, data.axes.omega[0], -3.14159012
;print, data.axes.omega[9], 2.51327205

end
