;pro read_data
;Written for commit ID from b8b80f3 to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

my_type_code = 'f'
my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

id_type ='12345678910'
id_in = id_type
;type for block id...

io_type =0.0
io_check = 3.0/32.0
io_in=io_type
;type for io verification const.

commit_type ='123456789112345'
commit_in = commit_type
;type for commit id

int_type = 1

openr, 1, filename
;"Tmp.txt"
;open file

readu, 1, id_in
readu, 1, io_in
readu, 1, commit_in

n_dims = int_type
readu, 1, n_dims

if(n_dims GT 0) THEN BEGIN
  dims = lonarr(n_dims)
  readu, 1, dims
  IF my_type_code EQ 'f' THEN BEGIN
    axes_list = {k:fltarr(dims[0]), omega:fltarr(dims[1])}
    data = {id:id_in,data:fltarr(dims[0], dims[1]), axes:axes_list}
  ENDIF ELSE BEGIN
    axes_list = {k:dblarr(dims[0]), omega:dblarr(dims[1])}
    data = {id:id_in,data:dblarr(dims[0], dims[1]), axes:axes_list}
  ENDELSE

  ;Now we do need the right majority...

  tmp2=fltarr(dims[0], dims[1])
  ;Can't read directly into anon structure field, so use tmp
  readu, 1, tmp2
  data.data = tmp2
  tmp2=0
ENDIF ELSE BEGIN
  print, "Array is ragged. Use read_ragged.pro"

ENDELSE

;EXPLICTLY 2-D here
tmpa=fltarr(dims[0])
;Can't read directly into anon structure field, so use tmp
readu, 1, tmpa
data.axes.k = tmpa

tmpa=fltarr(dims[1])
readu, 1, tmpa
data.axes.omega = tmpa

tmpa = 0

close, 1

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
