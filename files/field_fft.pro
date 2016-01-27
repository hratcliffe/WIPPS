
pro field_fft, input, output=output, axes=axes, res=res

;Created by H Ratcliffe, University of Warwick
;16/4/2015
;Does FFT of supplied field (minimising temporary creations etc)
;
;-------------Params----------------
;input: the array to transform
;output: returned fft
;axes: returned axes for each dimension
;res: array of axis resolutions in normed units (see make_axis). If not specified/wrong dimension then 1 is used for each
;Call example:
; field_fft, ex, output=ex_fft, axes=axes
;and for a 2-d input, e.g.:
; contour, alog(ex_fft),axes[0], axes[1], nlev=40, /fi, xrange=[-1, 1], yrange=[0, 2]

sz = SIZE(input)
n_dims=sz[0]

if(n_dims EQ 0) THEN BEGIN
  print, 'Input has zero dimension'
  return
ENDIF

if(n_dims GT 4) THEN BEGIN
  print, 'Array dimension too large'
  return
ENDIF

if(N_ELEMENTS(res) GT 0) THEN IF(SIZE(res, /dim) NE n_dims) THEN BEGIN
  print, 'WARNING: res dimension does not match input dimension. Using ones'
  res=replicate(1, n_dims)
ENDIF 
if(N_ELEMENTS(res) EQ 0) THEN res=replicate(1, n_dims)

;check dimensions of array and of params

IF(ARG_PRESENT(output) EQ 0) THEN BEGIN
  print, 'No output array reference provided, returning'
  return
ENDIF
;check for return array by reference

print, 'Performing fft'
output = fft(input)
output = (temporary(output))^2
output = abs(output)

;here the fft changes type to complex so temporaries can do nothing. The squaring does not type change and temporary therefore helps. Abs again changes type. In simple test, saving of approx 2*size(input) in step 2 due to not creating 2 copies with complex type

shft = intarr(sz[0])
for i=1, sz[0] DO shft[i-1] = sz[i]/2
output = shift(output, shft)
;shift fft 

IF(ARG_PRESENT(axes) EQ 1) THEN BEGIN
  print, 'Generating axes'
  tmp_ax = make_axis(sz[1], res=res[0])
  axes = LIST(tmp_ax, /no_copy)
  for i=2, sz[0] DO BEGIN
    tmp_ax = make_axis(sz[i], res=res[i-1])
    axes.add, tmp_ax, /no_copy
  END
  ;create axes as list (ragged array)
END

END