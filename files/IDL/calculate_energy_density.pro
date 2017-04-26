pro calculate_energy_density, input, axes, dispersion_relation, output=output, margin=margin

;Created by H Ratcliffe, University of Warwick
;16/4/2015
;Vertical integration of input over/pm margin percent around specified dispersion relation
;
;-------------Params----------------
;input: the data, as fft of required field. Must be 2-d array
;axes: the x and y axes of the data
;dispersion relation: curve to integrate along. I.e. omega at each k in supplied k axis in the units used in supplied axes. 
;output: returned energy density as function of x
;margin: percentage margin around dispersion curve. Default /pm 1%
;
;Call example:
;calculate_energy_density, ex_fft, axes, sqrt(1.+3.*axes[0]^2), output=output, margin=5


IF(ARG_PRESENT(output) EQ 0) THEN BEGIN
  print, 'No output array reference provided, returning'
  return
ENDIF
;check for return array by reference

if(N_ELEMENTS(margin) EQ 0) THEN margin = 0.01 ELSE margin = float(margin)/100.

sz = SIZE(input)

if(sz[0] NE 2) THEN BEGIN
  print, 'Input must be 2-d array, returning'
  return
ENDIF

IF(N_ELEMENTS(dispersion_relation) NE sz[1]) THEN BEGIN
  print, 'Dispersion relation size inconsistent, returning'
  return
ENDIF


om_min_ind=intarr(sz[1])
om_max_ind=intarr(sz[1])

FOR i=0, sz[1]-1 DO BEGIN
  om_min=(1-margin)*dispersion_relation[i]
  om_max=(1+margin)*dispersion_relation[i]

;  om_min_ind(i)=where(abs(axes[1]-om_min) eq min(abs(axes[1]-om_min)))
;  om_max_ind(i)=where(abs(axes[1]-om_max) eq min(abs(axes[1]-om_max)))
  om_min_ind[i]=min(where((axes[1] GT om_min)))
  om_max_ind[i]=min(where((axes[1] GT om_max)))
  ;calculate frequency band and corresponding indices
END

mask_arr=fltarr(sz[1:2])
output=fltarr(sz[1])

FOR i=0, sz[1]-1 DO BEGIN
  mask_arr[i, om_min_ind[i]:om_max_ind[i]] = 1. 
  output[i]=total(mask_arr[i,*] * input[i,*])
END

END