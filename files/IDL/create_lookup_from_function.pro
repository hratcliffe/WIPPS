pro create_lookup_from_function, n_x, n_y, function_name, p_x_lims=p_x_lims, p_y_lims=p_y_lims, outfile = outfile, _extra=_extra
;Create a lookup table distribution of size n_x * n_y, with axes running from p_{x,y}_lims[0] to p_{x,y}_lims[1] if supplied, with each element populated by calling function_name (given as a string) with (p_x, p_y) and any additional keyword params passed in _extra
;The function should be suitably normed etc
;Result is written in outfile, which defaults to 'lookup.dat'
;Growth code expects linear axes in v


COMPILE_OPT IDL2
;force long ints and proper brackets

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

  IF(N_ELEMENTS(p_x_lims) NE 2) THEN p_x_lims = [-2e-22, 2e-22] ;0.6 c electron
  IF(N_ELEMENTS(p_y_lims) NE 2) THEN p_y_lims = [-2e-22, 2e-22]
  IF(N_ELEMENTS(outfile) EQ 0) THEN outfile = 'lookup.dat'

  v_x = float(findgen(n_x)*(p_x_lims[1]/m0/gamma_from_p(p_x_lims[1]) - p_x_lims[0]/m0/gamma_from_p(p_x_lims[0]))/(n_x-1) + p_x_lims[0]/m0/gamma_from_p(p_x_lims[0]))
  v_y = float(findgen(n_y)*(p_y_lims[1]/m0/gamma_from_p(p_y_lims[1]) - p_y_lims[0]/m0/gamma_from_p(p_y_lims[0]))/(n_y-1) + p_y_lims[0]/m0/gamma_from_p(p_y_lims[0]))

  lookup_data = fltarr(n_x, n_y)

  FOR ix = 0, n_x - 1 DO BEGIN
    FOR iy = 0, n_y -1 DO BEGIN
      lookup_data[ix, iy] = call_function(function_name, v_x[ix]*m0*gamma(v_x[ix]/v0), v_y[iy]*m0*gamma(v_y[iy]/v0), _extra = _extra)
    END
  END

  ;Force types in case function_name did not return floats
  lookup_data = float(lookup_data)

  err = write_data(outfile, lookup_data, list(v_x, v_y), time = [0.0, 1.0], space = [0, 1], B_ref = 1.0, block = "IDLlookup")

  if(err.err) THEN PRINT, "Error writing file"

  ;Print log of params. Allow one layer of structs
  openw, filenum, outfile+'.pars', /GET_LUN
  printf, filenum, 'Function ' + function_name
  if(N_ELEMENTS(_extra) GT 0) THEN tags = TAG_NAMES(_extra) else tags = []
  FOR i=0, (size(tags))[1]-1 DO BEGIN
    IF(~ISA(_extra.(i), 'STRUCT')) THEN BEGIN
      printf, filenum, tags[i], _extra.(i)
    ENDIF ELSE BEGIN
      printf, filenum, tags[i], '{'
      sub_tags = tag_names(_extra.(i))
      FOR j=0, (size(sub_tags))[1]-1 DO printf, filenum, sub_tags[j], _extra.(i).(j)
      printf, filenum, '}'
    ENDELSE

  END

  FREE_LUN, filenum

end