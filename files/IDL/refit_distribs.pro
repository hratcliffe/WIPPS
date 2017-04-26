
pro smooth_distrib, filename, outfile=outfile, smth=smth
  ;Read a distribution from filename, apply smoothing and dump to outfile if given, else to filename with _smth inserted before extension
  ;NB we renormalise things so the axes and data are as expected by lookup reader in calculate_growth
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  IF(N_ELEMENTS(smth) EQ 0) THEN smth = 1
  IF(N_ELEMENTS(norm_extr) EQ 0) THEN norm_extr = 1
  ;Insert smth into string or onto end if no extension
  if(N_ELEMENTS(outfile) EQ 0) THEN BEGIN
    IF((strpos(filename, '.') NE -1)) THEN BEGIN outfile = strjoin((strsplit(filename, '.', /extract))[0:-2], '.')+'_smth.'+(strsplit(filename, '.', /extract))[-1]
    ENDIF ELSE outfile=filename+'_smth'
  ENDIF
  dir_in = (stregex(filename, '(.*/+)[^\/]+$', /extract, /subexpr))[1]
  print, 'Dir is ', dir_in
  dist_in = read_distribs(filename)
  deck_specs=read_deck_all(dir=dir_in)
  
  dist_in.data = gauss_smooth(dist_in.data, smth)
  norm_factor = double(1.0)/(deck_specs.x_max*deck_specs.y_max*deck_specs.dens)/((dist_in.axes.x[1]-dist_in.axes.x[0])/m0*((dist_in.axes.y[1]-dist_in.axes.y[0])/m0)^2)
  d_x_sz = (size(dist_in.data))[1]
  d_y_sz = (size(dist_in.data))[2]
  z_renorm = 1.0/dist_in.data[d_x_sz/2, d_y_sz/2]*total(dist_in.data[d_x_sz/2, *])
  norm_factor = float(norm_factor/z_renorm)
  print, "Norm is ", norm_factor, "and z renorm is", z_renorm
  err=write_data(outfile, dist_in.data*norm_factor, list(dist_in.axes.x/float(m0), dist_in.axes.y/float(m0)), time=dist_in.time, space=dist_in.space, B_ref=dist_in.B_ref, block=dist_in.block)
  PRINT, err
  if(err) THEN PRINT, "Error writing file"

end
