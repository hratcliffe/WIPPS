function check_energy, dump_time, E=E, H=H, str=str, dir=dir
;Check individual energy totals and summed conservation
;One of E, H and str should be given, last will be used, and filenames are assumed to have form 'x_p[xy]_[str]_dump_time.dat where str is EnergeticE (E=1), EnergeticEH (H=1) or the given string. If none are given E is used

filestr = 'EnergeticE'
IF(KEYWORD_SET(E)) THEN filestr = 'EnergeticE'
IF(KEYWORD_SET(H)) THEN filestr = 'EnergeticEH'
IF(N_ELEMENTS(str) GT 0) THEN filestr = str

dump_str = strtrim(string(dump_time, /print), 2)

m0 = 9.10938291d-31
v0 = 2.99792458d8

tmp_str = 'x_px_'+filestr+'_'+dump_str+'.dat'
print, tmp_str
x_dist = read_distribs(dir+tmp_str)
tmp_str = 'x_py_'+filestr+'_'+dump_str+'.dat'
print, tmp_str
y_dist = read_distribs(dir+tmp_str)

gamm_x = sqrt(1.0 + 1.0*(x_dist.axes.x/m0/v0)^2)

dpx = abs(x_dist.axes.x[0]-x_dist.axes.x[1])/m0
E_x = total(x_dist.data*m0*v0*v0*(gamm_x-1.0)) ;*dpx

gamm_y = sqrt(1.0 + 1.0*(y_dist.axes.x/m0/v0)^2)

dpy = abs(y_dist.axes.x[0]-y_dist.axes.x[1])/m0
E_y = total(y_dist.data*m0*v0*v0*(gamm_y-1.0)) ;*dpy

ret = create_struct('E_x', E_x, 'E_y',E_y, 'E_tot',E_x+E_y)
return, ret
END

