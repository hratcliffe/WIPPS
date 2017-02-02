function pancake_curves, om_ratio, x_c, x_ax=x_ax
  ;Calculate the pancake distribution as described in Meredith 1999 or Summers, 1998
  ;om_ratio is Om_ce/om_pe, x_c is desired curve and x_ax is the (normalised to Om_ce) resonant frequency required

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck

  IF(N_ELEMENTS(x_ax) LT 1) THEN x_ax= findgen(100)/100

  E_crit = m0*v0*v0/2.0*(om_ratio)^2

  C = sqrt(2.0*E_crit/m0)
  len = (size(x_ax))[1]
  v_perp = fltarr(len)
  v_par=fltarr(len)
  
  FOR i=0, len-1 DO BEGIN
    v_perp[i] = C*sqrt(pancake_f(x_c)-pancake_f(x_ax[i]))  
    v_par[i] = C* pancake_g(x_ax[i])
  END
  return, create_struct({v_perp:v_perp, v_par:v_par, x_ax:x_ax})
end

function pancake_f, x
;Positive and Monotonic decreasing from x=0 to 1

  return, ( (1.0-x)/x + 2.0*x - alog(x))
end
function pancake_g, x
;Negative and Monotonic increasing from x=0 to 1
  return, -sqrt((1.0-x)^3/x)
end

function pancake_f_wrapper, x, val=val
;Wrap pancake_f to equal val. Val MUST be present
  return, pancake_f(x) - val
end

function pancake_g_wrapper, x, val=val
;Wrap pancake_g to equal val. Val MUST be present

;  if(x LT 0 OR x GT 1) THEN return, [-1e30, -1e30, -1e30]
  return, pancake_g(x) - val
  ;avoid complex returns for x> 1
end

;Compile overriding copy of fx_root
;@ fx_root_splunge.pro
function binary_invert, val, func_name, precision=precision, _extra=ext
;Locate val in func_name. Assumes monotonic and non-constant on args from 0 to 1

;Check increasing or decreasing
  grad = 1
  IF( call_function(func_name,0.2, _extra=ext) LT call_function(func_name,0.1, _extra=ext)) THEN grad = -1

  ;Fractional precision to terminate
  IF(N_ELEMENTS(precision) EQ 0) THEN precision = 0.00001d0
  selection = 0.5d0
  tmp = call_function(func_name,selection, _extra=ext)
  current_increment = 0.25d0
  counter = 0
  max_it = 50 ; 0.5^50 = 8e-16
  while(abs((tmp - val)) GT precision AND counter LT max_it) DO BEGIN
    tmp = call_function(func_name,selection, _extra=ext)
    if(tmp GT val) THEN selection = selection - grad*current_increment
    if(tmp LT val) THEN selection = selection + grad*current_increment
    if(tmp EQ val) THEN break
    current_increment = current_increment * 0.5
    counter = counter + 1
  END

  IF(counter GT max_it -1) THEN return, -1
  
  return, selection

end
function pancake_distribution, om_ratio, px_ax, py_ax, f_py_0
  ;Create a pancake distribution. We know analytically the isolines for given resonant velocity. So we require supply of the f_py_0 values of density for each supplies ps in the x_axis supplied
  ;Then we have to invert the iso lines to work out which x_c the required point lies on
  ;and then we interpolate f_py_0 for this pt
  ;I don't think the functions are both analytically invertible
  ;So is easiest to create us a grid. We know g is negative and increasing, and f positive and increasing. And we know the rings move outwards for decreasing x_c. We can ignore quadrants and consider only e.g. v_par, v_perp > 0
  ;Then as v_par is independent of x_c we can "invert" g by binary bisect hunting to find x
  ;Then we do the same using f to get x_c
  
  ;Use fx_root (secant type method) The docs are awful so remember, the 3 element initial guess is the best guess and the bounds on the root. These mustn't be equal, and the second two must be ordered and it seems like it'll find one in this band first???? Seriously the docs suck. I have a local copy which allows passing the required val down to the function
  
  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  x_sz =(size(px_ax))[1]
  y_sz =(size(py_ax))[1]
  IF(N_ELEMENTS(f_py_0) NE x_sz) THEN BEGIN
    PRINT, 'Size if f_py_0 does not match x axis'
    RETURN, !NULL
  ENDIF

  E_crit = m0*v0*v0/2.0*(om_ratio)^2
  C = sqrt(2.0*E_crit/m0)

  ;Correct axes back to v rather than p
  ;TODO add relativistic correction
  vx_ax = px_ax/m0
  vy_ax = py_ax/m0
  distrib = fltarr(x_sz, y_sz)
  FOR i=0, x_sz-1 DO BEGIN
    IF(I MOD 10 EQ 0) THEN print, i
    FOR j=0, y_sz-1 DO BEGIN
      ;Do the inversion for this pt
      ;same guess every time
      x_current = binary_invert(0.0, 'pancake_g_wrapper', val = -abs(vx_ax[i])/C)
;       v_perp[i] = C*sqrt(pancake_f(x_c)-pancake_f(x_ax[i]))  
; --> (v_perp/C)^2 + pancake_f(x_current) = pancake_f(x_c)
      val = (vy_ax[j]/C)^2 + pancake_f(x_current)
      x_crit = binary_invert(0.0, 'pancake_f_wrapper', val=val)
      ;Now we know the x_crit value, so we can get the x_axis density for this pt
      v_par_0 = abs(C*pancake_g(x_crit))
      lower_cell_bnd = value_locate(vx_ax, v_par_0)
      IF(Lower_cell_bnd LT 0) THEN dens = f_py_0[0]
      IF(lower_cell_bnd GE 0 AND lower_cell_bnd LT x_sz-1) THEN dens = (f_py_0(lower_cell_bnd)*abs(v_par_0-vx_ax(lower_cell_bnd)) + f_py_0(lower_cell_bnd+1)*abs(vx_ax(lower_cell_bnd +1)-v_par_0))/abs(vx_ax(lower_cell_bnd +1)-vx_ax(lower_cell_bnd))
      IF(lower_cell_bnd EQ x_sz -1) THEN dens = f_py_0(lower_cell_bnd)
      distrib(i, j) = dens
    END
  END
  return, distrib
end


