function get_p_res, p_perp_ax, omega, om_ce=om_ce, om_pe=om_pe
;Calculate resonant p_parallel from cold plasma dispersion for omega

if(N_ELEMENTS(om_ce) EQ 0) THEN om_ce=30000
if(N_ELEMENTS(om_pe) EQ 0) THEN om_pe=om_ce*3
 v0 =  299792458.000

sz=size(p_perp_ax)
ret_val=p_perp_ax

k = get_dispersion(omega, om_ce=om_ce, om_pe=om_pe, /k)
ck_om_sq = (k*v0/omega)^2
gamma_denom = (ck_om_sq -1.0)*(omega/om_ce)
gamma_num_part = gamma_denom*(omega/om_ce)

FOR i=0, sz[1]-1 DO BEGIN
  tmp = gamma_num_part*(1.0+ (p_perp_ax[i]/v0)^2) + 1.0
  tmp=sqrt(tmp)
  gamma_res = (-1.0 + sqrt(ck_om_sq)*tmp)/gamma_denom
  ret_val[i] = (gamma_res*omega - om_ce)/k
print, gamma_res
END

return, ret_val
END

