function get_dispersion, input, om_ce = om_ce, om_pe = om_pe, theta=theta, k=k, appr=appr
;Get the corresponding value for input. If k is set then return k for omega=input, else, return omega = k=input. Reference om_ce and om_pe are required, else 10000 and 30000 are assumed. Theta in radians, appr flag sets to use simpler (fully reversible) form

  common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
  if(n_Elements(om_ce) EQ 0) THEN om_ce = 10000*2.0*!pi
  if(n_Elements(om_pe) EQ 0) THEN om_pe = 30000*2.0*!pi
  if(n_Elements(theta) EQ 0) THEN theta = 0

  if(keyword_set(k)) THEN BEGIN
    IF(KEYWORD_SET(appr)) THEN BEGIN tmp = -((om_pe*om_pe/(input*(input - om_ce * cos(theta)))))
    ENDIF ELSE BEGIN tmp = (1.0 - (om_pe*om_pe/(input*(input - om_ce * cos(theta)))))
    ENDELSE
    k_out = input/v0 *sqrt(tmp)
    return, k_out
  ENDIF ELSE BEGIN
    kcsq = (v0*input)^2
    om_out =  kcsq*om_ce*cos(theta)/(kcsq + om_pe^2)
    return, om_out

  ENDELSE

end

