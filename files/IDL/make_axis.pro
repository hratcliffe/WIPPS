function make_axis, n_pts, res=res, len=len

; Created by H Ratcliffe, University of Warwick
;16/4/2015
;Defines axes for field FFTs
;-------Params------
; n_pts: number of points
; Supply res OR len::
; res: resolution in normalised units
; len: total length in normalised units
; Returned axis is normalised to the inverse of the supplied res
;Call example:
; ax=make_axis(512, res=1)

IF((N_ELEMENTS(res) GT 0) AND (N_ELEMENTS(len) GT 0)) THEN BEGIN
  IF((res * n_pts)-len GT res) THEN print, 'WARNING: Supplied res and len inconsistent. Using res value'
ENDIF
; Check consistency of res and len if both supplied

IF((N_ELEMENTS(res) EQ 0) AND (N_ELEMENTS(len) GT 0)) THEN res = len/n_pts
;Calc res if only len supplied

IF((N_ELEMENTS(res) EQ 0) AND (N_ELEMENTS(len) EQ 0)) THEN BEGIN
  print, 'No length or resolution supplied'
  return, -1
ENDIF

n_x2=float(n_pts)/2.
ax=!pi*(findgen(n_pts)-n_x2)/n_x2/res
return, ax

;print, 'Make axis failed. Sorry about that.'


END
