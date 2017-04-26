COMPILE_OPT IDL2

function is_a_distrib, candidate, fail=fail
  ;Check given candidate is valid dsitribs structure. I.e. contains needed tags, has valid data type, and axis sizes match data. Returns 1 if so, else 0 and fail is set to 1, 2, 3 for respective fails
  ;In particular, this checks everything we assume in write_distribs
  ;Check has the needed tags to be a distrib
  IF(N_ELEMENTS(fail) EQ 0) THEN fail = 0
  fail = 1
  tag_expected=['data', 'axes', 'time', 'space', 'b_ref', 'block', 'filename']
  n_tags =(size(tag_expected))[1]
  ;n_tags_cand =(size(tag_names(candidate)))[1]
  indices = where(strupcase(tag_expected) EQ tag_names(candidate), nmatch)
  ;IF(n_tags_cand NE n_tags || nmatch NE n_tags) THEN RETURN, 0
  ;This would exclude excess tags, but we probably don't want to
  IF(nmatch NE n_tags) THEN RETURN, 0
  ;Missing tags, so no good
  fail = 2
  data_type = size(candidate.data, /type) 
  IF(data_type NE 4 && data_type NE 5) THEN RETURN, 0
  IF(size(candidate.axes, /type) NE 8) THEN RETURN, 0
  ;check axes is a struct
  n_dims = (size(candidate.data))[0]
  ax_tags = tag_names(candidate.axes)
  IF(TOTAL(ax_tags EQ (['X', 'Y', 'Z'])[0:n_dims-1]) NE n_dims) THEN RETURN, 0
  ;Wrong number of axes, or wrong subset.
  FOR i=0, n_dims-1 DO IF(total(size(candidate.axes.(i), /type) EQ [4, 5]) NE 1) THEN RETURN, 0
    
  fail = 3
  dat_sz = (size(candidate.data))
  ax_sz = (size(candidate.axes.x))[1]
  IF(dat_sz[0] GE 2) THEN ax_sz = [ax_sz, (size(candidate.axes.y))[1]]
  IF(dat_sz[0] GE 3) THEN ax_sz = [ax_sz, (size(candidate.axes.z))[1]]
  tmp=where(dat_sz[1:-3] EQ ax_sz, nmat)
  IF(nmat NE dat_sz[0]) THEN RETURN, 0
  
  fail = 0
  ;If we get here, we're good
  RETURN, 1
END

function write_distribs, filename, distribs
;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....
;Write an unknown number of distribs to a file. 
COMPILE_OPT IDL2

;Check input
sz = size(distribs)
n_distribs = sz[1]
IF(sz[0] NE 1) THEN BEGIN
  PRINT, "Please supply either a single or 1-D array of distributions"
  RETURN, 1
ENDIF
;Since we're taking an array all the elements must be the same struct...

IF(~is_a_distrib(distribs[0], fail=fail)) THEN BEGIN
  PRINT, "Invalid distribution (err "+string(format='(I1)', fail, /print)+')'
  RETURN,  1
ENDIF

openw, filenum, filename, /get_lun

for i=0, n_distribs-1 DO BEGIN
  axes_list = list(distribs[i].axes.x, distribs[i].axes.y)
  IF(i EQ n_distribs-1) THEN close_file = 1 ELSE close_file = 0
  err = write_data(filenum, distribs[i].data, axes_list, block=distribs[i].block, space=distribs[i].space, time=distribs[i].time, b_ref=distribs[i].b_ref, /usenum, close_file=close_file)
end

free_lun, filenum

RETURN, 0

end
