function read_block, filenum, my_type_code, id_in

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

;Check file is good
fileinf=FSTAT(filenum)
IF(NOT fileinf.open) THEN BEGIN
  PRINT, "Invalid file handle"
  RETURN, 1
ENDIF

n_dims = 1
readu, filenum, n_dims

if(n_dims GT 0) THEN BEGIN
  dims = lonarr(n_dims)
  readu, filenum, dims
  IF my_type_code EQ 'f' THEN BEGIN
    axes_list = {x:fltarr(dims[0])}
    if(n_dims GT 1) THEN axes_list = create_struct(axes_list, {Y:fltarr(dims[1])})
    if(n_dims GT 2) THEN axes_list = create_struct(axes_list, {Z:fltarr(dims[2])})

    data = {id:id_in,data:fltarr(dims), axes:axes_list}
    tmp2=fltarr(dims)

  ENDIF ELSE BEGIN
    axes_list = {k:dblarr(dims[0]), omega:dblarr(dims[1])}
    axes_list = {x:dblarr(dims[0])}
    if(n_dims GT 1) THEN axes_list = create_struct(axes_list, {Y:dblarr(dims[1])})
    if(n_dims GT 2) THEN axes_list = create_struct(axes_list, {Z:dblarr(dims[2])})

    data = {id:id_in,data:dblarr(dims), axes:axes_list}
    tmp2=dblarr(dims)

  ENDELSE

  ;Can't read directly into anon structure field, so use tmp
  readu, filenum, tmp2
  data.data = tmp2
  tmp2=0
ENDIF ELSE BEGIN
  print, "Array is ragged. Use read_ragged.pro"
  RETURN, !NULL
ENDELSE


IF my_type_code EQ 'f' THEN BEGIN
  tmpa=fltarr(dims[0])
  readu, filenum, tmpa
  data.axes.X = tmpa

  IF(n_dims GT 1) THEN BEGIN
    tmpa=fltarr(dims[1])
    readu, filenum, tmpa
    data.axes.Y = tmpa
  ENDIF
  IF(n_dims GT 2) THEN BEGIN
    tmpa=fltarr(dims[2])
    readu, filenum, tmpa
    data.axes.Z = tmpa
  ENDIF
ENDIF ELSE BEGIN
  tmpa=dblarr(dims[0])
  readu, filenum, tmpa
  data.axes.X = tmpa

  IF(n_dims GT 1) THEN BEGIN
    tmpa=dblarr(dims[1])
    readu, filenum, tmpa
    data.axes.Y = tmpa
  ENDIF
  IF(n_dims GT 2) THEN BEGIN
    tmpa=dblarr(dims[2])
    readu, filenum, tmpa
    data.axes.Z = tmpa
  ENDIF
ENDELSE

tmpa = 0

RETURN, data
END