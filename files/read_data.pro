function read_data, filename

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

my_type_code = 'f'
my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

int_type = 1
OPENR, filenum,  filename, /GET_LUN
;open file

error=read_header(filenum)
IF(error) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

n_dims = int_type
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
  FREE_LUN, filenum
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
FREE_LUN, filenum


return, data

end
