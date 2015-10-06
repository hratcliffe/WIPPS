function read_data, filename

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

my_type_code = 'f'
my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

id_type ='1234567891'
id_in = id_type
;type for block id...

io_type =0.0
io_check = 3.0/32.0
io_in=io_type
;type for io verification const.

commit_type ='123456789112345'
commit_in = commit_type
;type for commit id

int_type = 1

openr, 1, filename
;"Tmp.txt"
;open file

readu, 1, id_in
readu, 1, io_in
readu, 1, commit_in

n_dims = int_type
readu, 1, n_dims

if(n_dims GT 0) THEN BEGIN
  dims = lonarr(n_dims)
  readu, 1, dims
  IF my_type_code EQ 'f' THEN BEGIN
    axes_list = {k:fltarr(dims[0]), omega:fltarr(dims[1])}
    data = {id:id_in,data:fltarr(dims[0], dims[1]), axes:axes_list}
  ENDIF ELSE BEGIN
    axes_list = {k:dblarr(dims[0]), omega:dblarr(dims[1])}
    data = {id:id_in,data:dblarr(dims[0], dims[1]), axes:axes_list}
  ENDELSE

  ;Now we do need the right majority...

  tmp2=fltarr(dims[0], dims[1])
  ;Can't read directly into anon structure field, so use tmp
  readu, 1, tmp2
  data.data = tmp2
  tmp2=0
ENDIF ELSE BEGIN
  print, "Array is ragged. Use read_ragged.pro"

ENDELSE

;EXPLICTLY 2-D here
tmpa=fltarr(dims[0])
;Can't read directly into anon structure field, so use tmp
readu, 1, tmpa
data.axes.k = tmpa

tmpa=fltarr(dims[1])
readu, 1, tmpa
data.axes.omega = tmpa

tmpa = 0

close, 1

return, data

end
