;pro read_data

COMPILE_OPT IDL2
;force long ints and proper brackets

my_type_code = 'f'
my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

id_type ='12345678910'
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

openr, 1, "Tmp.txt"
;open file

readu, 1, id_type
readu, 1, io_in
readu, 1, commit_type

n_dims = int_type
readu, 1, n_dims

dims = lonarr(n_dims)

readu, 1, dims

IF my_type_code EQ 'f' THEN BEGIN
  axes_list = {k:fltarr(dims[0]), omega:fltarr(dims[1])}
  data = {id:id_in,data:fltarr(dims[1], dims[0]), axes:axes_list}
ENDIF ELSE BEGIN
  axes_list = {k:dblarr(dims[0]), omega:dblarr(dims[1])}
  data = {id:id_in,data:dblarr(dims[1], dims[0]), axes:axes_list}
ENDELSE
;Now we do need the right majority...

tmp2=fltarr(dims[1], dims[0])
;Can't read directly into anon structure field, so use tmp
readu, 1, tmp2
data.data = tmp2
tmp2=0

;FOR i=0, dims[1] -1 DO BEGIN
;  FOR j=0, dims[0] -1 DO BEGIN
;    readu, 1, tmp
;    data.data[i,j] = tmp
;  END
;END
;file.read((char *) data_tmp , sizeof(my_type)*dims[0]*dims[1]);

;readu, 1, data.data

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

end
