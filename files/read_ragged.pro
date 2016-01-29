function read_ragged, filename

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
IF(abs(io_in - io_check) GT 0.001) THEN BEGIN
  close, 1
  print, "File read error!"
  ;return
  stop
ENDIF

readu, 1, commit_in

n_dims = int_type
readu, 1, n_dims
print, n_dims

if(n_dims LT 0) THEN BEGIN
  n_dims = abs(n_dims)
  ;now we make a list of each row eh?
  dims = lonarr(n_dims)
  readu, 1, dims
  lengths = lonarr(dims[1])
  readu, 1, lengths
  IF my_type_code EQ 'f' THEN BEGIN
    data_list = LIST(length = dims[1])
    axes_list = LIST(length = dims[1])
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = fltarr(lengths[i])
      readu, 1, tmp
;      print, minmax(tmp)
      data_list[i] = tmp
    END
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = fltarr(lengths[i])
      readu, 1, tmp
;      print, minmax(tmp)
      axes_list[i] = tmp
    END

  ENDIF ELSE BEGIN
    data_list = LIST(length = dims[1])
    axes_list = LIST(length = dims[1])
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = dblarr(lengths[i])
      readu, 1, tmp
      data_list[i] = tmp
    END
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = dblarr(lengths[i])
      readu, 1, tmp
      axes_list[i] = tmp
    END
  ENDELSE

  data = {id:id_in,data:data_list, axes:axes_list}

ENDIF ELSE BEGIN
  print, "Array is normal. Use read_data.pro"

ENDELSE

close, 1

return, data

end