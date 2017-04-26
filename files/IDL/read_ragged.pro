function read_ragged, filename

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

my_type_code = 'f'
my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

int_type = 1

openr, filenum, filename, /GET_LUN
;open file

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

n_dims = int_type
readu, filenum, n_dims
print, n_dims

if(n_dims LT 0) THEN BEGIN
  n_dims = abs(n_dims)
  ;now we make a list of each row eh?
  dims = lonarr(n_dims)
  readu, filenum, dims
  lengths = lonarr(dims[1])
  readu, filenum, lengths
  IF my_type_code EQ 'f' THEN BEGIN
    data_list = LIST(length = dims[1])
    axes_list = LIST(length = dims[1])
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = fltarr(lengths[i])
      readu, filenum, tmp
;      print, minmax(tmp)
      data_list[i] = tmp
    END
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = fltarr(lengths[i])
      readu, filenum, tmp
;      print, minmax(tmp)
      axes_list[i] = tmp
    END

  ENDIF ELSE BEGIN
    data_list = LIST(length = dims[1])
    axes_list = LIST(length = dims[1])
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = dblarr(lengths[i])
      readu, filenum, tmp
      data_list[i] = tmp
    END
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = dblarr(lengths[i])
      readu, filenum, tmp
      axes_list[i] = tmp
    END
  ENDELSE

  data = {id:id_in,data:data_list, axes:axes_list}

ENDIF ELSE BEGIN
  print, "Array is normal. Use read_data.pro"

ENDELSE

FREE_LUN, filenum

return, data

end
