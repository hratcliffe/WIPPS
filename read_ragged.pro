;pro read_ragged
;Written for commit ID from b8b80f3 to ... FILL IN IF IO CHNAGES....

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

openr, 1, filename
;"Tmp.txt"
;open file

readu, 1, id_in
readu, 1, io_in
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

    FOR i=0, dims[1]-1 DO BEGIN
      tmp = fltarr(lengths[i])
      readu, 1, tmp
      print, minmax(tmp)
      data_list[i] = tmp
    END
    axes_list = {k:fltarr(1), omega:fltarr(dims[1])}
  ENDIF ELSE BEGIN
    data_list = LIST(length = dims[1])
    FOR i=0, dims[1]-1 DO BEGIN
      tmp = dblarr(lengths[i])
      readu, 1, tmp
      data_list[i] = tmp
    END

    axes_list = {k:dblarr(dims[0]), omega:dblarr(dims[1])}
  ENDELSE

  data = {id:id_in,data:data_list, axes:axes_list}


  ;Now we do need the right majority...

ENDIF ELSE BEGIN
  print, "Array is normal. Use read_data.pro"

ENDELSE

close, 1

end
