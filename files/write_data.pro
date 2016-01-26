
function write_data, filename, data, axes, id=id

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....
;The id string is done as exactly 10 chars ending in null. To faffy to put the null terminator within the string
COMPILE_OPT IDL2
;force long ints and proper brackets

sz = size(data)
n_dims = sz[0]
dims= sz[1:n_dims]
;Note sz is already LONG ints

IF((sz[-2] NE 4 && sz[-2] NE 5) || ~ISA(data, /array)) THEN BEGIN
  PRINT, "Data is not float or double array, try write_ragged"
  RETURN, 1
END
;Check we have float or double not List

IF(n_dims EQ 0) THEN RETURN, 1
;And a sane dimensionality

my_type_code = 'f'
my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

id_out = '         '
IF(N_ELEMENTS(id) GT 0) THEN strput, id_out, id, 0
;copy id into correct length string
;type for block id...

io_check = 3.0/32.0
;io verification const.

commit_out ='IDL data write'
;type for commit id

null_byte = BYTE(0)

openw, 1, filename
;open file
writeu, 1, id_out
writeu, 1, null_byte
writeu, 1, io_check
writeu, 1, commit_out
writeu, 1, null_byte

writeu, 1, n_dims
writeu, 1, dims

writeu, 1, data

FOR i=0, n_dims-1 DO BEGIN
  writeu, 1, axes[i]
END

close, 1

RETURN, 0

end
