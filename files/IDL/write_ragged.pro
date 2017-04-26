function write_ragged, filename, data, axes, id=id

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....
;The id string is done as exactly 10 chars ending in null. To faffy to put the null terminator within the string
COMPILE_OPT IDL2
;force long ints and proper brackets

sz = size(data)

IF((sz[-2] NE 11 ) || ~ISA(data, "list")) THEN BEGIN
  PRINT, "Data is not a list, try write_data"
  RETURN, 1
END
;Check we have list type. The two checks are redundant, but hey, short circuiting :)

n_dims = sz[1]
;Not strictly number of dims but of rows, but is what we need
dims=lonarr(n_dims)
FOR i=0, sz[1]-1 DO dims[i] = (size(data[i]))[1]
;Note sz is already LONG ints

IF(n_dims EQ 0) THEN RETURN, 1
;And a sane dimensionality

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

writeu, 1, -n_dims
writeu, 1, dims

FOR i=0, n_dims-1 DO BEGIN
  writeu, 1, data[i]
END

FOR i=0, n_dims-1 DO BEGIN
  writeu, 1, axes[i]
END

close, 1

RETURN, 0

end


