
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
;IF((sz[-2] EQ 5) my_type_code ='d'
sz_sz = 8
;Intsz will be 8 bytes at the moment, matching size_t on my system
my_sz = 4
IF((sz[-2] EQ 5)) THEN my_sz = 8
;Floats are 4 or 8 bytes
;this matches the type of the C code my_type...

id_out = '          '
IF(N_ELEMENTS(id) GT 0) THEN strput, id_out, id, 0
;copy id into correct length string
;type for block id...

io_check = 3.0/32.0
if(my_sz EQ 8) THEN io_check = 3.0d0/32.0d0
;io verification const.

commit_out ='IDL data write'
;type for commit id NB 14 chars long!!!

null_byte = BYTE(0)

openw, 1, filename
;open file
writeu, 1, sz_sz
writeu, 1, my_sz

writeu, 1, io_check
writeu, 1, commit_out
writeu, 1, null_byte

next_loc = 0ull
next_loc = next_loc + 2*sz_sz+my_sz+15+sz_sz*(n_dims+2)
writeu, 1, next_loc
print, size(next_loc), next_loc

writeu, 1, ulong64(n_dims)
writeu, 1, ulong64(dims)

next_loc = next_loc + ulong64(product(dims))*my_sz+sz_sz
writeu, 1, next_loc
print, size(next_loc), next_loc

writeu, 1, data

next_loc = next_loc+sz_sz+ulong64(total(dims))*my_sz
writeu, 1, next_loc
FOR i=0, n_dims-1 DO BEGIN
  writeu, 1, axes[i]
END
ftr_start = next_loc
;Footer data
next_loc  = next_loc + my_sz*3 + sz_sz*3
writeu, 1, next_loc

next_loc = next_loc + sz_sz + 10
writeu, 1, next_loc
block_out = '10 chars '
writeu, 1, block_out
writeu, 1, null_byte
writeu, 1, ftr_start

close, 1

RETURN, 0

end
