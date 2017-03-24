function write_header, filenum, my_sz_type

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....
COMPILE_OPT IDL2
;force long ints and proper brackets

sz_sz = 8
;Intsz will be 8 bytes at the moment, matching size_t on my system
my_sz = 4
IF((my_sz_type EQ 5)) THEN my_sz = 8
;Floats are 4 or 8 bytes
;this matches the type of the C code my_type...

io_check = 3.0/32.0
if(my_sz EQ 8) THEN io_check = 3.0d0/32.0d0
;io verification const.

commit_out ='IDL data write'
;type for commit id NB 14 chars long!!!

null_byte = BYTE(0)

writeu, filenum, sz_sz
writeu, filenum, my_sz

writeu, filenum, io_check
writeu, filenum, commit_out
writeu, filenum, null_byte

IF(my_sz EQ 4) THEN BEGIN
  my_type = 0.0
ENDIF ELSE IF(my_sz EQ 8) THEN BEGIN
  my_type = 0.0d
END

IF(sz_sz EQ 2) THEN BEGIN
  block_type = 0U
ENDIF ELSE IF(sz_sz EQ 4) THEN BEGIN
  block_type = 0UL
ENDIF ELSE IF(sz_sz EQ 8) THEN BEGIN
  block_type = 0ULL
END

hdr_info = create_struct({err:0, int_sz:sz_sz, type_sz: my_sz, my_type:my_type,  block_type:block_type})

return, hdr_info

end
