function read_data, filename

;Written for commit ID from dc9e387 to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

OPENR, filenum,  filename, /GET_LUN
;open file

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

data = read_block(filenum, hdr.my_type, hdr.block_type)

readu, filenum, next_block

id_type ='1234567891'
id_in = id_type
;type for block_izes and block id...

readu, filenum, id_in
data=create_struct(data, {block:id_in})

FREE_LUN, filenum

return, data

end
