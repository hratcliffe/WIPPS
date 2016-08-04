function read_data, filename, is_d=is_d

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

tmp = FSTAT(filenum)
POINT_LUN, filenum, (tmp.size - hdr.int_sz)
start_pos = hdr.block_type
readu, filenum, start_pos
POINT_LUN, filenum, start_pos

next_block = start_pos
readu, filenum, next_block

id_type ='1234567891'
id_in = id_type
;type for block_izes and block id...

readu, filenum, id_in
data=create_struct(data, {block:id_in})

if(KEYWORD_SET(is_d)) THEN BEGIN
  wave_in = 1l
  readu, filenum, wave_in

  tag_in = id_type
  readu, filenum, tag_in
  data=create_struct(data, {wave_id: wave_in, tag:tag_in})

END

FREE_LUN, filenum

return, data

end
