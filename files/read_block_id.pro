function read_block_id, filename
;Read the block ID (first element)
;Written for commit ID from dc9e387 to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, ''

OPENR, filenum,  filename, /GET_LUN
;open file

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, ''
ENDIF

POINT_LUN, filenum, 0
int_sz=1
readu, filenum, int_sz

tmp = FSTAT(filenum)
POINT_LUN, filenum, (tmp.size - int_sz)
start_pos = hdr.block_type
readu, filenum, start_pos
POINT_LUN, filenum, start_pos
readu, filenum, start_pos
id_type ='1234567891'
id_in = id_type
;type for block id...

readu, filenum, id_in

FREE_LUN, filenum
return, id_in
END
