function read_spect, filename

;Written for commit ID from 3b24e7d to  FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

;my_type_code = 'f'
;my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

;int_type = 1
OPENR, filenum,  filename, /GET_LUN
;open file

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

hdr_end = hdr.block_type
POINT_LUN, -filenum, hdr_end

POINT_LUN, filenum, 0
int_sz=1
readu, filenum, int_sz

tmp = FSTAT(filenum)
POINT_LUN, filenum, (tmp.size - int_sz)
start_pos = hdr.block_type
readu, filenum, start_pos

POINT_LUN, filenum, hdr_end
;Grab int size and footer start

spect = {B:read_block(filenum, hdr.my_type, hdr.block_type)}
;Read one block

next_pos = hdr.block_type
readu, filenum, next_pos

PRINT, next_pos, start_pos, fstat(filenum)
if(next_pos EQ start_pos) THEN BEGIN
  print, "Insufficient arrays found, spectrum incomplete"
  return, spect
end

POINT_LUN, filenum, next_pos
;Skip on to the next block..

;Now the next complete array
hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

spect = create_struct({ang:read_block(filenum, hdr.my_type, hdr.block_type)}, spect)

readu, filenum, next_pos
if(next_pos NE start_pos) THEN BEGIN
  print, "Extra arrays in input file"
  return, spect
end

POINT_LUN, filenum, start_pos
readu, filenum, next_pos

id_type ='1234567891'
id_in = id_type
;type for block_izes and block id...

readu, filenum, id_in
PRINT, id_in
spect=create_struct(spect, {block:id_in})


return, spect
END
