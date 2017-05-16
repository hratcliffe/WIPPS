function read_distribs, filename
;Read an unknown number of identical structs into an array

;Written for commit ID from dc9e387 to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

id_type ='1234567891'
;type for block_izes and block id...

IF(~FILE_TEST(filename)) THEN BEGIN
  PRINT, 'File '+filename+' not found'
  RETURN, !NULL
END
OPENR, filenum,  filename, /GET_LUN
;open file

;First we grab the position of end block for termination condition
;Read header to get block_type
hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

POINT_LUN, filenum, 0
int_sz=1
readu, filenum, int_sz

;Read last element for final footer start
tmp = FSTAT(filenum)
POINT_LUN, filenum, (tmp.size - int_sz)
end_pos = hdr.block_type
readu, filenum, end_pos
;This is now where the final closer should start

;Rewind file
POINT_LUN, filenum, 0
next_block = hdr.block_type
POINT_LUN, -filenum, next_block

i=0

WHILE(~EOF(filenum) AND (next_block NE end_pos) ) DO BEGIN

  hdr=read_header(filenum)
  IF(hdr.err) THEN BEGIN
    FREE_LUN, filenum
    PRINT, "Error reading file header"
    RETURN, !NULL
  ENDIF

  data = read_block(filenum, hdr.my_type, hdr.block_type)
  IF(N_ELEMENTS(data) EQ 0) THEN continue
  next_block = hdr.block_type
  readu, filenum, next_block

  ;Read in the time and space fields and the B_ref
  space_in = [hdr.block_type, hdr.block_type]
  time_in = [hdr.my_type, hdr.my_type]
  B_ref = hdr.my_type

  readu, filenum, time_in
  readu, filenum, space_in
  readu, filenum, B_ref
  data=create_struct(data, {time: time_in, space:space_in, B_ref: B_ref})
  readu, filenum, next_block

  id_in = id_type

  readu, filenum, id_in
  data=create_struct(data, {block:id_in})
  data=create_struct(data, {filename: filename})

  POINT_LUN, filenum, next_block

  IF(i EQ 0) THEN BEGIN
    n_dats = FIX(end_pos/n_tags(data, /length))
    all_data = replicate(data, n_dats)

  ENDIF ELSE BEGIN
    all_data[i] = data
  ENDELSE
  i=i+1
END

;Skip values in footer, just read the last bit
readu, filenum, next_block
;If it seems to be a valid position, use it
IF next_block GT end_pos AND next_block LT tmp.size THEN BEGIN
  POINT_LUN, filenum, next_block
END
readu, filenum, next_block
POINT_LUN, -filenum, next_block
;Final ID may or may not be there
IF(next_block LT tmp.size) THEN BEGIN
  id_in = id_type
  readu, filenum, id_in
ENDIF ELSE BEGIN
  id_in =""
  IF((size(all_data))[1] GT 0) THEN id_in = all_data[0].block
ENDELSE
print, "Read " + string(format='(I3)', i, /print)+ " blocks of "+id_in
FREE_LUN, filenum

return, all_data

end
