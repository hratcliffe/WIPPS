function read_data, filename, is_d=is_d, noext=noext
;The ext keyword is a temporary as some files have additional info in (commit after 3e3b01c)

;Written for commit ID from dc9e387 to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

IF(~FILE_TEST(filename)) THEN BEGIN
  PRINT, "File " + filename+" not found"
  RETURN, !NULL
END
OPENR, filenum,  filename, /GET_LUN
;open file

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

data = read_block(filenum, hdr.my_type, hdr.block_type)

next_block = hdr.block_type
readu, filenum, next_block

IF(~KEYWORD_SET(NOEXT)) THEN BEGIN
  ;Read in the time and space fields and the B_ref
  space_in = [hdr.block_type, hdr.block_type]
  time_in = [hdr.my_type, hdr.my_type]
  B_ref = hdr.my_type

  readu, filenum, time_in
  readu, filenum, space_in
  readu, filenum, B_ref
  data=create_struct(data, {time: time_in, space:space_in, B_ref: B_ref})
readu, filenum, next_block

ENDIF

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

data=create_struct(data, {filename: filename})
return, data

end
