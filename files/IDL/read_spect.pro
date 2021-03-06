function read_spect, filename

;Written for commit ID from 3b24e7d to  FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

IF(~FILE_TEST(filename)) THEN BEGIN
  PRint, "No spectra " + filename+" found"
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

hdr_end = hdr.block_type
POINT_LUN, -filenum, hdr_end
POINT_LUN, filenum, 0
int_sz=1
readu, filenum, int_sz

file_end = hdr.block_type
start_pos = read_footer_start(filenum, hdr_info=hdr, file_end=file_end)

POINT_LUN, filenum, hdr_end
;Grab int size and footer start

B = read_block(filenum, hdr.my_type, hdr.block_type)
;Read one block

next_pos = hdr.block_type
readu, filenum, next_pos

;Read in the time and space fields and the B_ref
space_in = [hdr.block_type, hdr.block_type]
time_in = [hdr.my_type, hdr.my_type]
B_ref = hdr.my_type

readu, filenum, time_in
readu, filenum, space_in
readu, filenum, B_ref
B = create_struct(B, {time: time_in, space:space_in, B_ref: B_ref})
spect = {B:B}

readu, filenum, next_block

if(next_pos EQ start_pos) THEN BEGIN
  print, "Insufficient arrays found, spectrum incomplete"
  return, spect
end

POINT_LUN, filenum, next_pos
;Skip on to the next block..
readu, filenum, next_pos
POINT_LUN, filenum, next_pos

;Now the next complete array
hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading second block header"
  RETURN, !NULL
ENDIF

ang = read_block(filenum, hdr.my_type, hdr.block_type)

;spect = create_struct({ang:read_block(filenum, hdr.my_type, hdr.block_type)}, spect)
readu, filenum, next_block
readu, filenum, next_block

;Read in the time and space fields and the B_ref
space_in = [hdr.block_type, hdr.block_type]
time_in = [hdr.my_type, hdr.my_type]
B_ref = hdr.my_type

readu, filenum, time_in
readu, filenum, space_in
readu, filenum, B_ref
ang = create_struct(ang, {time: time_in, space:space_in, B_ref: B_ref})
spect = create_struct({ang:ang}, spect)

readu, filenum, next_pos
;This should be the header start pos now. If not, print a warning then skip to what should be header start

if(next_pos NE start_pos) THEN BEGIN
  print, "Extra arrays in input file"
end

POINT_LUN, filenum, start_pos
readu, filenum, next_pos

id_type ='1234567891'
id_in = id_type
;type for block_izes and block id...

readu, filenum, id_in
spect=create_struct(spect, {block:id_in})
POINT_LUN, -filenum, next_pos

;If there is space for EXACTLY ONE int, should be smooth param
IF next_pos EQ file_end - int_sz OR next_pos EQ file_end - int_sz - 1 THEN BEGIN
  smth=hdr.block_type
  readu, filenum, smth
  spect=create_struct(spect, {smooth:smth})
END

;If there is more space then we expect the following
IF next_pos LT file_end - int_sz THEN BEGIN
  ;Read in the time and space fields and the B_ref
  space_in = [hdr.block_type, hdr.block_type]
  time_in = [hdr.my_type, hdr.my_type]
  B_ref = hdr.my_type

  readu, filenum, time_in
  readu, filenum, space_in
  spect=create_struct(spect, {time: time_in, space:space_in})

  tmp_id=hdr.block_type
  readu, filenum, tmp_id
  spect = create_struct(spect, {wave_id:tmp_id})
  readu, filenum, tmp_id
  spect = create_struct(spect, {function_type:tmp_id})
  smth=hdr.block_type
  readu, filenum, smth
  spect=create_struct(spect, {smooth:smth})

  readu, filenum, next_block

END
spect=create_struct(spect, {filename: filename})

FREE_LUN, filenum
return, spect
END
