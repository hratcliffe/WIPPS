function read_spect, filename

;Written for commit ID from 00db4e7 to dc9e387 FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

my_type_code = 'f'
my_type = 0.0
;this matches the type of the C code my_type...
;Usually float or double. Set code to f for float, d for double...

int_type = 1
OPENR, filenum,  filename, /GET_LUN
;open file

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

n_blocks = int_type
readu, filenum, n_blocks
IF(n_blocks NE 2) THEN BEGIN
  PRINT, "Spectrum should contain two blocks"
  RETURN, !NULL
ENDIF

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

spect = {B:read_block(filenum, my_type_code, hdr.block)}

hdr=read_header(filenum)
IF(hdr.err) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF
spect = create_struct({ang:read_block(filenum, my_type_code, hdr.block)}, spect)

return, spect
END
