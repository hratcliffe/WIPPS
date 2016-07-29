function read_spect, filename

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....

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

error=read_header(filenum)
IF(error) THEN BEGIN
  FREE_LUN, filenum
  PRINT, "Error reading file header"
  RETURN, !NULL
ENDIF

n_blocks = int_type
readu, filenum, n_blocks
