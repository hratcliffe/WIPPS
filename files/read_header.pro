function read_header, filenum

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

;Check file is good
fileinf=FSTAT(filenum)
IF(NOT fileinf.open) THEN BEGIN
  PRINT, "Invalid file handle"
  RETURN, 1
ENDIF

id_type ='1234567891'
id_in = id_type
;type for block id...

io_type =0.0
io_check = 3.0/32.0
io_in=io_type
;type for io verification const.

commit_type ='123456789112345'
commit_in = commit_type
;type for commit id

readu, filenum, id_in
readu, filenum, io_in
readu, filenum, commit_in
print, id_in," ",  io_in," ",  commit_in

IF(io_in NE io_check) THEN BEGIN
  print, "File read error!"
  return, 1
ENDIF

return, 0
END
