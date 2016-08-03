function read_header, filenum, report=report

;Written for commit ID from dc9e387 to ... FILL IN IF IO CHNAGES....
;* (int)sizeof(size_t) (int)sizeof(my_type) io_verification_code Version string
;*Next_block n_dims dims[n_dims]
;*Next_block data

COMPILE_OPT IDL2
;force long ints and proper brackets

commit_type ='123456789112345'
commit_in = commit_type
;type for commit id

hdr_info = {err:1}

;Check file is good
fileinf=FSTAT(filenum)
IF(NOT fileinf.open) THEN BEGIN
  PRINT, "Invalid file handle"
  RETURN, hdr_info
ENDIF

;Read number sizes
sz_sz=0
my_sz=0
readu, filenum, sz_sz
readu, filenum, my_sz

IF(KEYWORD_SET(report)) THEN PRINT, sz_sz, my_sz
IF(my_sz EQ 4) THEN BEGIN
  hdr_info =create_struct(hdr_info,{my_type:0.0})
ENDIF ELSE IF(my_sz EQ 8) THEN BEGIN
  hdr_info =create_struct(hdr_info,{my_type:0.0d})
ENDIF ELSE BEGIN
  return, hdr_info
END

IF(sz_sz EQ 2) THEN BEGIN
  hdr_info =create_struct(hdr_info,{block_type:0U})
ENDIF ELSE IF(sz_sz EQ 4) THEN BEGIN
  hdr_info =create_struct(hdr_info,{block_type:0UL})
ENDIF ELSE IF(sz_sz EQ 8) THEN BEGIN
  hdr_info =create_struct(hdr_info,{block_type:0ULL})
ENDIF ELSE BEGIN
  return, hdr_info
END

io_type =hdr_info.my_type
io_check = 3.0/32.0
io_in=io_type
;type for io verification const.

next_block = hdr_info.block_type

readu, filenum, io_in
readu, filenum, commit_in
IF(KEYWORD_SET(report)) THEN print, io_in," ",  commit_in

IF(io_in NE io_check) THEN BEGIN
  print, "File read error!"
ENDIF ELSE BEGIN
  hdr_info.err = 0
ENDELSE
return, hdr_info

END
