
function write_data_by_num, filenum, data, axes, close_file=close_file, _extra=extr

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....
;The id string is done as exactly 10 chars ending in null. To faffy to put the null terminator within the string
;If file is not closed, make sure to 
COMPILE_OPT IDL2
;force long ints and proper brackets

sz = size(data)
n_dims = sz[0]
dims= sz[1:n_dims]
;Note sz is already LONG ints

IF((sz[-2] NE 4 && sz[-2] NE 5) || ~ISA(data, /array)) THEN BEGIN
  PRINT, "Data is not float or double array, try write_ragged"
  RETURN, 1
END
;Check we have float or double not List

IF(n_dims EQ 0) THEN RETURN, 1
;And a sane dimensionality

ret_hdr = write_header(filenum, sz[-2])

sz_sz = ret_hdr.int_sz
my_sz = ret_hdr.type_sz

next_loc = 0ull
POINT_LUN, -filenum, next_loc
next_loc = next_loc + ulong64(sz_sz*(n_dims+2))
writeu, filenum, next_loc

writeu, filenum, ulong64(n_dims)
writeu, filenum, ulong64(dims)

next_loc = next_loc + ulong64(product(dims))*my_sz+sz_sz
writeu, filenum, next_loc
writeu, filenum, data

next_loc = next_loc+sz_sz+ulong64(total(dims))*my_sz
writeu, filenum, next_loc
FOR i=0, n_dims-1 DO BEGIN
  writeu, filenum, axes[i]
END

next_loc  = next_loc + my_sz*3 + sz_sz*3
writeu, filenum, next_loc
IF(N_ELEMENTS(extr) GT 0) THEN extra_tags=tag_names(extr) else extra_tags=[]

time=[0.0*ret_hdr.my_type, 1.0]
IF(where(extra_tags EQ 'TIME') NE -1) THEN time[*] = extr.time
;Use io_check as we know it has correct type
writeu, filenum, time

space = [0*ret_hdr.block_type, 1]
IF(where(extra_tags EQ 'SPACE') NE -1) THEN space[*] = extr.space
writeu, filenum, space

IF(where(extra_tags EQ 'B_REF') NE -1) THEN b_ref = fix(extr.b_ref, type=size(ret_hdr.my_type, /type)) ELSE b_ref = 0.0*ret_hdr.my_type
writeu, filenum, b_ref

next_loc = next_loc + sz_sz + 10
writeu, filenum, next_loc

id_out = '10 chars '
IF(where(extra_tags EQ 'BLOCK') NE -1) THEN strput, id_out, extr.block, 0
writeu, filenum, id_out
null_byte = BYTE(0)
writeu, filenum, null_byte

IF(KEYWORD_SET(close_file) || where(extra_tags EQ 'CLOSE_FILE') NE -1) THEN BEGIN
  ftr_start = 0ull
  point_lun, -filenum, ftr_start
  ftr_start = ulong64(ftr_start)
  next_loc = next_loc + sz_sz + 10
  writeu, filenum, next_loc
  writeu, filenum, id_out
  writeu, filenum, null_byte
  writeu, filenum, ftr_start
END

RETURN, ret_hdr

end

function write_data, filelabel, data, axes, usenum=usenum, _extra=extr
;Write data wrapper. If usenum is set, assume filelabel is a lun. Otherwise assume it's a string name
;Note that a string-named file will always be closed
  IF(~KEYWORD_SET(usenum)) THEN BEGIN
    filenum=0
    print, filelabel
    openw, filenum, filelabel, /get_lun
    err=write_data_by_num(filenum, data, axes, /close_file, _extra=extr)
    free_lun, filenum
    return, err
  ENDIF ELSE return, write_data_by_num(filelabel, data, axes, _extra=extr)
end

