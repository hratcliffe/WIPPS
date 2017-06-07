function read_footer_start, filenum, hdr_info=hdr_info, file_end=file_end, file_sz=file_sz

tmp = FSTAT(filenum)
POINT_LUN, filenum, (tmp.size - hdr_info.int_sz)
file_sz = tmp.size
file_end = tmp.size - hdr_info.int_sz
start_pos = hdr_info.block_type
readu, filenum, start_pos

;If start_pos is larger than the file_size something has gone wrong. Most likely that is a trailing newline char. We try reading from 1 char back and see if that makes sense
IF start_pos GT tmp.size THEN BEGIN
  POINT_LUN, filenum, (tmp.size - hdr_info.int_sz - 1)
  file_end = tmp.size - hdr_info.int_sz
  start_pos = hdr_info.block_type
  readu, filenum, start_pos
END
;If it's still bad, fail
IF start_pos GT tmp.size THEN BEGIN
  PRINT, "Footer position invalid"
  RETURN, -1
END

RETURN, start_pos
END
