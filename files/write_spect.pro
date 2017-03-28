function write_spect, filename, spectrum

;Written for commit ID from 3b24e7d to  FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

OPENW, filenum, filename, /GET_LUN
;open file

ax_list = list(spectrum.B.axes.x)
ret_hdr = write_data(filenum, spectrum.B.data, ax_list, space=spectrum.B.space, time=spectrum.B.time, b_ref=spectrum.B.b_ref, /usenum)

ax_list = list(spectrum.ang.axes.x, spectrum.ang.axes.y)
ret_hdr = write_data(filenum, spectrum.ang.data, ax_list, space=spectrum.ang.space, time=spectrum.ang.time, b_ref=spectrum.ang.b_ref, /usenum)

next_loc = ret_hdr.block_type
ftr_start = ret_hdr.block_type

POINT_LUN, -filenum, ftr_start
ftr_start = fix(ftr_start, type=size(ret_hdr.block_type, /type))
next_loc  = ftr_start + ret_hdr.type_sz*2 + ret_hdr.int_sz*6 + 10
writeu, filenum, next_loc

id_out = 'IDL spect'
null_byte = BYTE(0)
writeu, filenum, id_out
writeu, filenum, null_byte

writeu, filenum, fix(spectrum.B.time[0], type=size(ret_hdr.my_type, /type))
writeu, filenum, fix(spectrum.B.time[1], type=size(ret_hdr.my_type, /type))

writeu, filenum, fix(spectrum.B.space[0], type=size(ret_hdr.block_type, /type))
writeu, filenum, fix(spectrum.B.space[1], type=size(ret_hdr.block_type, /type))

writeu, filenum, ret_hdr.block_type
writeu, filenum, ret_hdr.block_type
writeu, filenum, ret_hdr.block_type ;Other stuff that we don't have in the IDL'

writeu, filenum, ftr_start

FREE_LUN, filenum
return, 0
END
