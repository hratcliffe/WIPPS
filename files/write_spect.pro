function write_spect, filename, spectrum

;Written for commit ID from 3b24e7d to  FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

OPENW, filenum, filename, /GET_LUN
;open file

ax_list = list(spectrum.B.axes.x)
ret = write_data(filenum, spectrum.B.data, ax_list, space=spectrum.B.space, time=spectrum.B.time, b_ref=spectrum.B.b_ref, /usenum)

ax_list = list(spectrum.ang.axes.x, spectrum.ang.axes.y)
ret = write_data(filenum, spectrum.ang.data, ax_list, space=spectrum.ang.space, time=spectrum.ang.time, b_ref=spectrum.ang.b_ref, /usenum, /close_file)

FREE_LUN, filenum
return, 0
END
