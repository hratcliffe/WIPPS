function read_padie_data, dir, file_start = file_start, n = n

COMPILE_OPT IDL2
;force long ints and proper brackets

IF(N_ELEMENTS(dir) EQ 0) THEN dir = './'
IF(N_ELEMENTS(file_start) EQ 0) THEN file_start = "Figure_2_fpefce_2.5"
n_part = ""
IF(N_ELEMENTS(n) NE 0) THEN BEGIN
  n_part = '_n_' + strtrim(string(n, /print), 2)
ENDIF ELSE n = 0

filename = dir + file_start + n_part +'.out'

IF(~FILE_TEST(filename)) THEN BEGIN
  PRINT, "File " + filename+" not found"
  RETURN, !NULL
END

OPENR, filenum,  filename, /GET_LUN
;open file


;This extracts the integer part from a string like "Number of pitch angles:  100"
tag = 'Number of pitch angles:'
sz = strlen(tag)
str=''
READF, filenum, str
str = strtrim(str, 1)
WHILE (~EOF(filenum) AND ~(strcmp(tag, str, sz))) DO BEGIN
  READF, filenum, str
  str = strtrim(str, 1)
ENDWHILE
;Attempt to read an integer from this line
str_secs = strsplit(str, ':', /extract)
n_angs = fix(str_secs[1])

;Spin on until this which precedes data
tag = 'Pitch-Angle'
sz=strlen(tag)
WHILE (~EOF(filenum) AND ~(strcmp(tag, str, sz))) DO BEGIN
  READF, filenum, str
  str = strtrim(str, 1)
ENDWHILE
columns = (strsplit(strcompress(str), ' ', /extract))[0:5] ;Extract titles from 6 columns

;Now the next lines should be data, 6 columns, n_angs lines
;Double data
data = dblarr(6, n_angs)
READF, filenum, data

;Particle velocity:  2.821E+08 m/s
;Spin on for a bit more info, assumed in order
tag = 'Particle velocity'
sz=strlen(tag)
WHILE (~EOF(filenum) AND ~(strcmp(tag, str, sz))) DO BEGIN
  READF, filenum, str
  str = strtrim(str, 1)
ENDWHILE
;Extract particle velocity part
vel = float((strsplit((strsplit(str, ':', /extract))[1], ' ', /extract))[0])

free_lun, filenum
struc = {n:n, axis:reform(data[0, *]), columns:strjoin(columns, ' '), particle_velocity: vel}
FOR i = 1, 5 DO BEGIN
  struc = create_struct(struc, strmid(columns[i], 0, 3), reform(data[i, *]))
END
struc = create_struct(name='PadieStruct', struc)
return, struc

end