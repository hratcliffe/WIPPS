function read_constants_dump, dir=dir, pref=pref

;Read from constants.dump into a struct of the available constants

COMPILE_OPT IDL2
;force long ints and proper brackets
IF(N_ELEMENTS(pref) EQ 0) THEN pref = ""
IF(N_ELEMENTS(dir) EQ 0) THEN BEGIN
  filename = get_wkdir()+"/"+pref+"constants.dump"
ENDIF ELSE BEGIN 
  filename = dir[0]+'/'+pref+'constants.dump'
ENDELSE

OPENR, filenum, filename, /GET_LUN

const=0
str=''
WHILE (~EOF(filenum)) DO BEGIN
  READF, filenum, str
  nv=parse_name_val_space(str)
  if((SIZE(nv))[1] LT 2) THEN CONTINUE
  IF( ISA(const, 'struct')) THEN BEGIN
    const = create_struct(const, nv[0], DOUBLE(nv[1]))
  ENDIF ELSE BEGIN
    const = create_struct(nv[0], DOUBLE(nv[1])) 
  ENDELSE

ENDWHILE

close, filenum
RETURN, const

END

function parse_name_val_space, str
;Actually very simple but make function in case of complication eh?
  
  name_val =strtrim(strsplit(str, ' ', /extract), 2)

  return, name_val
END

