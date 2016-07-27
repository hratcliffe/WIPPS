function read_deck, dir=dir, pref=pref

;Read from deck.status into a struct of the available constants

COMPILE_OPT IDL2
;force long ints and proper brackets
IF(N_ELEMENTS(pref) EQ 0) THEN pref = ""
IF(N_ELEMENTS(dir) EQ 0) THEN BEGIN
  filename = get_wkdir()+"/"+pref+"deck.status"
ENDIF ELSE BEGIN 
  filename = dir[0]+'/'+pref+'deck.status'
ENDELSE

OPENR, 1, filename

str=''
;Spin through lines until
;
tag = ' Constant block values after '
sz=strlen(tag)
WHILE (~EOF(1) AND ~(strcmp(tag, str, sz))) DO BEGIN
  READF, 1, str
ENDWHILE

tag = 'Deck state:'
sz=strlen(tag)

;const=create_struct('blank', 0)
const=0
WHILE (~EOF(1) AND ~(strcmp(tag, str, sz))) DO BEGIN
  READF, 1, str
  nv=parse_name_val(str)
  if((SIZE(nv))[1] LT 2) THEN CONTINUE
  IF( ISA(const, 'struct')) THEN BEGIN
    const = create_struct(const, nv[0], DOUBLE(nv[1]))
  ENDIF ELSE BEGIN
    const = create_struct(nv[0], DOUBLE(nv[1])) 
  ENDELSE

ENDWHILE

close, 1
RETURN, const

END

function parse_name_val, str
;Actually very simple but make function in case of complication eh?
  
  name_val =strtrim(strsplit(str, '=', /extract), 2)

  return, name_val
END

