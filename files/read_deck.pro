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

OPENR, filenum, filename, /GET_LUN

str=''
;Spin through lines until
;
tag = ' Constant block values after '
sz=strlen(tag)
WHILE (~EOF(filenum) AND ~(strcmp(tag, str, sz))) DO BEGIN
  READF, filenum, str
ENDWHILE

tag = 'Deck state:'
sz=strlen(tag)

;const=create_struct('blank', 0)
const=0
WHILE (~EOF(filenum) AND ~(strcmp(tag, str, sz))) DO BEGIN
  READF, filenum, str
  nv=parse_name_val(str)
  IF(~ ISA(nv, 'LIST')) THEN CONTINUE
  IF((SIZE(nv))[1] LT 2) THEN CONTINUE

  IF( ISA(const, 'struct')) THEN BEGIN
    const = create_struct(const, nv[0], nv[1])
  ENDIF ELSE BEGIN
    const = create_struct(nv[0], nv[1])
  ENDELSE

ENDWHILE

close, filenum
RETURN, const

END
