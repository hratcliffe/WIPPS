function read_deck, dir=dir, pref=pref, file=file

;Read from deck.status into a struct of the available constants
;/todo add file exists check
COMPILE_OPT IDL2
;force long ints and proper brackets
IF(N_ELEMENTS(pref) EQ 0) THEN pref = ""
IF(N_ELEMENTS(file) GT 0) THEN BEGIN
  name = file
ENDIF ELSE BEGIN
  name = pref +'deck.status'
ENDELSE
IF(N_ELEMENTS(dir) EQ 0) THEN BEGIN
  filename = get_wkdir()+"/"+name
ENDIF ELSE BEGIN 
  filename = dir[0]+'/'+name
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
const = create_struct(const, "file", filename)
close, filenum

;Add the wpe amended by density
tags=TAG_NAMES(const)
rat = 1.0
IF( total(tags EQ 'DENS_RAT') EQ 1) THEN rat = rat + const.dens_rat
IF( total(tags EQ 'DENS_RATH') EQ 1) THEN rat = rat + const.dens_rath
const = create_struct(const, "wpe_fix", const.wpe*sqrt(rat))

RETURN, const

END
