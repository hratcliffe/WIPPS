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
free_lun, filenum

;Add the wpe amended by density
tags=TAG_NAMES(const)
rat = 1.0
IF( total(tags EQ 'DENS_RAT') EQ 1) THEN rat = rat + const.dens_rat
IF( total(tags EQ 'DENS_RATH') EQ 1) THEN rat = rat + const.dens_rath
const = create_struct(const, "wpe_fix", const.wpe*sqrt(rat))

RETURN, const

END

function read_deck_all, dir=dir, pref=pref, file=file, include_strings=include_strings

;Read all deck specs into struct
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
IF(N_ELEMENTS(include_strings) EQ 0) THEN include_strings = 0
openr, filenum, filename, /get_lun

str=''
deck_specs=0
WHILE (~EOF(filenum)) DO BEGIN
  READF, filenum, str
  IF(~strmatch(str, '*=*')) THEN CONTINUE
  ;IF(strmatch(str, '*Element*handled OK')) THEN CONTINUE
  match = stregex(str, 'Element (.*) handled OK', /extract, /subexpr)
  IF(match[0] NE "") THEN str = match[1]

  nv=parse_name_val(str, include_strings=include_strings)
  IF(~ ISA(nv, 'LIST')) THEN CONTINUE
  IF((SIZE(nv))[1] LT 2) THEN CONTINUE
  IF( ISA(deck_specs, 'struct')) THEN BEGIN
    if(where(tag_names(deck_specs) EQ strupcase(nv[0])) EQ -1) THEN deck_specs = create_struct(deck_specs, nv[0], nv[1])
  ENDIF ELSE BEGIN
    deck_specs = create_struct(nv[0], nv[1])
  ENDELSE

ENDWHILE
deck_specs = create_struct(deck_specs, "file", filename)


free_lun, filenum
return, deck_specs
end
