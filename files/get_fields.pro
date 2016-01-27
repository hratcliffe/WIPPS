
pro get_fields, start, stp, dir=dir, n_zeros=n_zeros, ex=ex, ey=ey, ez=ez, bx=bx, by=by, bz=bz, skip=skip

;Created by H Ratcliffe, University of Warwick
;16/4/2015
;Reads fields from EPOCH data dumps using supplied SDF read routines 
;
;-------------Params----------------
; start: first file number to read
; stp: last file number to read
;...........optional params......
; dir: directory to read from, otherwise uses current
; n_zeros: manually set number of zeros, otherwise is determined from the filenames
; skip: cadence of file reading
; ex, ey, ez, bx, by, bz: arrays to contain data. Only those requested will be read, and only if they are present in the file
;Call example
; get_fields, 0, 500, ex=ex, ey=ey, bz=bz, skip=5

COMPILE_OPT idl2

IF(N_ELEMENTS(dir) GT 0) THEN BEGIN
  old_dir= get_wkdir()
  set_wkdir, dir
ENDIF
;set working directory

IF(N_ELEMENTS(skip) EQ 0) THEN skip =1
;set increment of file read

IF(N_ELEMENTS(n_zeros) EQ 0) THEN BEGIN
  n_zeros=-1
  FOR i=3, 9 DO BEGIN
    fmt = '(I0'+string(format='(I01)', i, /print)+')'
    if(STRCMP(string(format=fmt, start, /print), '********', i) EQ 1) THEN CONTINUE
    ;check we have enough digits to print start
    success=file_test(get_wkdir()+string(format=fmt, start,/print) +'.sdf')
    IF(success EQ 1) THEN BEGIN
      n_zeros=i
      BREAK
    ENDIF
  END
  ;work out number of zeros in data files by testing for file numbered start

  if(n_zeros EQ -1) THEN BEGIN
    print, 'Cannot set file name length'
    return
  endif
  ;Fail if no first file found
ENDIF

len=(stp-start)/skip +1
q=getdata(start, n_zeros=n_zeros)
;get first data file for sizes and tags

strs = ['EX', 'EY', 'EZ', 'BX', 'BY', 'BZ']
get_vars = 0
tags = TAG_NAMES(q)
args = [ARG_PRESENT(ex), ARG_PRESENT(ey),ARG_PRESENT(ez),ARG_PRESENT(bx),ARG_PRESENT(by),ARG_PRESENT(bz)]

FOR i=0, N_ELEMENTS(args)-1 DO BEGIN
  bit = 2^i
  IF(args[i] EQ 1) THEN IF(where(tags EQ strs[i]) NE -1) THEN get_vars = get_vars+bit ELSE print, 'Warning:', strs[i], ' not present in file'
END

if(get_vars EQ 0) THEN BEGIN
  print, 'No data requested!'
  return
ENDIF
;set flags to get only asked for fields. These are arrays for return so must be passed by reference, hence ARG_PRESENT returns 1. We also check variable is available in file and print warning if not
  

ind = floor(alog(get_vars AND -get_vars)/alog(2.0))
sz = SIZE(q.(where(tags EQ strs[ind])))
;set sizes according to data in file, using least significant bit in get_vars for tag

sz_out = [sz[1:sz[0]], len]
len_offset=sz[-1]

IF((get_vars AND 1) EQ 1) THEN ex=fltarr(sz_out)
IF((get_vars AND 2) EQ 2) THEN ey=fltarr(sz_out)
IF((get_vars AND 4) EQ 4) THEN ez=fltarr(sz_out)
IF((get_vars AND 8) EQ 8)  THEN bx=fltarr(sz_out)
IF((get_vars AND 16) EQ 16)  THEN by=fltarr(sz_out)
IF((get_vars AND 32) EQ 32)  THEN bz=fltarr(sz_out)
;resize arrays for output

tmp_get = get_vars
FOR i=0, N_ELEMENTS(args)-1 DO BEGIN
  bit = 2^i
  IF((get_vars AND bit) EQ bit) THEN BEGIN
    IF(tmp_get EQ get_vars) THEN extr=CREATE_STRUCT(strs[i], 1) ELSE extr=CREATE_STRUCT(extr, strs[i], 1)
    tmp_get = tmp_get - bit
  ENDIF
END
;print, extr
;assemble extra struct for call to getdata

FOR i=start, stp, skip DO BEGIN

  j = (i-start)/skip

  q=getdata(i, n_zeros=n_zeros, _extra=extr)
  ;read only requested fields

  IF((get_vars AND 1) EQ 1) THEN ex[len_offset*j:len_offset*(j+1)-1]=q.ex[0:len_offset-1]
  ;read into array without knowing on code-writing how many dimensions it has, only that time is last

  IF((get_vars AND 2) EQ 2) THEN ey[len_offset*j:len_offset*(j+1)-1]=q.ey[0:len_offset-1]
  IF((get_vars AND 4) EQ 4) THEN ez[len_offset*j:len_offset*(j+1)-1]=q.ez[0:len_offset-1]
  IF((get_vars AND 8) EQ 8) THEN bx[len_offset*j:len_offset*(j+1)-1]=q.bx[0:len_offset-1]
  IF((get_vars AND 16) EQ 16) THEN by[len_offset*j:len_offset*(j+1)-1]=q.by[0:len_offset-1]
  IF((get_vars AND 32) EQ 32) THEN bz[len_offset*j:len_offset*(j+1)-1]=q.bz[0:len_offset-1]

  if(i mod 100 EQ 0) then print, 'Read file  '+ string(i, format=fmt, /print)+'.sdf'
END

IF(N_ELEMENTS(old_dir) GT 0) THEN set_wkdir, old_dir
;reset wk_dir

CATCH, Err

IF(Err NE 0) THEN print, 'Oopsie!'

END

