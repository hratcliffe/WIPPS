
function get_array_of_type, type_code, sz
;Follows http://www.harrisgeospatial.com/docs/IDL_Data_Types.html
  IF(type_code EQ 2) THEN BEGIN return,  INTARR(sz)
  ENDIF ELSE IF(type_code EQ 12) THEN BEGIN return, UINTARR(sz)
  ENDIF ELSE IF(type_code EQ 3) THEN BEGIN return, LONARR(sz)
  ENDIF ELSE IF(type_code EQ 13) THEN BEGIN return, ULONARR(sz)
  ENDIF ELSE IF(type_code EQ 14) THEN BEGIN return, LON64ARR(sz)
  ENDIF ELSE IF(type_code EQ 15) THEN BEGIN return, ULON64ARR(sz)
  ENDIF
end

function read_block, filenum, my_type, block_type
;data = read_block(filenum, hdr.block, hdr.my_type, hdr.block_type)
;\todo More than 2-dimension???
;*Next_block n_dims dims[n_dims]
;*Next_block data
;* Next_block axes
;* Next_block Block_id

;Written for commit ID from 801d57e to dc9e387 FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

;Check file is good
fileinf=FSTAT(filenum)
IF(NOT fileinf.open) THEN BEGIN
  PRINT, "Invalid file handle"
  RETURN, 1
ENDIF

next_block = block_type
;type for block_sizes

n_dims = block_type

readu, filenum, next_block
readu, filenum, n_dims

if(n_dims GT 0) THEN BEGIN

  dims = get_array_of_type(SIZE(block_type, /TYPE), n_dims)
  readu, filenum, dims
  readu, filenum, next_block
  IF(SIZE(my_type, /TYPE) EQ 4) THEN BEGIN
    axes_list = {x:fltarr(dims[0])}
    if(n_dims GT 1) THEN axes_list = create_struct(axes_list, {Y:fltarr(dims[1])})
    if(n_dims GT 2) THEN axes_list = create_struct(axes_list, {Z:fltarr(dims[2])})

    data = {data:fltarr(dims), axes:axes_list}
    tmp2=fltarr(dims)

  ENDIF ELSE IF(SIZE(my_type, /TYPE) EQ 5) THEN BEGIN
    axes_list = {x:dblarr(dims[0])}
    if(n_dims GT 1) THEN axes_list = create_struct(axes_list, {Y:dblarr(dims[1])})
    if(n_dims GT 2) THEN axes_list = create_struct(axes_list, {Z:dblarr(dims[2])})

    data = {id:id_in,data:dblarr(dims), axes:axes_list}
    tmp2=dblarr(dims)

  ENDIF ELSE BEGIN
    print, "Unknown data type"
    RETURN, !NULL
  ENDELSE

  ;Can't read directly into anon structure field, so use tmp
  readu, filenum, tmp2
  data.data = tmp2
  tmp2=0
ENDIF ELSE BEGIN
  print, "Array is ragged. Use read_ragged.pro"
  RETURN, !NULL
ENDELSE

readu, filenum, next_block


IF(SIZE(my_type, /TYPE) EQ 4) THEN BEGIN
  tmpa=fltarr(dims[0])
  readu, filenum, tmpa
  data.axes.X = tmpa

  IF(n_dims GT 1) THEN BEGIN
    tmpa=fltarr(dims[1])
    readu, filenum, tmpa
    data.axes.Y = tmpa
  ENDIF
  IF(n_dims GT 2) THEN BEGIN
    tmpa=fltarr(dims[2])
    readu, filenum, tmpa
    data.axes.Z = tmpa
  ENDIF

ENDIF ELSE IF(SIZE(my_type, /TYPE) EQ 5) THEN BEGIN
  tmpa=dblarr(dims[0])
  readu, filenum, tmpa
  data.axes.X = tmpa

  IF(n_dims GT 1) THEN BEGIN
    tmpa=dblarr(dims[1])
    readu, filenum, tmpa
    data.axes.Y = tmpa
  ENDIF
  IF(n_dims GT 2) THEN BEGIN
    tmpa=dblarr(dims[2])
    readu, filenum, tmpa
    data.axes.Z = tmpa
  ENDIF
ENDIF

tmpa = 0

RETURN, data
END