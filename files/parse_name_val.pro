function parse_name_val, str, delim=delim
;Actually very simple but make function in case of complication eh?

  IF(n_elements(delim) EQ 0) THEN delim='='
  IF(strpos(str, delim) EQ -1) THEN RETURN, !NULL
  name_val =strtrim(strsplit(str, delim, /extract), 2)

  IF((size(name_val))[1] GE 2) THEN BEGIN
    ;We have two fields, try making the second into number
    ;Default is double, but if no decimal point, make long

    IF( STRPOS(name_val[1], '.') EQ -1) THEN BEGIN
      val = LONG(name_val[1])
    ENDIF ELSE BEGIN
      val = DOUBLE(name_val[1])
    ENDELSE
  return, LIST(name_val[0], val)
  ;Return first value as string and second as correct type

  ENDIF ELSE BEGIN
    RETURN, !NULL
  ENDELSE
END

