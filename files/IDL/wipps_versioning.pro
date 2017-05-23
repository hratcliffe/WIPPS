

function read_wipps_version, filename
;Read the wipps version string from file

  openr, filenum, filename, /GET_LUN
  hdr=read_header(filenum)
  return, hdr.version
end

function extract_version, str
;Extract the version part, vx.y from str and return struct of {major:x, minor:y, err:false}

  match = stregex(str, 'v([0123456789]*).([0123456789]*)(.*)', /extract, /subexpr)

  IF (size(match))[1] GT 2  AND match[0] NE "" THEN BEGIN
    return, {major:fix(match[1]), minor:fix(match[2]), err:0}
  ENDIF ELSE BEGIN
    return, {major:-1, minor:-1, err:1}
  ENDELSE
end

function compare_as_version, str, vers_str, check_minor=check_minor, force_numeric = force_numeric
;Compare strings as version strings.
;If check_minor is set, we check minor version also
;As in the cpp code, returns 0 if equal, -1 if str is before vers_str, 1 if str is after vers_str unless strings don't have numeric version parts in which case return 0 for equal, 1 for unequal

  in_vers = extract_version(str)
  comp_vers = extract_version(vers_str)

  IF(~KEYWORD_SET(force_numeric) AND (in_vers.err EQ 1 OR comp_vers.err EQ 1)) THEN BEGIN
    return, ~strcmp(str, vers_str)
  ENDIF
  IF(KEYWORD_SET(check_minor)) THEN BEGIN
    IF(in_vers.major GT comp_vers.major) THEN RETURN, 1
    IF(in_vers.major LT comp_vers.major) THEN RETURN, -1
    ;Now major vers are equal so:
    IF(in_vers.minor GT comp_vers.minor) THEN RETURN, 1
    IF(in_vers.minor LT comp_vers.minor) THEN RETURN, -1
    RETURN, 0
  ENDIF ELSE BEGIN
    IF(in_vers.major GT comp_vers.major) THEN RETURN, 1
    IF(in_vers.major LT comp_vers.major) THEN RETURN, -1
    RETURN, 0
  ENDELSE

end
