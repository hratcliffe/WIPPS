function read_spectra_from_file, dir, spec_file

  IF(~FILE_TEST(dir+spec_file)) THEN BEGIN
    PRINT, 'Spectra file '+dir+spec_file+' not found'
    RETURN, !NULL
  END
  openr, filenum, dir+spec_file, /get_lun
  tmp='a'
  readf, filenum, tmp
  names = [tmp]
  while(~EOF(filenum)) DO BEGIN
    readf, filenum, tmp
    names=[names, tmp]
  endwhile

  FOR i=0, (size(names))[1]-1 DO BEGIN
    tmp=read_spect(dir+strtrim(names[i]))
    if(i EQ 0 && N_ELEMENTS(tmp) GT 0) THEN all_spect = [tmp] ELSE all_spect = [all_spect, [tmp]]
  END
return, all_spect
end
