pro setup_fftw

COMPILE_OPT IDL2
;force long ints and proper brackets

FFTW_dir = '/Users/heatherratcliffe/FFTW_IDL/IDLfftw3'
!PATH =!PATH + ':'+FFTW_dir

cd, current = tmp_dir
cd, FFTW_dir
link_fftw
cd, tmp_dir

test_idlfftw3


end
