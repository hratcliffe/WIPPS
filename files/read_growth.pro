function read_growth, filename

;Written for commit ID from 801d57e to ... FILL IN IF IO CHNAGES....

COMPILE_OPT IDL2
;force long ints and proper brackets

IF((N_ELEMENTS(filename) EQ 0)) THEN return, !NULL

OPENR, 1, filename

n_trials = 1
READF, 1, n_trials

str=''
;Spin through lines until
;"BEGIN"

WHILE (~EOF(1) AND ~(str EQ 'BEGIN')) DO BEGIN
  READF, 1, str
ENDWHILE

om=fltarr(n_trials)
gr=fltarr(n_trials)

count=0
WHILE (~EOF(1)) DO BEGIN

  READF, 1, omega, growth
  om[count] = omega
  gr[count] = growth
  count = count + 1
ENDWHILE

CLOSE, 1

RETURN,  {omega:om, growth:gr}

END




