;Scripts to plot things for paper
;Mostly not reuseable
pro plot_all_tests, png=png

common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
!p.background=255
window, /free
growth_dir = '/Volumes/Seagate Backup Plus Drive/DiracWhistlers/GrowthTests/'
names=['fu', 'mod1', 'mod2', 'mod3', 'mod4', 'mod5', 'mod6']

n_files =(size(names))[1]
const_1 = read_deck(dir=growth_dir, file='deck_'+names[0]+'.status')
all_consts = list(const_1)
FOR i=1, n_files -1 DO all_consts.add, read_deck(dir=growth_dir, file='deck_'+names[i]+'.status')

growth_1 = read_all_growth(growth_dir+'growth_'+names[0]+'.dat')
all_growth = list(growth_1)
FOR i=1, n_files -1 DO all_growth.add, read_all_growth(growth_dir+'growth_'+names[i]+'.dat')

colours=[0, 0, 80, 150, 230, 200, 120, 30]  
leg_text = strarr(n_files+1)
leg_text[0] = 'Run Wpar Wperp Hpar Hperp'

plot, growth_1.axes.x/const_1.wce, growth_1.data/const_1.wce, /ylog, yrange=[1e-5, 0.1], xrange=[0, 1.5], xtitle='!4x/x!3!Dce!N', ytitle='!4c/x!3!Dce!N', color=colours[1]
leg_text[1] = 'Fu et al '+string(format='(F5.2)',const_1.vtherm_par/v0)+' '+string(format='(F5.2)', const_1.vtherm_perp/v0)+' '+string(format='(F5.2)',const_1.vtherm_parh/v0)+' '+string(format='(F5.2)', const_1.vtherm_perph/v0) 

for i=0, n_files -1 DO BEGIN
  const_1 = all_consts[i]
  growth_1 = all_growth[i]
  leg_text[i+1] = 'M'+string(format='(I1)', i, /print)+' '+string(format='(F5.2)',const_1.vtherm_par/v0)+' '+string(format='(F5.2)', const_1.vtherm_perp/v0)+' '+string(format='(F5.2)',const_1.vtherm_parh/v0)+' '+string(format='(F5.2)', const_1.vtherm_perph/v0) 
  oplot, growth_1.axes.x/const_1.wce, growth_1.data/const_1.wce, color=colours[i+1]
    
END
old_charsize = !p.charsize
!p.charsize=1.0
legend, leg_text, textcolors=colours, outline_color=0, /right, /top
!p.charsize=old_charsize

if(keyword_set(png)) THEN BEGIN
  im=tvrd(true=1)
  details = {inputs:names}
  write_png_l, growth_dir+'all_tests.png', im, details=details
end

end
pro plot_time_evol, dirs, outfile=outfile, return_data=return_data, _extra=extr

base_dir = '/Volumes/Seagate Backup Plus Drive/DiracWhistlers/'

IF(N_ELEMENTS(dirs) EQ 0) THEN return
if(N_ELEMENTS(outfile) EQ 0) THEN outfile = base_dir+'time_evol_all.png'

n_dirs=(size(dirs))[1]

times_arr=list(fltarr(get_n_files(dirs[0])))
for i=1, n_dirs-1 DO times_arr.add, fltarr(get_n_files(dirs[i]))

cols = findgen(n_dirs)/float(n_dirs)*200+54
;This is how to copy a list content, rather than get a new reference to it
ens_arr = times_arr[*]
ens_arr_max= ens_arr[*]

for i=0, n_dirs-1 DO BEGIN
 ; while(file_test(base_dir+dirs[i]+'/'+string(format='(I04)', j, /print)+'.sdf')) DO BEGIN
  FOR j=0, (size(times_arr[i]))[1]-1 DO BEGIN
; FOR j=0, 5 DO BEGIN
    if(j mod 10 EQ 0) then print, j
    q=getdata(j, wkdir=base_dir+dirs[i], /ey, /ez)
    times_arr[i, j] = q.time
    ens_arr_max[i, j] = max(gauss_smooth(abs(fft(q.ey^2+q.ez^2)), 1))
    ens_arr[i, j] = total(sqrt(q.ey^2+q.ez^2))
    ;This is IDl list of arrays syntax. first index is list index, second is the array index into that list element. 
  END
END
  !p.background=255
  window, /free
  maxe = 1d-60
  foreach en, ens_arr DO maxe = max([max(en), maxe])
;  maxe=ceil(alog10(maxe))
  plot, times_arr[0], smooth(ens_arr[0], 3), xtitle='Time [s]', ytitle ='Wave Power [arb units]', color=0, yrange=[0, maxe], _extra=extr 
  for i=1, n_dirs-1 DO oplot, times_arr[i], smooth(ens_arr[i], 3), color=cols[i]
  
  im=tvrd(true=1)
  write_png_l, outfile, im, details={smooth:3}
  return_data = {times:times_arr, ens:ens_arr, dirs:dirs}
END
