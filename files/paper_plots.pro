;Scripts to plot things for paper
;Mostly not reuseable
pro plot_all_tests, png=png

common consts, q0, m0, v0, kb, mu0, epsilon0, h_planck
!p.background=255
window, 0
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

