function build_mask, ax_x, ax_y, ax_z, fzz=fzz, om_ce=om_ce, om_pe=om_pe

sz=[(size(ax_x))[1], (size(ax_y))[1], (size(ax_z))[1]]
mask=intarr(sz)
for i=0, sz[0]-1 DO BEGIN
  for j=0, sz[1]-1 DO BEGIN
    modk=sqrt(ax_x[i]^2+ax_y[j]^2)
    om=get_dispersion(modk, om_ce=om_ce, om_pe=om_pe, /appr)
    om_min=min(where(ax_z GT om*(1.0-fzz/100.0)))
    om_max=min(where(ax_z GT om*(1.0+fzz/100.0)))
    if(om_min EQ -1 OR om_max EQ -1) then continue
  for k=om_min, om_max DO mask[i, j, k] = 1 
  end
end
return, mask
end

