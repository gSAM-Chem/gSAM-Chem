files = systemfunc("ls -1 ../../OUT_2D/AGNES_1024x1024x72_3km_10s_CONTROL_TKEDISS_PRECDISS_1972-06-24-00-00-00_*.2D_atm.nc ../../OUT_2D/AGNES_1024x1024x72_3km_10s_CONTROL_TKEDISS_PRECDISS_ENS*1972-06-24-00-00-00_*.2D_atm.nc")
print(""+files)
nfiles = dimsizes(files)

case = "1812"

wks = gsn_open_wks("ps","ps")

if(ens.eq.ens_max-1) then
  ncfiles=systemfunc("ls -1 "+filename+"_19*.2D_atm.nc")
else
  ncfiles=systemfunc("ls -1 "+filename+"_ENS"+nens+"_19*.2D_atm.nc")
end if

f=addfiles(files,"r")

lat = f[0]->lat
lon = f[0]->lon
print("lon: "+min(lon)+" "+max(lon))
print("lat: "+min(lat)+" "+max(lat))

res = True
res@cnFillOn              = True
res@cnLinesOn             = False
res@cnLineLabelsOn = False                    ; turn the line labels off
res@cnFillMode            = "RasterFill"

res@gsnAddCyclic          = False
res@mpLimitMode            = "LatLon"

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = False           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
  resP@gsnPanelLabelBar                 = True                    ; add common colorbar
  resP@gsnPanelFigureStringsFontHeightF = 0.015                   
  resP@gsnMaximize = True  

  plot    = new((/nfiles+1/),graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   

  resP@gsnPanelMainString  = case+"  Precipitation Accumulation, inch"

  res@mpOutlineOn           = True
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1" 

  colors = (/ (/255,255,255/), (/248,185,125/), \
              (/243,120,41/), (/240,47,29/), (/120,24,15/), \
              (/126,58,226/), (/24,68,124/)/) / 255.
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels             = (/1.,3.,5.,7.,10.,15./)
  res@cnFillPalette        = colors

do mmm = 0,1

if(mmm.eq.0) then
  
lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)

else

lonmin = -83. + 360.
lonmax = -75. + 360.
latmin = 36.
latmax = 40.

end if
res@mpMinLatF              = latmin  
res@mpMaxLatF              = latmax 
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax 

do i=0,nfiles-1
  fld = f[i]->Precac(0,{latmin:latmax},{lonmin:lonmax})
  if(i.eq.0) then
   fld_mean = fld
   fld_mean = 0.
  end if
  time = f[i]->time
  dates = tointeger(cd_calendar(time, 0))
  year = dates(0,0)
  month = sprinti("%0.2d", dates(0,1))
  day = sprinti("%0.2d", dates(0,2))
  mins = sprinti("%0.2d", dates(0,4))
  if(mins.eq.59) then
   dates(0,3) = dates(0,3)+1
  end if
  hour = sprinti("%0.2d", dates(0,3))
  print(year + "-" + month + "-" + day + "-" + hour)
  res@gsnRightString = year + "-" + month + "-" + day + "-" + hour
  if(i.eq.0) then
    res@gsnLeftString = "ENS 0"
  else
    res@gsnLeftString = "ENS "+i
  end if
  fld = fld/25.4
  print("Precac: "+min(fld)+"  "+max(fld))
  fld_mean = fld_mean + fld
  plot(i) = gsn_csm_contour_map(wks,fld,res) 
end do

;===============================
; ERA5 Precip

fprec_obs = "/glade/work/marat/ERA5_Hourly/AGNES/AGNES_SFC.nc"
g=addfile(fprec_obs,"r")
latobs = g->latitude
lonobs = g->longitude
print("lonobs: "+(min(lonobs)+360.)+" "+(max(lonobs)+360.))
print("latobs: "+min(latobs)+" "+max(latobs))
if(any(lonobs.lt.0.)) then
   lonmin = lonmin-360.
   lonmax = lonmax-360.
   print("lonmin/max: "+lonmin+"  "+lonmax) 
end if
fobs1 = g->avg_tprate(12:132,{latmin:latmax},{lonmin:lonmax})
fobs = dim_sum_n_Wrap(fobs1,0)
fobs = fobs * 3600./25.4
print("ERA5 prec: "+min(fobs)+"  "+max(fobs))

;================================================

;res@gsnLeftString = "ENS MEAN"
;fld_mean = fld_mean/tofloat(nfiles)
;print("Precac (Ens mean): "+min(fld_mean)+"  "+max(fld_mean))

res@gsnLeftString = "ERA5"
plot(nfiles) = gsn_csm_contour_map(wks,fobs,res) 

gsn_panel(wks,plot,(/3,4/),resP)                     ; now draw as one plot

delete(fld)
delete(fld_mean)
delete(fobs)
delete(fobs1)

end do ; mmm
