files = systemfunc("ls -1 ../../OUT_2D/TEXAS_3072x1536x72_2.5km_10s_JUL3_2025-07-05-00-00-00_*.2D_atm.nc")
print(""+files)
nfiles = dimsizes(files)

cases = (/"NONUDGE"/)

wks = gsn_open_wks("ps","ps")

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

  plot    = new((/nfiles/),graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   

  res@mpOutlineOn           = True
  res@mpGeophysicalLineColor       = "white"
  res@mpCountyLineColor       = "white"
  res@mpUSStateLineColor       = "white"
  res@mpGeophysicalLineThicknessF       = 1.
  res@mpUSStateLineThicknessF       = 3.
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..2" 

  colors = (/ (/0,0,0/), (/248,185,125/), \
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

lonmin = -105. + 360.
lonmax = -90. + 360.
latmin = 24.
latmax = 34.

end if
res@mpMinLatF              = latmin
res@mpMaxLatF              = latmax
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax


do i=0,nfiles-1
  fobs = f[i]->Precac(0,:,:)
  if(i.eq.0) then
   fobs_mean = fobs
   fobs_mean = 0.
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
  res@gsnLeftString = cases(i)
  fobs = fobs/25.4
  print("Precac: "+min(fobs)+"  "+max(fobs))
  fobs_mean = fobs_mean + fobs
  plot(i) = gsn_csm_contour_map(wks,fobs,res) 
end do


resP@gsnPanelMainString  = "Precipitation Accumulation, inch"
gsn_panel(wks,plot,(/2,2/),resP)                     ; now draw as one plot

end do ; mmm


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

fobs = f[0]->ZSFC(0,:,:)
res@cnLevelSelectionMode = "AutomaticLevels"   ; set explicit contour levels
delete(res@cnFillPalette)
res@cnFillPalette        = "BlAqGrYeOrReVi200"
res@gsnLeftString = ""
res@gsnRightString = ""
res@cnLevelSpacingF      = 50.

plot(0) = gsn_csm_contour_map(wks,fobs,res) 
resP@gsnPanelMainString  = "Model Topography, m"
gsn_panel(wks,plot(0),(/1,1/),resP)                     ; now draw as one plot

end do ; mmm
