files = systemfunc("ls -1 ../../OUT_2DL/NOR_1536x1536x72_2.5km_10s_JAN3_2018-01-05-12-00-00_0*.2D_lnd.nc")
files0 = systemfunc("ls -1 ../../OUT_2DL/NOR_1536x1536x72_2.5km_10s_JAN3_2018-01-03-12-00-00_0*.2D_lnd.nc")
print(""+files)
nfiles = dimsizes(files)

cases = (/"NUDGE 6h"/)

wks = gsn_open_wks("ps","ps")

f=addfiles(files,"r")
f0=addfiles(files0,"r")

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
;  res@mpGeophysicalLineColor       = "white"
;  res@mpUSStateLineColor       = "white"
  res@mpGeophysicalLineThicknessF       = 1.
  res@mpUSStateLineThicknessF       = 1.
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1" 

  colors = (/ (/255,255,255/), (/232,239,244/), \
              (/196,214,228/), (/132,173,209/), (/82,132,181/), \
              (/41,87,150/), (/19,54,139/), (/253,248,165/), \
              (/241,196,83/), (/230,144,70/), (/189,74,48/), \
              (/134,52,39/), (/89,29,19/), (/35,13,48/), (/65,38,114/)/) /255.
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels             = (/0.1,1.,2.,3.,4.,6.,8.,12.,18.,24.,30.,36.,48./)
  res@cnFillPalette        = colors

do mmm = 0,1

if(mmm.eq.0) then

lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)

else

lonmin = -88. + 360.
lonmax = -65. + 360.
latmin = 30.
latmax = 48.

end if
res@mpMinLatF              = latmin
res@mpMaxLatF              = latmax
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax


do i=0,nfiles-1
  fobs = f[i]->SNOWD(0,:,:)
  fobs = fobs - f0[i]->SNOWD(0,:,:)
  fobs = fobs * 39.3701 ; convert to inches.
  fobs = fobs * 1.5 ; fudge
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
  print("Precac: "+min(fobs)+"  "+max(fobs))
  fobs_mean = fobs_mean + fobs
  plot(i) = gsn_csm_contour_map(wks,fobs,res) 
end do


resP@gsnPanelMainString  = "Snow Accumulation, inches"
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
