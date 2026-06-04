files = systemfunc("ls -1 ../../OUT_2D/FIONA_1024x1024x72_2km_10s_JUN1612*2022-09-21-00-00-00_*.2D_atm.nc")
files1 = systemfunc("ls -1 ../../OUT_2D/FIONA_1024x1024x72_2km_10s_JUN1612*2022-09-18-00-00-00_*.2D_atm.nc")
print(""+files)
print(""+files1)
nfiles = dimsizes(files)

cases = (/"1612", "1612 NUDGEUV 12h","1612 NUDGEUV 6h"/)


wks = gsn_open_wks("ps","ps")

f=addfiles(files,"r")
f1=addfiles(files1,"r")

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
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..4" 

  colors = (/ (/0,0,0/),(/247,214,132/), (/239,135,50/), \
              (/234,51,36/), (/139,27,17/), (/133,90,198/), \
              (/38,81,133/), (/76,169,239/),(/0,26,245/),(/158,251,78/),(/80,125,40/)/) / 255.
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels             = (/1.,3.,5.,7.,10.,15.,20.,25.,30.,35./)
  res@cnFillPalette        = colors

do mmm = 0,1

if(mmm.eq.0) then

lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)

else

lonmin = -68. + 360.
lonmax = -65. + 360.
latmin = 17.
latmax = 20.

end if
res@mpMinLatF              = latmin
res@mpMaxLatF              = latmax
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax


do i=0,nfiles-1
  fobs = f[i]->Precac(0,:,:)
  fobs1 = f1[i]->Precac(0,:,:)
  fobs = fobs - fobs1
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
gsn_panel(wks,plot,(/2,3/),resP)                     ; now draw as one plot

end do ; mmm


do mmm = 0,1

fobs = f[0]->ZSFC(0,:,:)

delete(res@cnFillPalette)
res@cnFillPalette        = "BlAqGrYeOrReVi200"
res@gsnLeftString = ""
res@gsnRightString = ""

if(mmm.eq.0) then

lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)
res@cnLevelSelectionMode = "AutomaticLevels"   ; set explicit contour levels
res@cnLevelSpacingF      = 50.

else

lonmin = -68. + 360.
lonmax = -65. + 360.
latmin = 17.
latmax = 20.
res@cnLevelSelectionMode =  "ManualLevels"
res@cnMinLevelValF       = 0.
res@cnMaxLevelValF       = 1000.
res@cnLevelSpacingF      = 50.


end if
res@mpMinLatF              = latmin
res@mpMaxLatF              = latmax
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax


plot(0) = gsn_csm_contour_map(wks,fobs,res) 
resP@gsnPanelMainString  = "Model Topography, m"
gsn_panel(wks,plot(0),(/1,1/),resP)                     ; now draw as one plot

end do ; mmm
