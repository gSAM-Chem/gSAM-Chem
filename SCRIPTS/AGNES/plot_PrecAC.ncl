;files = systemfunc("ls -1 ../../OUT_2D/AGNES_2048x2048x72_1.5km_10s_CONTROL_TKEDISS_PRECDISS_ENS08_9*.2D_atm.nc")
files = systemfunc("ls -1 ../../OUT_2D/AGNES_1024x1024x72_3km_10s_CONTROL_TKEDISS_PRECDISS_DELTA_JUN1812_19*.2D_atm.nc")
nfiles = dimsizes(files)

wks = gsn_open_wks("png","image")
;wks = gsn_open_wks("ps","ps")

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

lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)

res@gsnAddCyclic          = False
res@mpLimitMode            = "LatLon"
res@mpMinLatF              = min(lat)  
res@mpMaxLatF              = max(lat) 
res@mpMinLonF              = min(lon)
res@mpMaxLonF              = max(lon) 

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = False           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
  resP@gsnPanelLabelBar                 = True                    ; add common colorbar
  resP@gsnPanelFigureStringsFontHeightF = 0.015                   

  plot    = new(1,graphic)  

  resP@gsnPanelMainString  = "Precipitation"

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


fobs1 = f[0]->Precac(0,:,:)
fobs1 = 0.
do i=0,nfiles-1
  fobs = f[i]->Precac(0,:,:)
  fobs1 = fobs1 + f[i]->Prec(0,:,:)/24.
  print("Precac: "+min(fobs)+"  "+max(fobs))
  print("Precac1: "+min(fobs1)+"  "+max(fobs1))
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
  res@gsnLeftString = ""
  fobs = fobs/25.4
  fobs2 = fobs1/25.4
  copy_VarAtts(fobs, fobs2)
  copy_VarCoords(fobs, fobs2)
  plot(0) = gsn_csm_contour_map(wks,fobs,res) 

  gsn_panel(wks,plot,(/1,1/),resP)                     ; now draw as one plot
end do
