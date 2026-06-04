files = systemfunc("ls -1 ../../OUT_2D/AGNES_1024x1024x72_3km_10s_CONTROL_TKEDISS_PRECDISS_*1972-06-24-00-00-00_0000047520.2D_atm.nc")
print(""+files)
nfiles = dimsizes(files)

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
  resP@gsnMaximize = True  

  plot    = new((/nfiles+1/),graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   

  resP@gsnPanelMainString  = "Precipitation Accumulation, inches"

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
  if(i.eq.0) then
    res@gsnLeftString = "REANALYSIS"
  else
    res@gsnLeftString = "ENS "+i
  end if
  fobs = fobs/25.4
  print("Precac: "+min(fobs)+"  "+max(fobs))
  fobs_mean = fobs_mean + fobs
  plot(i) = gsn_csm_contour_map(wks,fobs,res) 
end do

res@gsnLeftString = "ENS MEAN"
fobs_mean = fobs_mean/tofloat(nfiles)
print("Precac (Ens mean): "+min(fobs_mean)+"  "+max(fobs_mean))
plot(nfiles) = gsn_csm_contour_map(wks,fobs_mean,res) 

gsn_panel(wks,plot,(/3,4/),resP)                     ; now draw as one plot
