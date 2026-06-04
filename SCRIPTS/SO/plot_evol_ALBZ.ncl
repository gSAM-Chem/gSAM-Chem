files = systemfunc("ls -1 ../../OUT_MISC/SO_3072x3072x72_1km_10s_M2005_MISC_2018*.2D_atm.nc")
nfiles = dimsizes(files)


wks = gsn_open_wks("png","plot")
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

res@lbLabelBarOn        = False           ; turn off individual cb's

;************************************************
  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnMinLevelValF       = 0.05
  res@cnMaxLevelValF       = 0.75
  res@cnLevelSpacingF      = 0.02 

  res@cnFillPalette       = "MPL_bone"

do i=0,nfiles-1
  fld = f[i]->ALBZ(0,:,:)
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
  res@gsnLeftString = "gSAM"
  plot = gsn_csm_contour_map(wks,fld,res) 

end do
