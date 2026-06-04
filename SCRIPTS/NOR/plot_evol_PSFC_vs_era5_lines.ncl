files = systemfunc("ls -1 ../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_2018*.2D_atm.nc")
nfiles = dimsizes(files)

g = addfile("/glade/work/marat/ERA5_Hourly/NOR/ERA5_SFC_JAN3-5_2018_NOR_hourly.nc","r")
nstart = 0 


wks = gsn_open_wks("png","image")
;wks = gsn_open_wks("ps","ps")

f=addfiles(files,"r")

lat = f[0]->lat
lon = f[0]->lon
print("lon: "+min(lon)+" "+max(lon))
print("lat: "+min(lat)+" "+max(lat))

res = True

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

resm = res

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = False           ; turn off individual cb's
res@gsnDraw             = False
res@gsnFrame            = False
res@cnFillOn            = False
res@cnLinesOn           = True
res@cnLineThicknessF    = 4.0
res@cnInfoLabelOn       = False
res@lbLabelBarOn        = False
res@cnLineLabelsOn = False                    ; turn the line labels off


;************************************************

  plot    = new(2,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnMinLevelValF       = 950.
  res@cnMaxLevelValF       = 1020.
  res@cnLevelSpacingF      = 2. 
  res@cnLineLabelBackgroundColor = -1
  res@mpOutlineOn          = True


;  res@mpOutlineBoundarySets = "AllBoundaries"
;  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
;  res@mpDataSetName = "Earth..1"

color = (/"black","red","blue","green","magenta"/)

;do i=0,nfiles-1
do i=0,0
  fld = f[i]->PSFC(0,:,:)
  print("fld="+min(fld)+"  "+max(fld))
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
  res@gsnLeftString = " "
  res@cnLineLabelsOn = False
  res@cnLineColor           = color(0)
  base = gsn_csm_map(wks, resm)
  plot(0) = gsn_csm_contour_map(wks,fld,res) 
  overlay(base, plot(0))
  delete(fld)

  fld = g->sp(nstart+i,:,:)
  fld = fld * 0.01
  print("fldobs="+min(fld)+"  "+max(fld))
  res@gsnLeftString = " "
  res@cnLineLabelsOn = False
  res@cnLineColor           = color(1)
  plot(1) = gsn_csm_contour_map(wks,fld,res) 
  overlay(base, plot(1))
  delete(fld)

  draw(base)
  frame(wks)

  delete(base)

end do
