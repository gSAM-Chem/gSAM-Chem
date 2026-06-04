files = "/glade/work/marat/ERA5_Hourly/AGNES/AGNES_SP.nc"
f = addfile(files,"r")
nfiles = dimsizes(f->valid_time)

wks = gsn_open_wks("png","image")
;wks = gsn_open_wks("ps","ps")

f=addfile(files,"r")

lat = f->latitude
lon = f->longitude
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

  resP@gsnPanelMainString  = "ERA5 Surface Pressure Perturbation, mb"

  res@mpOutlineOn           = True
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1" 

  res@cnLevelSelectionMode =  "ManualLevels"
  res@cnFillPalette       = "MPL_jet"

do i=0,nfiles-1
  fobs = f->sp(i,:,:)
  if(i.eq.0) then
   fmean = fobs
   fmean = 0.
  end if
  fmean = fmean + fobs
end do
fmean = fmean/tofloat(nfiles)

  res@cnMinLevelValF       = -30.
  res@cnMaxLevelValF       = 0.
  res@cnLevelSpacingF      = 1.

time = todouble(f->valid_time)
copy_VarAtts(f->valid_time,time)
do i=0,nfiles-1
  fobs = f->sp(i,:,:)
  fobs = fobs - fmean
  fobs = fobs*0.01
  print("PS: "+min(fobs)+"  "+max(fobs))
  dates = tointeger(cd_calendar(time(i), 0))
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
  plot(0) = gsn_csm_contour_map(wks,fobs,res) 



  gsn_panel(wks,plot,(/1,1/),resP)                     ; now draw as one plot
end do
