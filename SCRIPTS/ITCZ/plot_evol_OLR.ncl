files = systemfunc("ls -1 ../../OUT_2D/ITCZ_6144x3072x72_3.75km_10s_20*.2D_atm.nc")
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

lonmin = -64.+360.
lonmax = -8.+360.
latmin = -4.
latmax = 24.
;lonmin = min(lon)
;lonmax = max(lon)
;latmin = min(lat)
;latmax = max(lat)
;
res@gsnAddCyclic          = False
res@mpLimitMode            = "LatLon"
res@mpMinLatF              = latmin 
res@mpMaxLatF              = latmax 
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax 

res@lbLabelBarOn        = True           ; turn off individual cb's

;************************************************

  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnMinLevelValF       = 100.
  res@cnMaxLevelValF       = 300.
  res@cnLevelSpacingF      = 1.0 
  res@lbLabelStride = 40
  res@lbBoxLinesOn         = False 
  res@lbOrientation        = "Vertical"
  res@lbTitlePosition      = "Right" 

  res@tiMainString  = "Cloud Water Path, mm"
  res@cnFillPalette       = "MPL_Greys"

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"
  res@mpGeophysicalLineColor       = "white"


do i=0,nfiles-1
  fld = f[i]->LWNT(0,:,:)
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
  res@gsnRightString = year + "-" + month + "-" + day + "-" + hour + "-"+mins
  fld1 = dim_avg_Wrap(fld)
  if(i.eq.0) then
    lat = f[i]->lat
    w = cos(lat * 0.0174533)
    w2D = conform(fld1, w, 0)
  end if
  fld_mean = dim_avg_wgt_n(fld1,w2D,0,0)
  print("mean="+fld_mean)
  res@gsnLeftString = "gSAM "+fld_mean
  plot = gsn_csm_contour_map(wks,fld,res) 
end do

