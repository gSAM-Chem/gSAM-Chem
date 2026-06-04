files = systemfunc("ls -1 ../../OUT_2D/LI_3456x3456x83_100m_3s_Sep03_SPNUDGE013610_20*.2D_atm.nc")
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
res@lbLabelBarOn        = True           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
  resP@gsnPanelLabelBar                 = False                    ; add common colorbar
;  resP@gsnPanelFigureStrings            = (/"gSAM","ERA5"/)           ; add strings to panel
;  resP@gsnPanelFigureStringsFontHeightF = 0.015                   

  plot    = new(2,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   

  resP@gsnPanelMainString  = "  "
  res@cnFillPalette       = "MPL_jet"

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"


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
  res@gsnRightString = year + "-" + month + "-" + day + "-" + hour
  fld1 = dim_avg_Wrap(fld)
  if(i.eq.0) then
    lat = f[i]->lat
    w = cos(lat * 0.0174533)
    w2D = conform(fld1, w, 0)
  end if
  fld_mean = dim_avg_wgt_n(fld1,w2D,0,0)
  print("mean="+fld_mean)
  res@gsnLeftString = "LWNT "+fld_mean
  res@cnMinLevelValF       = 260.
  res@cnMaxLevelValF       = 360.
  res@cnLevelSpacingF      = 5. 
  print("LWNT: "+min(fld)+"  "+max(fld))
  plot(0) = gsn_csm_contour_map(wks,fld,res) 

  fobs = f[i]->SWNT(0,:,:)
  fobs2 = f[i]->SOLIN(0,:,:)
  fobs = fobs2-fobs
  fobs = fobs/(fobs2+1.)
  fld1 = dim_avg_Wrap(fobs)
  if(i.eq.0) then
    lat = f[i]->lat
    w = cos(lat * 0.0174533)
    w2D = conform(fld1, w, 0)
  end if
  fld_mean = dim_avg_wgt_n(fld1,w2D,0,0)
  print("mean_obs="+fld_mean)
  res@gsnLeftString = "SWNT "+fld_mean
  res@cnMinLevelValF       = 0.
  res@cnMaxLevelValF       = 0.7
  res@cnLevelSpacingF      = 0.05 
  print("ALB: "+min(fobs)+"  "+max(fobs))
  plot(1) = gsn_csm_contour_map(wks,fobs,res) 

  gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot
end do
