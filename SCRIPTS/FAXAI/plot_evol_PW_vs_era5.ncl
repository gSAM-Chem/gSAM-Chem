files = systemfunc("ls -1 ../../OUT_2D/FAXAI_1536x1536x72_3km_10s_NONUDGE_2019*.2D_atm.nc")
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
;  resP@gsnPanelFigureStrings            = (/"gSAM","ERA5"/)           ; add strings to panel
  resP@gsnPanelFigureStringsFontHeightF = 0.015                   

  plot    = new(2,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnMinLevelValF       = 3.
  res@cnMaxLevelValF       = 70.
  res@cnLevelSpacingF      = 2. 

  resP@gsnPanelMainString  = "Precipitable Water, mm"
  res@cnFillPalette       = "MPL_jet"

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"


do i=0,nfiles-1
  fld = f[i]->PW(0,:,:)
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
  res@gsnLeftString = "gSAM "+fld_mean
  plot(0) = gsn_csm_contour_map(wks,fld,res) 

  fobs = f[i]->PWOBS(0,:,:)
  fld1 = dim_avg_Wrap(fobs)
  if(i.eq.0) then
    lat = f[i]->lat
    w = cos(lat * 0.0174533)
    w2D = conform(fld1, w, 0)
  end if
  fld_mean = dim_avg_wgt_n(fld1,w2D,0,0)
  print("mean_obs="+fld_mean)
  res@gsnLeftString = "ERA5 "+fld_mean
  plot(1) = gsn_csm_contour_map(wks,fobs,res) 

  gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot
end do
