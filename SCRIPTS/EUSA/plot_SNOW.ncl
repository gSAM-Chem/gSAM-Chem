files = systemfunc("ls -1 ../../OUT_2D/EUSA_384x384x72_9km_10s_SNOWOUT_NOGR_PRMIN_GFS_20*.2D_atm.nc")
;files = systemfunc("ls -1 ../../OUT_2D/EUSA_768x768x72_4.5km_10s_SNOWOUT_NOGR_PRMIN_20*.2D_atm.nc")
nfiles = dimsizes(files)
print("nfiles = "+nfiles)

g = addfile("/glade/work/marat/ERA5_Hourly/NOR/ERA5_SFC_JAN3-5_2018_NOR_hourly.nc","r")
nstart = 0 
nend = 24*6+nstart-1 


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

;lonmin = -88. + 360.
;lonmax = -65. + 360.
;latmin = 30.
;latmax = 48.

res@gsnAddCyclic          = False
res@mpLimitMode            = "LatLon"
res@mpMinLatF              = latmin
res@mpMaxLatF              = latmax
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = False           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
  resP@gsnPanelLabelBar                 = True                    ; add common colorbar
;  resP@gsnPanelFigureStrings            = (/"gSAM","ERA5"/)           ; add strings to panel
  resP@gsnPanelFigureStringsFontHeightF = 0.015                   
  resP@gsnMaximize = True 

  plot    = new(1,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnMinLevelValF       = 0.
  res@cnMaxLevelValF       = 25.
  res@cnLevelSpacingF      = 1. 

  resP@gsnPanelMainString  = "Accumulated LWE Snow, mm"
  res@cnFillPalette       = "BkBlAqGrYeOrReViWh200"
  res@mpGeophysicalLineColor       = "white"
  res@mpUSStateLineColor       = "white"


  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"


  fld = f[nend]->Precsac(0,:,:)
  fld = fld - f[nstart]->Precsac(0,:,:)
;  fld = fld * 1000.0
  print("fld="+min(fld)+"  "+max(fld))
  time = f[nend]->time
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
  lat = f[0]->lat
  w = cos(lat * 0.0174533)
  w2D = conform(fld1, w, 0)
  fld_mean = dim_avg_wgt_n(fld1,w2D,0,0)
  print("mean="+fld_mean)
  res@gsnLeftString = "gSAM "+fld_mean
  plot(0) = gsn_csm_contour_map(wks,fld,res) 
  delete(fld)
  delete(fld1)
  delete(fld_mean)


  gsn_panel(wks,plot,(/1,1/),resP)                     ; now draw as one plot
