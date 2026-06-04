f = addfile("../../OUT_2D/CA_3072x3072x72_1km_10s_KK_2021-10-27-00-00-00_0000047520.2D_atm.nc","r")

wks = gsn_open_wks("ps","ps")


lat = f->lat
lon = f->lon
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

;res@lbLabelBarOn        = False           ; turn off individual cb's

;************************************************
;  res@gsnPanelLabelBar                 = True                    ; add common colorbar
;  resP@gsnPanelFigureStrings            = (/"gSAM","ERA5"/)           ; add strings to panel
;  res@gsnPanelFigureStringsFontHeightF = 0.015                   

  plot    = new(2,graphic)  

;  res@cnLevelSelectionMode =  "ManualLevels"   
;  res@cnMinLevelValF       = 3.
;  res@cnMaxLevelValF       = 60.
;  res@cnLevelSpacingF      = 3. 

  res@gsnCenterString  = "Accum. Precip."
  res@cnFillPalette       = "prcp_2"

  fld = f->Precac(0,:,:)
  time = f->time
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
;  res@gsnRightString = year + "-" + month + "-" + day + "-" + hour
  res@gsnLeftString = "gSAM"
  plot = gsn_csm_contour_map(wks,fld,res) 

