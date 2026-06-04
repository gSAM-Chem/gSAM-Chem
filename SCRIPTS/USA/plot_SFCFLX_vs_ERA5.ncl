; plot dayly average surface fluxes and compare them to ERA5

files = systemfunc("ls -1 ../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_NUDGEUV12_2019-04-*.2D_atm.nc")
; ERA5 data:
g = addfile("/glade/work/marat/ERA5_Hourly/USA/data_stream-oper_stepType-avg.nc","r")

nfiles = dimsizes(files)

ntime_start = 24*1+12
ntime_end = ntime_start+6
;ntime_end = ntime_start+23

ntime_start_obs = 24*1+12
ntime_end_obs = ntime_start_obs+6
;ntime_end_obs = ntime_start_obs+23

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

 res@tmXBTickSpacingF = 10              ; change maj lat tm spacing
 res@tmYLTickSpacingF = 10              ; change maj lon tm spacing


res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = True           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
  resP@gsnPanelLabelBar                 = False                    ; add common colorbar
;  resP@gsnPanelFigureStrings            = (/"gSAM","ERA5"/)           ; add strings to panel
  resP@gsnPanelFigureStringsFontHeightF = 0.015                   

  plot    = new(2,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnMinLevelValF       = 0.
  res@cnMaxLevelValF       = 2.
  res@cnLevelSpacingF      = 0.1 

  resP@gsnPanelMainString  = "Surface Stress, N/m2"
  res@cnFillPalette       = "MPL_jet"

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"


fld = f[0]->TAUX(0,:,:)
fld = 0.
do i=ntime_start,ntime_end
  fld1 = f[i]->TAUX(0,:,:)
  fld2 = f[i]->TAUY(0,:,:)
  fld = fld + sqrt(fld1^2 + fld2^2)
end do
fld = fld / tofloat(ntime_end-ntime_start+1)
print("fld="+min(fld)+"  "+max(fld))

fldobs = g->avg_iews(0,:,:)
fldobs = 0.
do i=ntime_start_obs,ntime_end_obs
  fldobs1 = g->avg_iews(i,:,:)
  fldobs2 = g->avg_inss(i,:,:)
  fldobs = fldobs + sqrt(fldobs1^2 + fldobs2^2)
end do
fldobs = fldobs / tofloat(ntime_end_obs-ntime_start_obs+1)
print("fldobs="+min(fldobs)+"  "+max(fldobs))

time = dim_avg_Wrap(f[ntime_start-1:ntime_end-1]->time)
dates = tointeger(cd_calendar(time, 0))
year = dates(0,0)
month = sprinti("%0.2d", dates(0,1))
day = sprinti("%0.2d", dates(0,2))
mins = sprinti("%0.2d", dates(0,4))
if(mins.eq.59) then
 dates(0,3) = dates(0,3)+1
end if
hour = sprinti("%0.2d", dates(0,3))
print(year + "-" + month + "-" + day)
res@gsnRightString = year + "-" + month + "-" + day 
res@gsnLeftString = "gSAM "
plot(0) = gsn_csm_contour_map(wks,fld,res) 
res@gsnLeftString = "ERA5 "
plot(1) = gsn_csm_contour_map(wks,fldobs,res) 
;wgt = cos(0.008726646*fld&lat)
;fldd = fldobs({latmin:latmax},{lonmin-360.:lonmax-360.})
;print("fldd="+min(fldd)+"  "+max(fldd))
;lono = fldd&longitude
;lono = lono+360.
;lato = fldd&latitude
;lato = lato
;fld_lowres = area_hi2lores_Wrap(lon,lat,fld,False,wgt,lono,lato,False)
;res@gsnLeftString = "gSAM (ERA5)"
;print("fld_lowres="+min(fld_lowres)+"  "+max(fld_lowres))
;plot(2) = gsn_csm_contour_map(wks,fld_lowres,res) 
;fld_lowres = fld_lowres - fldd
;res@gsnLeftString = "gSAM - ERA5"
;res@cnLevelSelectionMode =  "ManualLevels"   
;res@cnMinLevelValF       = -1.
;res@cnMaxLevelValF       = 1.
;res@cnLevelSpacingF      = 0.1 
;plot(3) = gsn_csm_contour_map(wks,fld_lowres,res) 

gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot

;================================================================================================

res@cnFillPalette       = "NCV_blu_red"

res@cnMinLevelValF       = -300.
res@cnMaxLevelValF       = 300.
res@cnLevelSpacingF      = 20. 
resP@gsnPanelMainString  = "Latent Heat Flux, W/m2"

fld = 0.
do i=ntime_start,ntime_end
  fld1 = f[i]->LHF(0,:,:)
  fld = fld + fld1
end do
fld = fld / tofloat(ntime_end-ntime_start+1)
print("fld="+min(fld)+"  "+max(fld))

fldobs = 0.
do i=ntime_start_obs,ntime_end_obs
  fldobs1 = g->avg_slhtf(i,:,:)
  fldobs = fldobs - fldobs1
end do
fldobs = fldobs / tofloat(ntime_end_obs-ntime_start_obs+1)
print("fld="+min(fldobs)+"  "+max(fldobs))

time = dim_avg_Wrap(f[ntime_start:ntime_end]->time)
dates = tointeger(cd_calendar(time, 0))
year = dates(0,0)
month = sprinti("%0.2d", dates(0,1))
day = sprinti("%0.2d", dates(0,2))
mins = sprinti("%0.2d", dates(0,4))
if(mins.eq.59) then
 dates(0,3) = dates(0,3)+1
end if
hour = sprinti("%0.2d", dates(0,3))
print(year + "-" + month + "-" + day)
res@gsnRightString = year + "-" + month + "-" + day
res@gsnLeftString = "gSAM "
plot(0) = gsn_csm_contour_map(wks,fld,res)
res@gsnLeftString = "ERA5 "
plot(1) = gsn_csm_contour_map(wks,fldobs,res)
lwf = fld
lwfobs = fldobs
;wgt = cos(0.008726646*fld&lat)
;fldd = fldobs({latmin:latmax},{lonmin-360.:lonmax-360.})
;print("fldd="+min(fldd)+"  "+max(fldd))
;lono = fldd&longitude
;lono = lono+360.
;lato = fldd&latitude
;lato = lato
;fld_lowres = area_hi2lores_Wrap(lon,lat,fld,False,wgt,lono,lato,False)
;res@gsnLeftString = "gSAM (ERA5)"
;print("fld_lowres="+min(fld_lowres)+"  "+max(fld_lowres))
;plot(2) = gsn_csm_contour_map(wks,fld_lowres,res)
;fld_lowres = fld_lowres - fldd
;res@gsnLeftString = "gSAM - ERA5"
;res@cnLevelSelectionMode =  "ManualLevels"
;res@cnMinLevelValF       = -100.
;res@cnMaxLevelValF       = 100.
;res@cnLevelSpacingF      = 10.
;plot(3) = gsn_csm_contour_map(wks,fld_lowres,res)


gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot

;================================================================================================

res@cnMinLevelValF       = -200.
res@cnMaxLevelValF       = 200.
res@cnLevelSpacingF      = 20.
resP@gsnPanelMainString  = "Sensible Heat Flux, W/m2"

fld = 0.
do i=ntime_start,ntime_end
  fld1 = f[i]->SHF(0,:,:)
  fld = fld + fld1
end do
fld = fld / tofloat(ntime_end-ntime_start+1)
print("fld="+min(fld)+"  "+max(fld))

fldobs = 0.
do i=ntime_start_obs,ntime_end_obs
  fldobs1 = g->avg_ishf(i,:,:)
  fldobs = fldobs - fldobs1
end do
fldobs = fldobs / tofloat(ntime_end_obs-ntime_start_obs+1)
print("fld="+min(fldobs)+"  "+max(fldobs))

time = dim_avg_Wrap(f[ntime_start:ntime_end]->time)
dates = tointeger(cd_calendar(time, 0))
year = dates(0,0)
month = sprinti("%0.2d", dates(0,1))
day = sprinti("%0.2d", dates(0,2))
mins = sprinti("%0.2d", dates(0,4))
if(mins.eq.59) then
 dates(0,3) = dates(0,3)+1
end if
hour = sprinti("%0.2d", dates(0,3))
print(year + "-" + month + "-" + day)
res@gsnRightString = year + "-" + month + "-" + day
res@gsnLeftString = "gSAM "
plot(0) = gsn_csm_contour_map(wks,fld,res)
res@gsnLeftString = "ERA5 "
plot(1) = gsn_csm_contour_map(wks,fldobs,res)
;wgt = cos(0.008726646*fld&lat)
;fldd = fldobs({latmin:latmax},{lonmin-360.:lonmax-360.})
;print("fldd="+min(fldd)+"  "+max(fldd))
;lono = fldd&longitude
;lono = lono+360.
;lato = fldd&latitude
;lato = lato
;fld_lowres = area_hi2lores_Wrap(lon,lat,fld,False,wgt,lono,lato,False)
;res@gsnLeftString = "gSAM (ERA5)"
;print("fld_lowres="+min(fld_lowres)+"  "+max(fld_lowres))
;plot(2) = gsn_csm_contour_map(wks,fld_lowres,res)
;fld_lowres = fld_lowres - fldd
;res@gsnLeftString = "gSAM - ERA5"
;res@cnLevelSelectionMode =  "ManualLevels"
;res@cnMinLevelValF       = -100.
;res@cnMaxLevelValF       = 100.
;res@cnLevelSpacingF      = 10.
;plot(3) = gsn_csm_contour_map(wks,fld_lowres,res)

gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot

;============================================================================
res@cnMinLevelValF       = -400.
res@cnMaxLevelValF       = 400.
res@cnLevelSpacingF      = 40.
resP@gsnPanelMainString  = "Surface Enthalpy Flux, W/m2"

fld = 0.
do i=ntime_start,ntime_end
  fld1 = f[i]->SHF(0,:,:)
  fld = fld + fld1
end do
fld = fld / tofloat(ntime_end-ntime_start+1)
print("fld="+min(fld)+"  "+max(fld))

fldobs = 0.
do i=ntime_start_obs,ntime_end_obs
  fldobs1 = g->avg_ishf(i,:,:)
  fldobs = fldobs - fldobs1
end do
fldobs = fldobs / tofloat(ntime_end_obs-ntime_start_obs+1)
print("fld="+min(fldobs)+"  "+max(fldobs))

time = dim_avg_Wrap(f[ntime_start:ntime_end]->time)
dates = tointeger(cd_calendar(time, 0))
year = dates(0,0)         
month = sprinti("%0.2d", dates(0,1))
day = sprinti("%0.2d", dates(0,2))            
mins = sprinti("%0.2d", dates(0,4))
if(mins.eq.59) then
 dates(0,3) = dates(0,3)+1
end if
hour = sprinti("%0.2d", dates(0,3))
print(year + "-" + month + "-" + day)
res@gsnRightString = year + "-" + month + "-" + day
res@gsnLeftString = "gSAM "
fld = fld + lwf
plot(0) = gsn_csm_contour_map(wks,fld,res)
res@gsnLeftString = "ERA5 "
fldobs = fldobs + lwfobs
plot(1) = gsn_csm_contour_map(wks,fldobs,res)
;wgt = cos(0.008726646*fld&lat)
;fldd = fldobs({latmin:latmax},{lonmin-360.:lonmax-360.})
;print("fldd="+min(fldd)+"  "+max(fldd))
;lono = fldd&longitude
;lono = lono+360.
;lato = fldd&latitude
;lato = lato
;fld_lowres = area_hi2lores_Wrap(lon,lat,fld,False,wgt,lono,lato,False)
;res@gsnLeftString = "gSAM (ERA5)"
;print("fld_lowres="+min(fld_lowres)+"  "+max(fld_lowres))
;plot(2) = gsn_csm_contour_map(wks,fld_lowres,res)
;fld_lowres = fld_lowres - fldd
;res@gsnLeftString = "gSAM - ERA5"
;res@cnLevelSelectionMode =  "ManualLevels"
;res@cnMinLevelValF       = -100.
;res@cnMaxLevelValF       = 100.
;res@cnLevelSpacingF      = 10.
;plot(3) = gsn_csm_contour_map(wks,fld_lowres,res)

gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot

