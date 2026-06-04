; plot dayly average surface fluxes and compare them to ERA5

files = systemfunc("ls -1 ../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR13_2019-04-14-*.2D_atm.nc")
; ERA5 data:
g = addfile("/glade/work/marat/ERA5_Hourly/USA/data_stream-oper_stepType-instant.nc","r")

nfiles = dimsizes(files)

ntime_start_obs = 24*4
ntime_end_obs = ntime_start_obs+23

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

  plot    = new(2,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"

;==========================================================================
res@cnMinLevelValF       = 0.
res@cnMaxLevelValF       = 20.
res@cnLevelSpacingF      = 1. 
resP@gsnPanelMainString  = "10m Wind, m/s"
res@cnFillPalette       = "MPL_jet"


fld = f[0]->T2m(0,:,:)
fld = 0.
do i=0,nfiles-1
  fld1 = f[i]->U10m(0,:,:)
  fld2 = f[i]->V10m(0,:,:)
  fld = fld + sqrt(fld1^2 + fld2^2)
end do
fld = fld / tofloat(nfiles)
print("fld="+min(fld)+"  "+max(fld))

fldobs = g->u10(0,:,:)
fldobs = 0.
do i=ntime_start_obs,ntime_end_obs
  fldobs1 = g->u10(i,:,:)
  fldobs2 = g->v10(i,:,:)
  fldobs = fldobs + sqrt(fldobs1^2 + fldobs2^2)
end do
fldobs = fldobs / tofloat(ntime_end_obs-ntime_start_obs+1)
print("fld="+min(fldobs)+"  "+max(fldobs))

time = f[0]->time
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

gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot

;================================================================================================


res@cnMinLevelValF       = 260.
res@cnMaxLevelValF       = 310.
res@cnLevelSpacingF      = 1. 
resP@gsnPanelMainString  = "2m Temparature, K"

fld = 0.
do i=0,nfiles-1
  fld1 = f[i]->T2m(0,:,:)
  fld = fld + fld1
end do
fld = fld / tofloat(nfiles)
print("fld="+min(fld)+"  "+max(fld))

fldobs = 0.
do i=ntime_start_obs,ntime_end_obs
  fldobs1 = g->t2m(i,:,:)
  fldobs = fldobs + fldobs1
end do
fldobs = fldobs / tofloat(ntime_end_obs-ntime_start_obs+1)
print("fld="+min(fldobs)+"  "+max(fldobs))

time = f[0]->time
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

gsn_panel(wks,plot,(/2,1/),resP)                     ; now draw as one plot

