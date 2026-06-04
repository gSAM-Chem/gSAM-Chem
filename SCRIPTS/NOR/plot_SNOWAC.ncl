files = systemfunc("ls -1 ../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_2018*.2D_atm.nc")
nfiles = dimsizes(files)

dt = 1./24.  ; time step between data sample in days

g = addfile("/glade/work/marat/ERA5_Hourly/NOR/ERA5_SFC_JAN3-5_2018_NOR_hourly.nc","r")
nstart = 12 
nend = 24*2+nstart 


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

;lonmin = min(lon)
;lonmax = max(lon)
;latmin = min(lat)
;latmax = max(lat)

lonmin = -88. + 360.
lonmax = -65. + 360.
latmin = 30.
latmax = 48.

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

  plot    = new(2,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnMinLevelValF       = 0.
  res@cnMaxLevelValF       = 50.
  res@cnLevelSpacingF      = 1. 

  resP@gsnPanelMainString  = "LWE Accumulated Snow, mm"

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"

  colors = (/ (/255,255,255/), (/232,239,244/), \
              (/196,214,228/), (/132,173,209/), (/82,132,181/), \
              (/41,87,150/), (/19,54,139/), (/253,248,165/), \
              (/241,196,83/), (/230,144,70/), (/189,74,48/), \
              (/134,52,39/), (/89,29,19/), (/35,13,48/), (/65,38,114/)/) /255.
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels             = (/0.1,1.,2.,3.,4.,6.,8.,12.,18.,24.,30.,36.,48./)
  res@cnFillPalette        = colors


do i=nstart,nend
  tmp = f[i]->TSFC(0,:,:)
  land = f[i]->LAND(0,:,:)
  fld = f[i]->Preca(0,:,:)
  if(i.eq.nstart) then
   fldac = fld
   fldac = 0.
  end if
  fld = where(tmp.lt.273..and.land.ne.0.,fld,0.)
  fldac = fldac + fld 
  print(i+"  fldac="+min(fldac)+"  "+max(fldac))
end do
coef = dt * 15./25.4
fldac = fldac * coef
print("  fldac="+min(fldac)+"  "+max(fldac))
plot(0) = gsn_csm_contour_map(wks,fldac,res) 
delete(fld)
delete(fldac)

fld = g->sd(nend,:,:)
fld = fld - g->sd(nstart,:,:)
fld = fld * 1000.0 * 15./25.4
print("fldobs="+min(fld)+"  "+max(fld))
fld1 = dim_avg_Wrap(fld)
lat1 = g->latitude
w1 = cos(lat1 * 0.0174533)
w2D1 = conform(fld1, w1, 0)
fld_mean = dim_avg_wgt_n(fld1,w2D1,0,0)
print("mean_obs="+fld_mean)
res@gsnLeftString = "ERA5 "+fld_mean
plot(1) = gsn_csm_contour_map(wks,fld,res) 
delete(fld)
delete(fld1)
delete(fld_mean)

 gsn_panel(wks,plot,(/1,2/),resP)                     ; now draw as one plot
