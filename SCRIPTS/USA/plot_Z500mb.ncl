; modle output:
files = (/"../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR13_2019-04-14-00-00-00_0000008640.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR12_2019-04-14-00-00-00_0000017280.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR11_2019-04-14-00-00-00_0000025920.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc"/)
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR11_NUDGEUV12_2019-04-14-00-00-00_0000025920.2D_atm.nc"/)

; ERA5 dataset
fobs = addfile("/glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/201904/e5.oper.an.pl.128_129_z.ll025sc.2019041400_2019041423.nc","r")

f = addfiles(files,"r")
nfiles = dimsizes(files)

wks = gsn_open_wks("ps","ps")

res = True
res@gsnDraw             = False
res@gsnFrame            = False
res@cnFillOn            = False
res@cnLinesOn           = True
res@cnLineThicknessF    = 4.0
res@cnInfoLabelOn       = False
res@lbLabelBarOn        = False
res@cnLevelSelectionMode = "ManualLevels"


color = (/"black","red","blue","green","magenta"/)

;---Base map (empty), so we can overlay multiple contours of different grids
res@gsnMaximize          = True
res@mpOutlineOn          = True
res@mpGeophysicalLineThicknessF = 1.5

plot    = new(nfiles+1,graphic)

; plot observations first:


do i=0,nfiles-1

 fld = f[i]->PHI500(0,:,:)
 print("fld: "+min(fld)+" "+max(fld))
 
 lat = f[i]->lat
 lon = f[i]->lon
 if(i.eq.0) then
  latmin = min(lat)
  latmax = max(lat)
  lonmin = min(lon)
  lonmax = max(lon)
  res@mpMinLatF   = latmin
  res@mpMaxLatF   = latmax
  res@mpMinLonF   = lonmin
  res@mpMaxLonF   = lonmax
 end if
 ; find domain average:
 fld1 = dim_avg_Wrap(fld)
 w = cos(lat * 0.0174533)
 w2D = conform(fld1, w, 0)
 fld_mean = dim_avg_wgt_n(fld1,w2D,0,0)
 fld = fld - fld_mean
 print("fld perturb: "+min(fld)+" "+max(fld))

; plot pbservation:
 if(i.eq.0) then
   fldobs1 = fobs->Z(0,{500},{latmin:latmax},{lonmin:lonmax}) ; 500-mb geopotential
   fldobs = fldobs1(::-1,:)
   delete(fldobs1)
   fldobs = fldobs/9.797  ; from geopotential to height
   print("fldobs: "+min(fldobs)+" "+max(fldobs))
   ; find domain average:
   latobs = fldobs&latitude
   lonobs = fldobs&longitude
   fldobs1 = dim_avg_Wrap(fldobs)
   wobs = cos(latobs * 0.0174533)
   w2Dobs = conform(fldobs1, wobs, 0)
   fldobs_mean = dim_avg_wgt_n(fldobs1,w2Dobs,0,0)
   fldobs = fldobs - fldobs_mean
   print("fldobs perturb: "+min(fldobs)+" "+max(fldobs))
   res = True
   base = gsn_csm_map(wks, res)
   res@cnLineColor           = color(0)
   res@cnMinLevelValF       = -400.
   res@cnMaxLevelValF       = 400.
   res@cnLevelSpacingF      = 40.
   res@cnLineLabelBackgroundColor = -1
   plot(0) = gsn_csm_contour(wks, fldobs, res)
   overlay(base, plot(0))
 end if

 res@cnLineLabelsOn = False
 res@cnLineColor           = color(i+1)
 ww = cos(lat*0.017453292519943)
 fld_lo = area_hi2lores_Wrap(lon, lat, fld, False, ww, lonobs, latobs, False)
 plot(i+1) = gsn_csm_contour(wks, fld_lo, res)
 overlay(base, plot(i+1))

; cleanup if next plot uses a different grid:
 delete(fld)
 delete(fld1)
 delete(lat)
 delete(lon)
 delete(w)
 delete(w2D)
 

end do

draw(base)
frame(wks)
