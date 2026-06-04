; modle output:
files = (/"../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_2018-01-03-12-00-00_0000004320.2D_atm.nc", \
          "../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_2018-01-04-00-00-00_0000008640.2D_atm.nc", \
          "../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_2018-01-04-12-00-00_0000012960.2D_atm.nc", \
          "../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_2018-01-05-00-00-00_0000017280.2D_atm.nc", \
          "../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_2018-01-05-12-00-00_0000021600.2D_atm.nc"/)

; ERA5 dataset
fobs = addfile("/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/201801/e5.oper.an.sfc.128_134_sp.ll025sc.2018010100_2018013123.nc","r")
nobs = (/12+2*24,0+3*24,12+3*24,0+4*24,12+4*24/)
nobs0 = 48

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

 fld = f[i]->PSFC(0,:,:)
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

; plot pbservation:
; if(i.eq.0) then
   fldobs1 = fobs->SP(nobs(i),{latmin:latmax},{lonmin:lonmax}) 
   fldobs1 = fldobs1 * 0.01
   fldobs = fldobs1(::-1,:)
   delete(fldobs1)
   print("fldobs: "+min(fldobs)+" "+max(fldobs))
   ; find domain average:
   latobs = fldobs&latitude
   lonobs = fldobs&longitude
   fldobs1 = dim_avg_Wrap(fldobs)
   if(i.eq.0) then
     fldobs0 = fobs->SP(nobs0,{latmin:latmax},{lonmin:lonmax}) ; 500-mb geopotential
     fldobs0 = fldobs0 *0.01  ; from geopotential to height
     fldobs00 = fldobs0(::-1,:)
     fldobs000 = dim_avg_Wrap(fldobs00)
     wobs = cos(latobs * 0.0174533)
     w2Dobs = conform(fldobs000, wobs, 0)
     fldobs_mean = dim_avg_wgt_n(fldobs000,w2Dobs,0,0)
   end if
   delete(fldobs1)
;   fldobs = fldobs - fldobs_mean
   print("fldobs_mean= "+fldobs_mean )
   print("fldobs perturb: "+min(fldobs)+" "+max(fldobs))
; end if
 base = gsn_csm_map(wks, res)
 res@cnLineColor           = color(0)
 res@cnMinLevelValF       = 950.
 res@cnMaxLevelValF       = 1020.
 res@cnLevelSpacingF      = 5.
 res@cnLineLabelBackgroundColor = -1
 fldobs@long_name = ""
 fldobs@units = ""
 plot(0) = gsn_csm_contour(wks, fldobs, res)
 overlay(base, plot(0))

 res@cnLineLabelsOn = False
 res@cnLineColor           = color(1)
 ww = cos(lat*0.017453292519943)
 fld_lo = area_hi2lores_Wrap(lon, lat, fld, False, ww, lonobs(::2), latobs(::2), False)
 fld_lo = fld_lo + fldobs_mean
 print("fld_lo: "+min(fld_lo)+" "+max(fld_lo))
 fld_lo@long_name = files(i)
 fld_lo@units = "m"
 plot(1) = gsn_csm_contour(wks, fld_lo, res)
 overlay(base, plot(1))

; cleanup if next plot uses a different grid:
 delete(fld)
 delete(fld1)
 delete(lat)
 delete(lon)
 delete(w)
 delete(w2D)
 
draw(base)
frame(wks)

end do

