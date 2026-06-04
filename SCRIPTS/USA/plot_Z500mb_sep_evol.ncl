; modle output:

;files = systemfunc("ls -1 ../../OUT_2D/USA_6144x3072x72_1.25km_10s_APR10_2019-04-1[1234]-00-00-00_00*.2D_atm.nc")
files = systemfunc("ls -1 ../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_2019-04-1[1234]-00-00-00_00*.2D_atm.nc")
print(""+files)

;files = (/"../../OUT_2D/USA_SM_768x768x72_3.75km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc"/)

;files = (/"../../OUT_2D/USA_1536x1536x72_2.5km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc", \
;          "../../OUT_2D/USA_MID_1536x1536x72_2.5km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc", \
;          "../../OUT_2D/USA_MID_1536x1536x72_2.5km_10s_APR10_NUDGEUV12_2019-04-14-00-00-00_0000034560.2D_atm.nc" /)

;files = (/"../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR13_2019-04-14-00-00-00_0000008640.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR12_2019-04-14-00-00-00_0000017280.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR11_2019-04-14-00-00-00_0000025920.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_NUDGEUV12_2019-04-14-00-00-00_0000034560.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_NUDGEUV24_2019-04-14-00-00-00_0000034560.2D_atm.nc"/)

; ERA5 dataset
filesobs = systemfunc(" ls -1 /glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/201904/e5.oper.an.pl.128_129_z.ll025sc.2019041[1234]00_2019041[1234]23.nc")
fobs = addfiles(filesobs,"r")
print(""+filesobs)

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

 fldobs1 = fobs[i]->Z(0,{500},{latmin:latmax},{lonmin:lonmax}) ; 500-mb geopotential
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
;   fldobs = fldobs - fldobs_mean
 delete(fldobs1)
 print("fldobs_mean= "+fldobs_mean )
 print("fldobs perturb: "+min(fldobs)+" "+max(fldobs))
 base = gsn_csm_map(wks, res)
 res@cnLineColor           = color(0)
 res@cnMinLevelValF       = 5200.
 res@cnMaxLevelValF       = 5900.
 res@cnLevelSpacingF      = 25.
 res@cnLineLabelBackgroundColor = -1
 fldobs@long_name = ""
 fldobs@units = ""
 plot(0) = gsn_csm_contour(wks, fldobs, res)
 overlay(base, plot(0))

 res@cnLineLabelsOn = False
 res@cnLineColor           = color(1)
 ww = cos(lat*0.017453292519943)
 fld_lo = area_hi2lores_Wrap(lon, lat, fld, False, ww, lonobs, latobs, False)
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

