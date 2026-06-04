; modle output:

;files = (/"../../OUT_2D/USA_SM_768x768x72_3.75km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc" /)

;files = (/"../../OUT_2D/USA_MID_1536x1536x72_2.5km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc", \
;         "../../OUT_2D/USA_MID_1536x1536x72_2.5km_10s_APR10_NUDGEUV12_2019-04-14-00-00-00_0000034560.2D_atm.nc"/)
;files = (/"../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR13_2019-04-14-00-00-00_0000008640.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR12_2019-04-14-00-00-00_0000017280.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR11_2019-04-14-00-00-00_0000025920.2D_atm.nc", \
;          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc"/)
files = (/"../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR13_2019-04-14-00-00-00_0000008640.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR12_2019-04-14-00-00-00_0000017280.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR11_2019-04-14-00-00-00_0000025920.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_2019-04-14-00-00-00_0000034560.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_NUDGEUV12_2019-04-14-00-00-00_0000034560.2D_atm.nc", \
          "../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR10_NUDGEUV24_2019-04-14-00-00-00_0000034560.2D_atm.nc"/)


; ERA5 dataset
fobs = addfile("/glade/work/marat/Radar_Reflectivity/uscomp_n0q_20190414_0000.nc","r")

f = addfiles(files,"r")
nfiles = dimsizes(files)

wks = gsn_open_wks("ps","ps")

res = True
res@cnFillOn            = True
res@cnLinesOn           = False
res@cnInfoLabelOn       = False
res@lbLabelBarOn        = True
res@cnLevelSelectionMode = "ManualLevels"
res@cnFillMode            = "RasterFill"

res@mpShapeMode = "FreeAspect"
;  res@vpXF                 = 0.1   ; Set X,Y location of upper left
;  res@vpYF                 = 0.98   ; of plot
  res@vpWidthF             = 0.80   ; Set width
  res@vpHeightF             = 0.53   ; Set height

 res@tmXBTickSpacingF = 10              ; change maj lat tm spacing
 res@tmYLTickSpacingF = 10              ; change maj lon tm spacing



color = (/"black","red","blue","green"/)

res@mpGeophysicalLineThicknessF = 1.5

; plot observations first:

res@cnMinLevelValF       = 0.
res@cnMaxLevelValF       = 50.
res@cnLevelSpacingF      = 2.
res@gsnAddCyclic = False

res@mpOutlineBoundarySets = "AllBoundaries"
res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
res@mpDataSetName = "Earth..1"

  colors = (/ (/255,255,255/), (/77,230,233/), \
              (/43,150,235/), (/0,35,232/), (/91,247,66/), \
              (/68,186,49/), (/47,128,34/), \
              (/255,252,72/),(/224,182,54/),(/245,136,45/), \
              (/239,47,29/),(/187,57,26/),(/134,27,16/), \
              (/237,60,245/),(/133,59,190/)/) / 255.
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels             = (/5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,70.,75./)
  res@cnFillPalette        = colors

do i=0,nfiles-1

 fld = f[i]->ZdBZ(0,:,:)
 print("fld: "+min(fld)+" "+max(fld))
 
; plot pbservation:
 if(i.eq.0) then
   lon = f[i]->lon
   fldobs = fobs->composite_n0q(::3,::3) 
   fldobs = where(ismissing(fldobs),-30.,fldobs)
   print("fldobs: "+min(fldobs)+" "+max(fldobs))

   latmin = min(fldobs&lat)
   latmax = max(fldobs&lat)
   lonmin = min(fldobs&lon)
   lonmax = max(fldobs&lon)
   print("lon:"+" "+lonmin+" "+lonmax)
   print("lat:"+" "+latmin+" "+latmax)
   res@mpLimitMode            = "LatLon"
   res@mpMinLatF   = latmin
   res@mpMaxLatF   = latmax
   res@mpMinLonF   = lonmin
   res@mpMaxLonF   = lonmax
   plot = gsn_csm_contour_map(wks, fldobs, res)
   res@mpMinLonF   = lonmin+360.
   res@mpMaxLonF   = lonmax+360.

 end if
 fld@long_name = files(i)
 res@gsnLeftStringFontHeightF           = 0.01
 plot = gsn_csm_contour_map(wks, fld, res)

end do

