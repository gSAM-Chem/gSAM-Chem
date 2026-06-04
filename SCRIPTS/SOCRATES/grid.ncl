load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; SAM datset
;filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_nonudge_4608_0000057120"
filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_nudge24h_SAM1MOM_gmg6_4608_0000057600"
;-----------------------------------------------------------------------------------------
; ERA5 datset

;date = 2018017023
date = 2018022100


;-----------------------------------------------------------------------------------------
; prepare initial settings for plots:

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = False                 ; turn on color fill
  res@cnLinesOn       = True
  res@cnFillPalette        = "BlueYellowRed"             ; choose colormap
  res@tmYRMode             = "Automatic"
  res@gsnYAxisIrregular2Linear = True
;  res@trYReverse        = True

;-----------------------------------------------------------------------------------------

; ------ U
f=addfile(filenameroot+"_U.nc","r")
fld = f->U(0,:,::36,0)
print(dimsizes(fld))
z = f->z
z = z/1000.
do k=0,dimsizes(z)-1
 fld(k,:) = z(k)
end do
fld&z = z
res@cnFillMode = "RasterFill"
res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
res@cnLevels = z
plot  = gsn_csm_contour(wks, fld({0.:3.},:), res )   ; place holder
