load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; SAM datset
filenameroot = "../../OUT_3D/SOCRATES_9216x4608x74_7.5s_4km_nudge24h_SAM1MOM_gmg6_4608_0000011520"
;-----------------------------------------------------------------------------------------
; ERA5 datset

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
  res                      = True                 ; plot mods desired
  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
;  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnLinesOn       = False
  res@cnFillPalette        = "BlueYellowRed"             ; choose colormap
  res@tmYRMode             = "Automatic"
  res@gsnYAxisIrregular2Linear = True
;  res@trYReverse        = True
  res@vpWidthF    = 0.8      ; Make plots wider than 
  res@vpHeightF   = 0.2      ; they are high.
;  res@tmYLOn               = False   ; Turn off left tickmarks
  res@tmYROn               = False    ; Turn on right tickmarks
  res@tmXTOn               = False   ; Turn off top tickmarks
  res@tmXBLabelsOn         = False   ; Turn off left labels
  res@gsnStringFontHeightF = 0.02
 


;controlling label bar
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = True         ; labelbar box lines.
  res@lbTitlePosition     = "Bottom"
  res@lbOrientation        = "Vertical"
  res@lbTitleFontHeightF = 0.02
  res@lbLabelFontHeightF = 0.015
  res@pmLabelBarHeightF = 0.2
  res@lbBoxMinorExtentF = 0.2
  res@lbBoxSeparatorLinesOn=False
  res@gsnStringFontHeightF = 0.02

  res@gsnDraw         = False               ; Do not draw plot
  res@gsnFrame        = False               ; Do not advance frame

plot = new(4,graphic)

;-----------------------------------------------------------------------------------------

zmax = 6000.
lon1 = 115.
lon2 = 160.
lat0 = -53.

f=addfile(filenameroot+"_TABS.nc","r")
fld = f->TABS(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
plot(0)  = gsn_csm_contour(wks, fld, res )   ; place holder
f=addfile(filenameroot+"_QV.nc","r")
fld = f->QV(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
plot(1)  = gsn_csm_contour(wks, fld, res )   ; place holder
f=addfile(filenameroot+"_U.nc","r")
fld = f->U(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
plot(2)  = gsn_csm_contour(wks, fld, res )   ; place holder
f=addfile(filenameroot+"_V.nc","r")
fld = f->V(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
plot(3)  = gsn_csm_contour(wks, fld, res )   ; place holder
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind"
gsn_panel(wks,plot,(/4,1/),resp)           ; draw: 2-across, 2-down

f=addfile(filenameroot+"_QC.nc","r")
fld = f->QC(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
plot(0)  = gsn_csm_contour(wks, fld, res )   ; place holder
f=addfile(filenameroot+"_QI.nc","r")
fld = f->QI(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
plot(1)  = gsn_csm_contour(wks, fld, res )   ; place holder
f=addfile(filenameroot+"_QP.nc","r")
fld = f->QP(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
plot(2)  = gsn_csm_contour(wks, fld, res )   ; place holder
f=addfile(filenameroot+"_W.nc","r")
fld = f->W(0,{0.:zmax},{lat0},{lon1:lon2})
print("t: "+min(fld)+"  "+max(fld))
;res@cnLevelSpacingF      = 5.                  ; contour spacing
;res@cnMinLevelValF       = -50.                ; min level
;res@cnMaxLevelValF       = 50.                 ; max level
res@cnFillPalette        = "BlueWhiteOrangeRed"             ; choose colormap
plot(3)  = gsn_csm_contour(wks, fld, res )   ; place holder
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Zonal Wind"
gsn_panel(wks,plot,(/4,1/),resp)           ; draw: 2-across, 2-down

