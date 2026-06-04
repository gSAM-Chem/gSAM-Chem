
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;----------------------------------------------------------------------
begin

filename = "../../OUT_2DZ/DYAMOND_9216x4608x74_7.5s_4km_4608_10-40d.com2DZ_1.nc"

  fin     = addfile (filename,"r")

  wks   = gsn_open_wks ("ps", "ps" )          ; send graphics to PNG file


  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnLinesOn       = False
  res@cnFillPalette        = "BlueYellowRed"             ; choose colormap
  res@tmYRMode             = "Automatic"
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse        = True

;controlling label bar
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = True         ; labelbar box lines.
  res@lbTitlePosition     = "Bottom"
  res@lbTitleFontHeightF = 0.02
  res@lbLabelFontHeightF = 0.015
  res@pmLabelBarHeightF = 0.1
  res@lbBoxMinorExtentF = 0.2
  res@lbBoxSeparatorLinesOn=False
  res@gsnStringFontHeightF = 0.02


  print(""+filename)


  fld     = fin->U(0,:,:)
  lev = fin->p
  lat = fin->lat
  fld!0 = "lev"
  fld&lev = lev
  fld!1 = "lat"
  fld&lat = lat
  res@cnLevelSpacingF      = 5.0                  ; contour spacing
  res@cnMinLevelValF       =-50.                 ; min level
  res@cnMaxLevelValF       = 50.                 ; max level
  fld@long_name="U"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  delete(fld)

  fld     = fin->V(0,:,:)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.5                  ; contour spacing
  res@cnMinLevelValF       = -5.5                 ; min level
  res@cnMaxLevelValF       =  5.5                 ; max level
  fld@long_name="V"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  delete(fld)

  fld     = fin->W(0,:,:)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.01                  ; contour spacing
  res@cnMinLevelValF       = -0.1                 ; min level
  res@cnMaxLevelValF       =  0.1                 ; max level
  fld@long_name="W"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  delete(fld)

  fld     = fin->QV(0,:,:)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 1.                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 18.                 ; max level
  fld@long_name="QV"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  delete(fld)

  fld     = fin->QRAD(0,:,:)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.1                  ; contour spacing
  res@cnMinLevelValF       = -2.0                 ; min level
  res@cnMaxLevelValF       =  2.0                 ; max level
  fld@long_name="QRAD"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  delete(fld)

  fld     = fin->QN(0,:,:)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.01                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 0.1                 ; max level
  fld@long_name="QN"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  delete(fld)

  fld     = fin->TABS(0,:,:)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 5.                  ; contour spacing
  res@cnMinLevelValF       = 190.                 ; min level
  res@cnMaxLevelValF       = 300.                 ; max level
  fld@long_name="TABS"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  delete(fld)



end
