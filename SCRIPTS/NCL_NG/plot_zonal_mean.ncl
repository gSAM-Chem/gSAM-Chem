
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;----------------------------------------------------------------------
begin

filename = "../../OUT_3D/HS94_288x144x74_20s_logtv_PB_144_mean.nc"
;filename = "../../OUT_3D/HS94_256x128x74_20s_128_mean.nc"

  wks   = gsn_open_wks ("ps", "ps" )          ; send graphics to PNG file


  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = True                 ; turn on line labels
  res@cnLinesOn       = True                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "BlAqGrYeOrReVi200"             ; choose colormap
  res@tmYRMode             = "Automatic"
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse        = True
;  res@trYMinF = 50.
  res@gsnDraw         = False                   ; don't draw yet
  res@gsnFrame        = False                   ; don't advance frame yet

  print(""+filename)


  fin     = addfile (filename,"r")

  f3d     = fin->U(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  lat = fin->lat
  fld!0 = "lev"
  fld&lev = lev
  fld!1 = "lat"
  fld&lat = lat
  res@cnLevelSpacingF      = 4.0                  ; contour spacing
  res@cnMinLevelValF       = -32.                 ; min level
  res@cnMaxLevelValF       =  32.                 ; max level
  fld@long_name="U"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)

  f3d     = fin->V(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.5                  ; contour spacing
  res@cnMinLevelValF       = -2.                 ; min level
  res@cnMaxLevelValF       =  2.                 ; max level
  fld@long_name="V"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)

  f3d     = fin->W(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.001                  ; contour spacing
  res@cnMinLevelValF       = -0.01                 ; min level
  res@cnMaxLevelValF       =  0.01                 ; max level
  fld@long_name="W"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)

  f3d     = fin->QV(0,:,:,:)*1.e3
  fld = dim_avg(f3d)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.3                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 3.                 ; max level
  fld@long_name="QV"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)

;  f3d     = fin->TT(0,:,:,:)
;  fld = dim_avg(f3d)
;  lev = fin->p
;  fld!0 = "lev"
;  fld&lev = lev
;  res@cnLevelSpacingF      = 5.                  ; contour spacing
;  res@cnMinLevelValF       = 300.                 ; min level
;  res@cnMaxLevelValF       = 500.                 ; max level
;  fld@long_name="TT"
;  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

;  f3d     = fin->QRAD(0,:,:,:)
;  fld = dim_avg(f3d)
;  lev = fin->p
;  fld!0 = "lev"
;  fld&lev = lev
;  res@cnLevelSpacingF      = 0.1                  ; contour spacing
;  res@cnMinLevelValF       = -2.0                 ; min level
;  res@cnMaxLevelValF       =  2.0                 ; max level
;  fld@long_name="QRAD"
;  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

;  f3d     = fin->QN(0,:,:,:)
;  fld = dim_avg(f3d)
;  lev = fin->p
;  fld!0 = "lev"
;  fld&lev = lev
;  res@cnLevelSpacingF      = 0.01                  ; contour spacing
;  res@cnMinLevelValF       = 0.                 ; min level
;  res@cnMaxLevelValF       = 0.1                 ; max level
;  fld@long_name="QN"
;;  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

  f3d     = fin->TABS(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 5.                  ; contour spacing
  res@cnMinLevelValF       = 180.                 ; min level
  res@cnMaxLevelValF       = 310.                 ; max level
  fld@long_name="TABS"
  plot = ZeroLineContour(plot)
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)



end
