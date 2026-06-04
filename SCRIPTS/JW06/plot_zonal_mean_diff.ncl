
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;----------------------------------------------------------------------
begin

filename1 = "../../OUT_3D/JW06_512x256x30_20s_r8_nadv2_nosgs_256_0000000000.atm.3D.nc"
filename2 = "../../OUT_3D/JW06_512x256x30_20s_r8_nadv2_nosgs_256_0000129600.atm.3D.nc"

  wks   = gsn_open_wks ("ps", "ps" )          ; send graphics to PNG file


  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "ViBlGrWhYeOrRe"             ; choose colormap
  res@tmYRMode             = "Automatic"
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse        = True
;  res@trYMinF = 50.
  res@gsnDraw         = False                   ; don't draw yet
  res@gsnFrame        = False                   ; don't advance frame yet

  print(""+filename1)
  print(""+filename2)


  fin1     = addfile (filename1,"r")
  fin2     = addfile (filename2,"r")

  nz1 = dimsizes(fin1->p)-4

  f3d     = fin2->U(0,:nz1,:,:)
  f3d     = f3d - fin1->U(0,:nz1,:,:)
  fld = dim_avg(f3d)
  lev = fin1->p(:nz1)
  lat = fin1->lat
  fld!0 = "lev"
  fld&lev = lev
  fld!1 = "lat"
  fld&lat = lat
  res@cnLevelSpacingF      = 0.01                  ; contour spacing
  res@cnMinLevelValF       = -0.2                 ; min level
  res@cnMaxLevelValF       = 0.2                 ; max level
  fld@long_name="U"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)

  f3d     = fin2->V(0,:nz1,:,:)
  fld = dim_avg(f3d)
  lev = fin1->p(:nz1)
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.01                  ; contour spacing
  res@cnMinLevelValF       = -0.2                 ; min level
  res@cnMaxLevelValF       =  0.2                 ; max level
  fld@long_name="V"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)

  f3d     = fin2->W(0,:nz1,:,:)
  fld = dim_avg(f3d)
  lev = fin1->p(:nz1)
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

  f3d     = fin1->QV(0,:nz1,:,:)*1.e3
  fld = dim_avg(f3d)
  lev = fin1->p(:nz1)
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.03                  ; contour spacing
  res@cnMinLevelValF       = 0.3                ; min level
  res@cnMaxLevelValF       = 0.3                 ; max level
  fld@long_name="QV"
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)

  f3d     = fin2->TABS(0,:nz1,:,:)
  f3d     = f3d - fin1->TABS(0,:nz1,:,:)
  fld = dim_avg(f3d)
  lev = fin1->p(:nz1)
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.01                  ; contour spacing
  res@cnMinLevelValF       = -0.2                 ; min level
  res@cnMaxLevelValF       = 0.2                 ; max level
  fld@long_name="TABS"
  plot = ZeroLineContour(plot)
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)



end
