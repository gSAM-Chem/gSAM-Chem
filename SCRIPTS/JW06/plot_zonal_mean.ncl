
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;----------------------------------------------------------------------
begin

filename = "../../OUT_3D/JW06_512x256x30_20s_r8_nadv2_nosgs_256_0000129600.atm.3D.nc"

  wks   = gsn_open_wks ("ps", "ps" )          ; send graphics to PNG file


  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = False                 ; turn on line labels
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
  res@cnLevelSpacingF      = 2.0                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 40.                 ; max level
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
  res@cnLevelSpacingF      = 0.01                  ; contour spacing
  res@cnMinLevelValF       = -0.1                 ; min level
  res@cnMaxLevelValF       =  0.1                 ; max level
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

  f3d     = fin->PP(0,:,:,:)
  print(min(f3d)+"  "+max(f3d))
  fld = dim_avg(f3d)
  rho = fin->rho
  z = fin->z
  do k=0,dimsizes(fin->p)-1
    fld(k,:) = fld(k,:)/rho(k)/9.81 + z(k)
  end do
  fld = fld*0.001
  pi = acos(-1.)
  wt = cos(lat*pi/180.)
  wt = wt/sum(wt)
  fld0 = 0. 
  do i=0,dimsizes(lat)-1
   fld0 = fld0+fld(0,i)*wt(i)
  end do
  print("fld0="+fld0)
  fld = fld - fld0
  print(min(fld)+"  "+max(fld))
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.5                  ; contour spacing
  res@cnMinLevelValF       = -0.5                 ; min level
  res@cnMaxLevelValF       = 40.                 ; max level
  fld@long_name="Geopotential Height"
  res@cnLineLabelsOn       = True
  res@cnLineLabelBackgroundColor = -1
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  res@gsnContourZeroLineThicknessF = 3
  plot = ZeroLineContour(plot)
  draw(plot)
  frame (wks)
  res@cnLineLabelsOn       = False
  delete(fld0)

  f3d     = fin->PP(0,:,:,:)
  print(min(f3d)+"  "+max(f3d))
  fld = dim_avg(f3d)
  rho = fin->rho
  z = fin->z
  p = fin->p
  pi = acos(-1.)
  wt = cos(lat*pi/180.)
  wt = wt/sum(wt)
  fld0 = z
  fld0 = 0.
  do i=0,dimsizes(z)-1
   fld0(i) = fld0(i)+sum(fld(i,:)*wt(:))
  end do
  print("fld0="+fld0)
  do i=0,dimsizes(z)-1
    fld(i,:) = 0.01*(fld(i,:) - fld0(i)) + p(i)
  end do
  print(min(fld)+"  "+max(fld))
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 25.                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 1100.                 ; max level
  fld@long_name="Pressure"
  res@cnLineLabelsOn       = True
  res@cnLineLabelBackgroundColor = -1
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  res@cnLineLabelsOn       = False



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
