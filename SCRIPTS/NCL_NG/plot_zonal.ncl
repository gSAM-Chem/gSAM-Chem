
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;----------------------------------------------------------------------
begin

filename_root = "../../OUT_3D/NG_2048x1024x34_19.5x19km_20s_CTRL_1024_"
ind_start = 259200
ind_end = 388800
ind_step = 4320

  wks   = gsn_open_wks ("ps", "ps" )          ; send graphics to PNG file


  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = False                 ; turn on line labels
  res@cnLinesOn       = False                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
  res@cnFillPalette        = "BlWhRe"             ; choose colormap
  res@tmYRMode             = "Automatic"
 res@gsnDraw     = False
 res@gsnFrame    = False

plot =new(4,graphic)


do i = ind_start, ind_end, ind_step

  if(i.lt.1000) then
   filename = filename_root+"0000000"+i 
  else if(i.lt.10000) then
   filename = filename_root+"000000"+i 
  else if(i.lt.100000) then
   filename = filename_root+"00000"+i 
  else if(i.lt.1000000) then
   filename = filename_root+"0000"+i 
  end if
  end if
  end if
  end if
  print(""+filename)


  fin     = addfile (filename+"_U.nc","r")
  f3d     = fin->U(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  lat = fin->lat
  fld!0 = "lev"
  fld&lev = lev
  fld!1 = "lat"
  fld&lat = lat
  res@cnLevelSpacingF      = 5.0                  ; contour spacing
  res@cnMinLevelValF       = -70.                 ; min level
  res@cnMaxLevelValF       =  70.                 ; max level
  fld@long_name="U"
  plot(0)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

  fin     = addfile (filename+"_V.nc","r")
  f3d     = fin->V(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 1.0                  ; contour spacing
  res@cnMinLevelValF       = -10.                 ; min level
  res@cnMaxLevelValF       =  10.                 ; max level
  fld@long_name="V"
  plot(1)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

;  fin     = addfile (filename+"_W.nc","r")
;  f3d     = fin->W(0,:,:,:)
;  fld = dim_avg(f3d)
;  lev = fin->p
;  fld!0 = "lev"
;  fld&lev = lev
;  res@cnLevelSpacingF      = 0.01                  ; contour spacing
;  res@cnMinLevelValF       = -0.1                 ; min level
;  res@cnMaxLevelValF       =  0.1                 ; max level
;  fld@long_name="W"
;  plot(2)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

;  fin     = addfile (filename+"_QV.nc","r")
;  f3d     = fin->QV(0,:,:,:)
;  fld = dim_avg(f3d)
;  lev = fin->p
;  fld!0 = "lev"
;  fld&lev = lev
;  res@cnLevelSpacingF      = 0.1                  ; contour spacing
;  res@cnMinLevelValF       = 0.                 ; min level
;  res@cnMaxLevelValF       = 2.                 ; max level
;  fld@long_name="QV"
;  plot(1)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

;  fin     = addfile (filename+"_TT.nc","r")
;  f3d     = fin->TT(0,:,:,:)
;  fld = dim_avg(f3d)
;  lev = fin->p
;  fld!0 = "lev"
;  fld&lev = lev
;  res@cnLevelSpacingF      = 5.                  ; contour spacing
;  res@cnMinLevelValF       = 300.                 ; min level
;  res@cnMaxLevelValF       = 500.                 ; max level
;  fld@long_name="TT"
;  plot(2)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

;  fin     = addfile (filename+"_QRAD.nc","r")
;  f3d     = fin->QRAD(0,:,:,:)
;  fld = dim_avg(f3d)
;  lev = fin->p
;  fld!0 = "lev"
;  fld&lev = lev
;  res@cnLevelSpacingF      = 0.1                  ; contour spacing
;  res@cnMinLevelValF       = -2.0                 ; min level
;  res@cnMaxLevelValF       =  2.0                 ; max level
;  fld@long_name="QRAD"
;  plot(2)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

  fin     = addfile (filename+"_QN.nc","r")
  f3d     = fin->QN(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.01                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 0.1                 ; max level
  fld@long_name="QN"
  plot(2)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

  fin     = addfile (filename+"_TABS.nc","r")
  f3d     = fin->TABS(0,:,:,:)
  fld = dim_avg(f3d)
  lev = fin->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 5.                  ; contour spacing
  res@cnMinLevelValF       = 180.                 ; min level
  res@cnMaxLevelValF       = 310.                 ; max level
  fld@long_name="TABS"
  plot(3)  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

gsn_panel(wks,plot,(/2,2/),False)

end do

end
