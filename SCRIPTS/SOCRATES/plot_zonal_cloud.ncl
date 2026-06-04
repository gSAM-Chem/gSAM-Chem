
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;----------------------------------------------------------------------

filename = "../../OUT_3D/SOCRATES_9216x4608x74_SAM1MOM_4608_20180217020000_0000009360"

  wks   = gsn_open_wks ("ps", "ps" )          ; send graphics to PNG file


  res                      = True                 ; plot mods desired
  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
;  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
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

  fin     = addfile (filename+"_QC.nc","r")
  fld     = dim_avg(fin->QC(0,:,:,:))
  fin     = addfile (filename+"_QI.nc","r")
  fld     = fld + dim_avg(fin->QI(0,:,:,:))

  print("qc+qi:"+min(fld)+max(fld))
  lev = fin->p
  lat = fin->lat
  fld!0 = "lev"
  fld&lev = lev
  fld!1 = "lat"
  fld&lat = lat
;  res@cnLevelSpacingF      = 5.0                  ; contour spacing
;  res@cnMinLevelValF       =-50.                 ; min level
;  res@cnMaxLevelValF       = 50.                 ; max level
  fld@long_name="U"
;  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder

  cld = fld
  cld = 0.
  delete(fld)


  fin1     = addfile (filename+"_QC.nc","r")
  fin2     = addfile (filename+"_QI.nc","r")
  p = fin1->p
  do k=0,dimsizes(p)-1
    if(p(k).lt.70.) then
         break
    end if
    fld     = fin1->QC(0,k,:,:) + fin2->QI(0,k,:,:)
    cld(k,:) = dim_avg(where(fld.gt.0., 1.,0.))
    print("k= "+k+"  "+p(k)+" "+max(cld(k,:)))
  end do

  lev = fin1->p
  lat = fin1->lat
  cld!0 = "lev"
  cld&lev = lev
  cld!1 = "lat"
  cld&lat = lat
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLevelSpacingF      = 0.02                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 0.12                 ; max level
  cld@long_name="CLD"
  plot  = gsn_csm_pres_hgt(wks, cld, res )   ; place holder
  delete(fld)




