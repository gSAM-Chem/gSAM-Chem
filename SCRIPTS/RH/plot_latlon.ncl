    f=addfile("../../OUT_2D/RH_576x288x30_30s_corz_576_0000069120.2D_atm.nc","r")

    wks  = gsn_open_wks("ps","ps")             ; send graphics to PNG file
    gsn_define_colormap(wks,"BlAqGrYeOrReVi200")        ; choose colormap
    
    res                     = True
    res@gsnSpreadColors     = True
    res@lbLabelBarOn        = True            ; turn off individual label bars
    res@mpFillOn            = False            ; no need   
    res@cnLevelSelectionMode= "ManualLevels"   ; manual set levels
 ;   res@cnLevelSelectionMode= "AutomaticLevels"   ; manual set levels
    res@cnFillOn            = True             ; color fill plot
    res@cnLinesOn           = False
    res@cnLineLabelsOn      = False
    res@cnInfoLabelOn       = False
  res@tmXBLabelFontHeightF = 0.02       ; label font height
  res@tmYLLabelFontHeightF = 0.02
  res@lbLabelFontHeightF      = 0.02     ; label bar font
    res@mpOutlineOn = False
    res@mpCenterLonF         = 180.               ; center plot at 180
;    res@mpMinLonF            = 30.               ; select a subregion
;    res@mpMaxLonF            = 330.
;    res@mpMinLatF            = 0. 
;    res@mpMaxLatF            = 90.

    fld =f->U850(0,:,:)
    res@cnMinLevelValF      = 0.
    res@cnMaxLevelValF      = 24.
    res@cnLevelSpacingF     = 2.              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    fld =f->V850(0,:,:)
    res@cnMinLevelValF      = -15.
    res@cnMaxLevelValF      = 15.
    res@cnLevelSpacingF     = 2.5              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    fld =f->PSFC(0,:,:)
    res@cnMinLevelValF      = 960.
    res@cnMaxLevelValF      = 1040.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels        
    fld@long_name = "surface pressure"
    fld@units = "hPa"
    plot = gsn_csm_contour_map_ce(wks,fld,res)

