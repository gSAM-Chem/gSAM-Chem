    fi=addfile("../../OUT_3D/JW06_512x256x30_20s_BARO_r8_nadv2_nosgs2_256_0000060480.atm.3D.nc","r")
;    fi=addfile("../../OUT_3D/JW06_512x256x30_20s_BARO_r8_nadv2_nosgs2_256_0000047520.atm.3D.nc","r")
    wks  = gsn_open_wks("ps","ps")             ; send graphics to PNG file
  ;  gsn_define_colormap(wks,"BlAqGrYeOrRe")        ; choose colormap
    gsn_define_colormap(wks,"grads_rainbow")        ; choose colormap
    
    res                     = True
    res@gsnSpreadColors     = True
    res@lbLabelBarOn        = True            ; turn off individual label bars
    res@mpFillOn            = False            ; no need   
    res@mpOutlineOn  = False
    res@cnLevelSelectionMode= "ManualLevels"   ; manual set levels
;    res@cnLevelSelectionMode= "AutomaticLevels"   ; manual set levels
;    res@cnFillOn            = True             ; color fill plot
;    res@cnLinesOn           = True
;    res@cnLineLabelsOn      = False
;    res@cnInfoLabelOn       = False
;    res@mpOutlineOn = False
    res@mpCenterLonF         = -90.               ; center plot at 180
;    res@mpMinLonF            = 45.               ; select a subregion
;    res@mpMaxLonF            = 360.
;    res@mpMinLatF            = 0. 
;    res@mpMaxLatF            = 90.

    fld =fi->PP(0,0,:,:)
    fld = fld*0.01
    fld&lon@units="degrees east"
    fld&lat@units="degrees north"
    fld0 = conform(fld,dim_avg_Wrap(fld(:,{0.:30.})),0)
  ;  fld0 = conform(fld,fld(:,0),0)
    fld = fld - fld0 + 1010.
    res@gsnPolar   = "NH"                          ; specify the hemisphere
    res@cnMinLevelValF      = 820.
    res@cnMaxLevelValF      = 1180.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels        
    plot = gsn_csm_contour_map_polar(wks,fld,res)    ; create the plot

