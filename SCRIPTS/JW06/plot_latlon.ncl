    fi=addfile("../../OUT_3D/JW06_512x256x30_20s_BARO_r8_nadv2_nosgs2_256_0000038880.atm.3D.nc","r")
;    fi=addfile("../../OUT_3D/JW06_512x256x30_20s_BARO_r8_nadv2_nosgs_256_0000038880.atm.3D.nc","r")

    n850mb = 6
    wks  = gsn_open_wks("ps","ps")             ; send graphics to PNG file
  ;  gsn_define_colormap(wks,"BlAqGrYeOrRe")        ; choose colormap
  ;  gsn_define_colormap(wks,"grads_rainbow")        ; choose colormap
  ;  gsn_define_colormap(wks,"cosam")        ; choose colormap
  ;  gsn_define_colormap(wks,"cmocean_balance")        ; choose colormap
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
    res@tmXBLabelFontHeightF = 0.015       ; label font height
    res@tmYLLabelFontHeightF = 0.015
    res@lbLabelFontHeightF      = 0.015     ; label bar font
    res@pmLabelBarOrthogonalPosF = 0.2

    res@mpOutlineOn = False
    res@mpCenterLonF         = 180.               ; center plot at 180
    res@mpMinLonF            = 60.               ; select a subregion
    res@mpMaxLonF            = 330.
    res@mpMinLatF            = 0. 
    res@mpMaxLatF            = 90.

    rho = fi->rho
    p = fi->p
    fld =fi->PP(0,0,:,:)
    fld = fld*0.01
    fld&lon@units="degrees east"
    fld&lat@units="degrees north"
;   fld0 = conform(fld,dim_avg_Wrap(fld(:,{0.:30.})),0)
    fld0 = conform(fld,fld(:,0),0)
    fld = fld - fld0 + 1000.
    res@cnMinLevelValF      = 940.
    res@cnMaxLevelValF      = 1020.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels        

    fld@units = "hPa"
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    fld =fi->TABS(0,n850mb,:,:)
    fld&lon@units="degrees east"
    fld&lat@units="degrees north"
    res@cnMinLevelValF      = 220.
    res@cnMaxLevelValF      = 300.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    fld =fi->VOR(0,n850mb,:,:)
    fld = fld*1.e5
    fld&lon@units="degrees east"
    fld&lat@units="degrees north"
    res@cnMinLevelValF      = -10.
    res@cnMaxLevelValF      = 40.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels        
    fld@units = "10^-5 1/s"
    plot = gsn_csm_contour_map_ce(wks,fld,res)
