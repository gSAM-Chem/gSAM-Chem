    fi=addfile("../../OUT_3D/JW06_576x288x74_20s_std_baro_heper256_288_0000038880.nc","r")
    wks  = gsn_open_wks("ps","ps")             ; send graphics to PNG file
  ;  gsn_define_colormap(wks,"BlAqGrYeOrRe")        ; choose colormap
    gsn_define_colormap(wks,"temp1")        ; choose colormap
    
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
    res@mpOutlineOn = False
    res@mpCenterLonF         = 0.               ; center plot at 180
    res@mpMinLonF            = 30.               ; select a subregion
    res@mpMaxLonF            = 330.
    res@mpMinLatF            = 0. 
    res@mpMaxLatF            = 90.

    fld =fi->PP(0,0,:,:)
    fld = fld*0.01
    fld0 = conform(fld,dim_avg_Wrap(fld(:,{0.:30.})),0)
  ;  fld0 = conform(fld,fld(:,0),0)
    fld = fld - fld0 + 1000.
    res@cnMinLevelValF      = 940.
    res@cnMaxLevelValF      = 1020.
    res@cnLevelSpacingF     = 10.              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    fld =fi->TABS(0,{1348.},:,:)
    res@cnMinLevelValF      = 220.
    res@cnMaxLevelValF      = 310.
    res@cnLevelSpacingF     = 2.              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    fld =fi->VOR(0,{1348.},:,:)
    fld = fld*1.e5
    res@cnMinLevelValF      = -10.
    res@cnMaxLevelValF      = 40.
    res@cnLevelSpacingF     = 2.              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)
