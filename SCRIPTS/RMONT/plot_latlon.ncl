;    fi=addfile("../../OUT_3D/RMONT_576x288x74_30s_nadvmom23_0.5_576_0000014400.atm.3D.nc","r")
;    fi=addfile("../../OUT_3D/RMONT_576x288x74_30s_nadvmom23_0.5_576_0000043200.atm.3D.nc","r")
    fi=addfile("../../OUT_3D/RMONT_576x288x74_30s_nadvmom23_0.5_576_0000072000.atm.3D.nc","r")

    z = 3000.

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
    res@mpOutlineOn = False
    res@mpCenterLonF         = 180.               ; center plot at 180
;    res@mpMinLonF            = 30.               ; select a subregion
;    res@mpMaxLonF            = 330.
;    res@mpMinLatF            = 0. 
;    res@mpMaxLatF            = 90.

    fld =fi->U(0,{z},:,:)
    res@cnMinLevelValF      = -10.
    res@cnMaxLevelValF      = 45.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    fld =fi->V(0,{z},:,:)
    res@cnMinLevelValF      = -30.
    res@cnMaxLevelValF      = 15.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels        
    plot = gsn_csm_contour_map_ce(wks,fld,res)

