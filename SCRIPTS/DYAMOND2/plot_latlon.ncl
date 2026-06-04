    fi=addfile("../../OUT_2D/DYAMOND2_1440x720x74_10s_npres0_p00_topo_720_20200120010000_0000000360.2D_atm.nc","r")
    fo=addfile("/glade/scratch/marat/ps.nc","r")
    fo1 = addfile("/glade/collections/rda/data/ds633.0/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc","r")

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
    res@cnFillMode = "RasterFill"
    res@cnInfoLabelOn       = False
    res@tmXBLabelFontHeightF = 0.015       ; label font height
    res@tmYLLabelFontHeightF = 0.015
    res@lbLabelFontHeightF      = 0.015     ; label bar font
    res@pmLabelBarOrthogonalPosF = 0.2

    res@mpOutlineOn = False
    res@mpCenterLonF         = 0.               ; center plot at 180
;    res@mpMinLonF            = 60.               ; select a subregion
;    res@mpMaxLonF            = 330.
;    res@mpMinLatF            = 0. 
;    res@mpMaxLatF            = 90.

    fld =fi->PSFC(0,:,:)
    fld&lon@units="degrees east"
    fld&lat@units="degrees north"
    res@cnMinLevelValF      = 440.
    res@cnMaxLevelValF      = 1050.
    res@cnLevelSpacingF     = 20.              ; 20 contour levels        

    fld@units = "hPa"
    plot = gsn_csm_contour_map_ce(wks,fld,res)

    nlato=dimsizes(fo->latitude)
    fld1 =fo->SP(0,nlato-2:0:1,:)
    fld1 =fld1 +fo->SP(0,nlato-1:1:1,:)
    fld1 = 0.5*0.01*fld1
    res@cnMinLevelValF      = 440.
    res@cnMaxLevelValF      = 1050.
    res@cnLevelSpacingF     = 20.              ; 20 contour levels

    fld@units = "hPa"
    plot = gsn_csm_contour_map_ce(wks,fld1,res)


    fld2 = fld-fld1
    fld2!0="lat"
    fld2!1="lon"
    fld2&lat = fld&lat
    fld2&lon = fld&lon
    res@cnMinLevelValF      = -50.
    res@cnMaxLevelValF      = 50.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels
    plot = gsn_csm_contour_map_ce(wks,fld2,res)


    fld =fi->ZSFC(0,:,:)
    fld&lon@units="degrees east"
    fld&lat@units="degrees north"
    res@cnMinLevelValF      = 0.
    res@cnMaxLevelValF      = 7000.
    res@cnLevelSpacingF     = 100.              ; 20 contour levels

    fld@units = "m"
    plot = gsn_csm_contour_map_ce(wks,fld,res)


    nlato=dimsizes(fo->latitude)
    fld1 =fo1->Z(0,nlato-2:0:1,:)
    fld1 =fld1 +fo1->Z(0,nlato-1:1:1,:)
    fld1 = 0.5*fld1/9.79764
    res@cnMinLevelValF      = 0.
    res@cnMaxLevelValF      = 7000.
    res@cnLevelSpacingF     = 100.              

    fld@units = "m"
    plot = gsn_csm_contour_map_ce(wks,fld1,res)


    fld2 = fld-fld1
    fld2!0="lat"
    fld2!1="lon"
    fld2&lat = fld&lat
    fld2&lon = fld&lon
    res@cnMinLevelValF      = -500.
    res@cnMaxLevelValF      = 500.
    res@cnLevelSpacingF     = 5.              ; 20 contour levels
    plot = gsn_csm_contour_map_ce(wks,fld2,res)

