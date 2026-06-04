    f=addfile("../../OUT_2D/out.nc","r")

    ny = dimsizes(f->lat)
    u = f->U850(:,ny/2,:)

    wks  = gsn_open_wks("ps","ps")             ; send graphics to PNG file

  res                      = True               ; plot mods desired
  res@cnFillOn             = True               ; turn on color fill
  res@cnFillPalette        = "BlWhRe"           ; set color map

  res@tiMainString         = "Pacific Region"   ; title

  res@cnLevelSelectionMode = "ManualLevels"     ; manual contour levels
  res@cnMinLevelValF       = 0.	        ; min level
  res@cnMaxLevelValF       = 24.               ; max level
  res@cnLevelSpacingF      = 2.                 ; contour level spacing
  
  plot = gsn_csm_hov(wks, u(:,:), res) 
  
