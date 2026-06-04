
  f = addfile("../../OUT_2D/DYAMOND2_1440x720x72_10s_dyvar_CUP_auto1.2_1440_JAN21-25_2020_AVG.2D_atm.nc","r")
  f_obs = addfile("../../../obs_data/CERES2_02_climo.nc","r")


  land = tointeger(f->LAND(0,:,:))
  swnt = f->SWNTC(0,:,:)
  print("swnt:"+min(swnt)+" "+max(swnt))

  fld1 = f_obs->FSNTOAC(0,:,:)
  swnt_obs = linint2_Wrap(fld1&lon,fld1&lat,fld1,True,swnt&lon,swnt&lat,0)
  print("swnt_obs:"+min(swnt_obs)+" "+max(swnt_obs))

  wks   = gsn_open_wks ("ps","ps")           ; send graphics to PNG file

  res                   = True                     ; plot mods desired
  res@gsnMaximize       = True                     ; maximize plot
  res@xyMarkLineMode    = "Markers"                ; choose to use markers
  res@xyMarkers         =  16                      ; choose type of marker  
  res@xyMarkerColor     = "NavyBlue"               ; Marker color
  res@xyMarkerSizeF     = 0.001                     ; Marker size (default 0.01)
  res@pmLegendDisplayMode    = "Never"  

  res@trXMinF = 0.
  res@trXMaxF = 500.
  res@trYMinF = 0.
  res@trYMaxF = 500.

  do i=0,16
    res@tiMainString      = "LAND TYPE = "+i           ; add title
    print("landtype"+i)
    fld = ndtooned(where(land.eq.i,swnt,-999.))
    fld_obs = ndtooned(where(land.eq.i,swnt_obs,-999.))
    fld@_FillValue = -999.
    fld_obs@_FillValue = -999.
    print("fld mean:"+avg(fld))
    print("fld_obs mean:"+avg(fld_obs))
    plot  = gsn_csm_xy (wks,fld(::10),fld_obs(::10),res)                    ; create plot
  end do
