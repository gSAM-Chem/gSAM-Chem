
ncfiles=systemfunc("ls ../../OUT_3D/HS94_288x144x30_20s_144_00*.nc")
nfiles=dimsizes(ncfiles)

f=addfiles(ncfiles,"r")

time=f[:]->time
lon=f[0]->lon
lat=f[0]->lat
z=f[0]->z

ntime=dimsizes(time)
nz=dimsizes(z)
nlat=dimsizes(lat)
nlon=dimsizes(lon)

print("ntime= "+ntime)
print("nz= "+nz)
print("nlat= "+nlat)
print("nlon= "+nlon)

counttime = 0
nfiles_start = 0
do ifile=nfiles_start,nfiles-1
   print ( (/ncfiles(ifile)/) )

   ut:=f[ifile]->U
    uprime:=dim_rmvmean_Wrap(ut) ; deviation from the zonal mean

    if ( ifile.eq.nfiles_start ) then 
      vptp=new( (/ntime,nz,nlat/),typeof(ut))
      upvp=new( (/ntime,nz,nlat/),typeof(ut))
      tptp=new( (/ntime,nz,nlat/),typeof(ut))
      u=new( (/ntime,nz,nlat/),typeof(ut))
      v=new( (/ntime,nz,nlat/),typeof(ut))
      w=new( (/ntime,nz,nlat/),typeof(ut))
      t=new( (/ntime,nz,nlat/),typeof(ut))
    end if

    vt:=f[ifile]->V
    vprime:=dim_rmvmean_Wrap(vt) ; deviation from the zonal mean

    wt:=f[ifile]->W
    wprime:=dim_rmvmean_Wrap(wt) ; deviation from the zonal mean

    tt:=f[ifile]->TABS
    tprime:=dim_rmvmean_Wrap(tt) ; deviation from the zonal mean

    time:=f[ifile]->time

    ;calculate v'T'
    vptptemp=uprime
    vptptemp=(/vprime*tprime/)

    ; calculate u'v'
    upvptemp=uprime
    upvptemp=(/uprime*vprime/)

    ; calculate T'T'
    tptptemp=uprime
    tptptemp=(/tprime*tprime/)


    ;calculate zonal mean of U and eddy fluxes
    vptpzm=dim_avg_Wrap(vptptemp(time|:, z|:, lat|:, lon|:))
    tptpzm=dim_avg_Wrap(tptptemp(time|:, z|:, lat|:, lon|:))
    upvpzm=dim_avg_Wrap(upvptemp(time|:, z|:, lat|:, lon|:))
    uzm=dim_avg_Wrap(ut(time|:, z|:, lat|:, lon|:))
    vzm=dim_avg_Wrap(vt(time|:, z|:, lat|:, lon|:))
    wzm=dim_avg_Wrap(wt(time|:, zi|:, lat|:, lon|:))
    tzm=dim_avg_Wrap(tt(time|:, z|:, lat|:, lon|:))

    delete(vptptemp)
    delete(tptptemp)
    delete(upvptemp)
    delete(ut)
    delete(vt)
    delete(wt)
    delete(tt)

    ; accumulate u and eddy fluxes from each file
    vptp(counttime:counttime+dimsizes(time)-1,:,:)=vptpzm(:,:,:)
    tptp(counttime:counttime+dimsizes(time)-1,:,:)=tptpzm(:,:,:)
    upvp(counttime:counttime+dimsizes(time)-1,:,:)=upvpzm(:,:,:)
    u(counttime:counttime+dimsizes(time)-1,:,:)=uzm(:,:,:)
    v(counttime:counttime+dimsizes(time)-1,:,:)=vzm(:,:,:)
    w(counttime:counttime+dimsizes(time)-1,:,:)=wzm(:,:,:)
    t(counttime:counttime+dimsizes(time)-1,:,:)=tzm(:,:,:)

    delete(vptpzm)
    delete(tptpzm)
    delete(upvpzm)
    delete(uzm)
    delete(vzm)
    delete(wzm)
    delete(tzm)

    counttime=counttime+dimsizes(time)

end do

; take the time average
vptpclim=dim_avg_Wrap(vptp(z|:, lat|:, time|:))
tptpclim=dim_avg_Wrap(tptp(z|:, lat|:, time|:))
upvpclim=dim_avg_Wrap(upvp(z|:, lat|:, time|:))
uclim=dim_avg_Wrap(u(z|:, lat|:, time|:))
vclim=dim_avg_Wrap(v(z|:, lat|:, time|:))
wclim=dim_avg_Wrap(w(zi|:, lat|:, time|:))
tclim=dim_avg_Wrap(t(z|:, lat|:, time|:))


  ;set units and long_name for plotting purposes
  vptpclim@long_name="Northward eddy heat flux"
  vptpclim@units="Km/s"

  tptpclim@long_name="Eddy temperature variance"
  tptpclim@units="K^2"

  upvpclim@long_name="Northward eddy momentum flux"
  upvpclim@units="m^2/s^2"

  uclim@long_name="Zonal mean zonal wind"
  uclim@units="m/s"

  vclim@long_name="Zonal mean meridioanal wind"
  vclim@units="m/s"

  wclim@long_name="Zonal mean vertical wind"
  wclim@units="m/s"

  tclim@long_name="Zonal mean temperature"
  tclim@units="K"

  wks   = gsn_open_wks ("ps", "ps" )          ; send graphics to PNG file

  res                      = True                 ; plot mods desired
;  res@cnLevelSelectionMode = "AutomaticLevels"       ; manually select levels
  res@cnLevelSelectionMode = "ManualLevels"       ; manually select levels
  res@cnLineLabelsOn       = True                 ; turn on line labels
  res@cnLinesOn       = True                 ; turn on line labels
  res@cnFillOn             = True                 ; turn on color fill
;  res@cnFillPalette        = "BlueYellowRed"             ; choose colormap
  res@cnFillPalette        = "BlueWhiteOrangeRed"             ; choose colormap
  res@tmYRMode             = "Automatic"
  res@gsnYAxisIrregular2Linear = True
  res@trYReverse        = True
  res@gsnContourZeroLineThicknessF=0. ;don't plot zero line
  res@gsnContourNegLineDashPattern = 2
  res@cnLineThicknessF=3.

;  res@trYMinF = 50.
  res@gsnDraw         = False                   ; don't draw yet
  res@gsnFrame        = False                   ; don't advance frame yet

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


  fld = uclim
  lev = f[0]->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 2.0                  ; contour spacing
  res@cnMinLevelValF       = -40.                 ; min level
  res@cnMaxLevelValF       =  40.                 ; max level
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  delete(fld)


  fld = vclim
  lev = f[0]->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.5                  ; contour spacing
  res@cnMinLevelValF       = -2.                 ; min level
  res@cnMaxLevelValF       =  2.                 ; max level
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  delete(fld)


  fld = wclim
  lev = f[0]->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 0.001                  ; contour spacing
  res@cnMinLevelValF       = -0.01                 ; min level
  res@cnMaxLevelValF       =  0.01                 ; max level
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  delete(fld)


  fld = tclim
  lev = f[0]->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 5.0                  ; contour spacing
  res@cnMinLevelValF       = 180.                 ; min level
  res@cnMaxLevelValF       = 310.                 ; max level
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  delete(fld)

  fld = upvpclim
  lev = f[0]->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 8.0                  ; contour spacing
  res@cnMinLevelValF       = -120.                 ; min level
  res@cnMaxLevelValF       = 120.                 ; max level
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  delete(fld)

  fld = tptpclim
  lev = f[0]->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 10.0                  ; contour spacing
  res@cnMinLevelValF       = -40.                 ; min level
  res@cnMaxLevelValF       = 40.                 ; max level
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  delete(fld)

  fld = vptpclim
  lev = f[0]->p
  fld!0 = "lev"
  fld&lev = lev
  res@cnLevelSpacingF      = 2.0                  ; contour spacing
  res@cnMinLevelValF       = -40.                 ; min level
  res@cnMaxLevelValF       = 40.                 ; max level
  plot  = gsn_csm_contour(wks, fld, res )   ; place holder
;  plot = ZeroLineContour(plot)
;  plot  = gsn_csm_pres_hgt(wks, fld, res )   ; place holder
  draw(plot)
  frame (wks)
  delete(fld)


