;filename_root = "IRMA_9216x4608x74_10s_4608_20170911000000_0000051840"
filename_root = "IRMA_9216x4608x74_10s_mom23_0.8_4608_20170911000000_0000051840"
lonmin = 264.
lonmax = 295.
latmin = 10.
latmax = 30.

f=addfile("../../OUT_2D/"+filename_root+".2D_atm.nc","r")
landmask = f->LANDMASK(0,{latmin:latmax},{lonmin:lonmax})
psfc = f->PSFC(0,{latmin:latmax},{lonmin:lonmax})
psfc = where(landmask.lt.1.,psfc,9999.)

p1D     = ndtooned(psfc)
dsizes = dimsizes(psfc)
inds  = ind_resolve(minind(p1D),dsizes)
ilat0 = inds(0,0)
ilon0 = inds(0,1)
psfc0 = psfc(ilat0,ilon0)
print("psfc0="+psfc0)

lat = psfc&lat
lon = psfc&lon
lat0 = lat(ilat0)
lon0 = lon(ilon0)
print("lat0 lon0:"+lat0+" "+lon0)

deg2rad = acos(-1.)/180.
x = lon*cos(lat0*deg2rad)*40000./360.
y = lat*40000./360.
x0 = x(ilon0)
y0 = y(ilat0)
print("x0 y0:"+x0+" "+y0)

x = x-x0
x@units = "km"
y = y-y0
y@units = "km"

u = f->U850(0,{latmin:latmax},{lonmin:lonmax})
print("u:"+min(u)+" "+max(u))
v = f->V850(0,{latmin:latmax},{lonmin:lonmax})
print("v:"+min(v)+" "+max(v))
p = f->PSFC(0,{latmin:latmax},{lonmin:lonmax})
tb = f->TB(0,{latmin:latmax},{lonmin:lonmax})
;p!0 = "y"
;p!1 = "x"
;p&x = x
;p&y = y

print("x:"+min(x)+" "+max(x))
print("y:"+min(y)+" "+max(y))

;u!0 = "y"
;u!1 = "x"
;u&x = x
;u&y = y
;v!0 = "y"
;v!1 = "x"
;v&x = x
;v&y = y

vel = sqrt(u^2+v^2)
vel!0 = "lat"
vel!1 = "lon"
vel&lon = p&lon
vel&lat = p&lat

 wks  = gsn_open_wks("ps","ps")   ; send graphics to PNG file
  res                     = True          ; plot mods desired
  res@cnFillOn            = True          ; turn on color
  res@vpWidthF            = 0.6           ; change aspect ratio of plot
  res@vpHeightF           = 0.6
  res@gsnAddCyclic                = False  
;  res@gsnYAxisIrregular2Linear = True
;  res@trYMinF=lat0-300./40000.*360.
;  res@trYMaxF=lat0+300./40000.*360.
;  res@trXMinF=lon0-300./40000.*360./cos(lat0*deg2rad)
;  res@trXMaxF=lon0+300./40000.*360./cos(lat0*deg2rad)
 res@mpLimitMode                 = "Corners"            ; method to zoom
 res@mpLeftCornerLatF            = lat0-500./40000.*360.
 res@mpLeftCornerLonF            = lon0-500./40000.*360./cos(lat0*deg2rad)
 res@mpRightCornerLatF           = lat0+500./40000.*360.
 res@mpRightCornerLonF           = lon0+500./40000.*360./cos(lat0*deg2rad)
 res@mpDataBaseVersion           = "MediumRes"

;  res@trYMinF=-300.
;  res@trYMaxF=300.
;  res@trXMinF=-300.
;  res@trXMaxF=300.
  res@cnLinesOn = False
  res@cnRasterModeOn = True

  res@cnLevelSelectionMode = "ManualLevels"

  res@cnFillPalette       = "NCV_bright" ; set color map
  res@cnLevelSpacingF      = 2.                  ; contour spacing
  res@cnMinLevelValF       = 210.                 ; min level
  res@cnMaxLevelValF       =  310.                 ; max level
  plot = gsn_csm_contour_map(wks,tb,res)       ; contour the variable

  res@cnLevelSpacingF      = 2.                  ; contour spacing
  res@cnMinLevelValF       = 0.                 ; min level
  res@cnMaxLevelValF       = 60.                 ; max level
  plot = gsn_csm_contour_map(wks,vel,res)       ; contour the variable

  res@cnLevelSpacingF      = 2.                  ; contour spacing
  res@cnMinLevelValF       = 930.                 ; min level
  res@cnMaxLevelValF       =  1010.                 ; max level
  plot = gsn_csm_contour_map(wks,p,res)       ; contour the variable





