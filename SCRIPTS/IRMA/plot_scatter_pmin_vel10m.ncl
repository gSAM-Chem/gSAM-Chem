ncfiles=systemfunc("ls -1 ../../OUT_2D/IRMA_9216x4608x72_10s_u25w30_KKau_qi5_icef3_eis_2017-09-0[5678]*_atm.nc")
f = addfiles(ncfiles,"r")
nfiles=dimsizes(ncfiles)
print(""+ncfiles)
lonmin = 274.
lonmax = 295.
latmin = 15.
latmax = 25.

;g = addfile("../../OUT_INV/IRMA_9216x4608x74_10s_4608_TERR_MASKS.atm.3D.nc","r")

pmin = new((/nfiles/),float)
vel10max = new((/nfiles/),float)
j=0
;do i=24,nfiles-1
do i=0,nfiles-1
 p = f[i]->PSFC(0,{latmin:latmax},{lonmin:lonmax}) 
 time = f[i]->time
 land = f[i]->LANDMASK(0,{latmin:latmax},{lonmin:lonmax}) 
 p = where(land.lt.1.,p,9999.)
 pmin(j) = min(p)
 u = f[i]->U10m(0,{latmin:latmax},{lonmin:lonmax}) 
 v = f[i]->V10m(0,{latmin:latmax},{lonmin:lonmax}) 
 vel = sqrt(u^2+v^2)
 vel = where(land.lt.1.,vel,0.)
 vel10max(j) = max(vel)
 print(""+ncfiles(i))
 print(j+" "+time+" "+pmin(j)+"  "+vel10max(j))
 j=j+1
end do

 wks  = gsn_open_wks("ps","ps")             ; send graphics to PNG file
 
 res                     = True                   ; plot mods desired
 res@gsnMaximize         = True                   ; maximize plot in frame
 res@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
 res@xyMarkers           = 16                     ; choose type of marker 
 res@xyMarkerColor       = "black"                  ; Marker color
 res@xyMarkerSizeF       = 0.006                  ; Marker size (default 0.01)
 res@xyDashPatterns      = 1                      ; solid line 
 res@xyLineThicknesses   = (/1,2/)                ; set second line to 2
 res@vpWidthF            = 0.5           ; change aspect ratio of plot
 res@vpHeightF           = 0.3


 efit  = exponential_curve_fit(vel10max(0:j-1),pmin(0:j-1), 0, False)
 A     = efit(0)                                         ; extract for clarity
 B     = efit(1)
 yfit  = A*exp(B*vel10max(0:j-1))
 data = new((/2,j/),float)
 data(0,0:j-1) = pmin(0:j-1)
 data(1,0:j-1) = yfit
 plot  = gsn_csm_xy (wks,vel10max(0:j-1),data(0,0:j-1),res)        ; create plot


