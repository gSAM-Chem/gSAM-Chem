
files = systemfunc("ls -1 ../../OUT_2D/IRMA_2304x1152x74_15s_mom23_0.8_ens00_0.4K_1152_2017*0.2D_atm.nc")
MICRO="IRMA"
nfiles = dimsizes(files)
npoints_max = 100000
dist = 300.
xp = new((/npoints_max/),float)
yp = xp
np = 0
do i=0,nfiles-1
  n = 0
  print(""+files(i))
  f=addfile(files(i),"r")
  pw = f->PSFC(0,{0.:30.},:)
  if(i.eq.0) then
   landmask = f->LANDMASK(0,{0.:30.},:) 
  end if
  pw = where(landmask.eq.0..and.pw.lt.995.,pw,9999999.)
  pw@_FillValue = 9999999.
  print("PW:"+min(pw)+" "+max(pw))
  locmin = local_min(pw,True,0.)
  if(locmin.eq.0) then
    continue
  end if
  lon = pw&lon
  lat = pw&lat
  x=lon(locmin@xi)                ; get lat/lon points of minima
  y=lat(locmin@yi)
  v=locmin@minval
  print("number of local minima:"+locmin)
  do k=0,locmin-1
   if(v(k).lt.995.) then
     xp(np) = x(k)
     yp(np) = y(k)
     if(np.gt.0) then
      if((xp(np)-xp(np-1))^2+(yp(np)-yp(np-1))^2.gt.300.) then
       print(k+"  "+v(k)+" "+x(k)+" "+y(k))
       np = np+1
       n = n+1
      end if
     else
      np = np+1
      n = n+1
     end if
   end if
  end do
  delete(locmin)
  delete(x)
  delete(y)
  delete(v)
  print("n="+n)
end do
print("np="+np)
print("mean="+np/nfiles)
;*****************************************
; plot original data
;*****************************************
  sst = f->SKT(0,:,:)
  wks  = gsn_open_wks("ps","ps")           ; send graphics to PNG file
  cmap = read_colormap_file("BlWhRe")           ; read color map file
  ncolors = dimsizes(cmap(:,0))                 ; get number of colors
  res  = True
  res@cnFillOn             = True               ; turn on color fill
  res@mpFillOn             = False              ; turn off gray continents
  res@cnLinesOn            = False
  res@mpFillOn = False
  res@gsnFrame             = False 
  sst@long_name = MICRO
  sst@units = np/(1.*nfiles)
  plot = gsn_csm_contour_map(wks,sst,res)
;*****************************************
; plot min extrema over original data
;*****************************************
  polyres               = True          ; poly marker mods desired 
  polyres@gsMarkerIndex = 16            ; choose circle as polymarker
  polyres@gsMarkerSizeF = 3.0           ; select size to avoid streaking
  polyres@gsMarkerColor = (/"black"/)   ; choose color
  gsn_polymarker(wks,plot,xp(0:np-1),yp(0:np-1),polyres)  ; draw polymarkers
  frame(wks)
