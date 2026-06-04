;files3D = systemfunc("ls -1 ../../OUT_3D/IRMA_9216x4608x72_10s_u25w30_KKau_qi5_icef3_eis_2017-09-08-00-00-00_0000025920_U.3D_atm.nc")
files3D = systemfunc("ls -1 ../../OUT_3D/IRMA_9216x4608x72_10s_u25w30_KKau_qi5_icef3_eis_2017-09-08*_U.3D_atm.nc")
print(files3D)
lonmin = 264.
lonmax = 295.
latmin = 10.
latmax = 30.
rad0 = 310.*1000. ;radius for plots
rad_out = 700.*1000. ; maximum radius for background field

nfiles = dimsizes(files3D)
;nfiles = 1

do m=0,nfiles-1

filename_root = str_get_cols(files3D(m),13,90)
print(""+filename_root)

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

drad = 5.*1000.
dphi = 5.
radmax = rad0
phimax = 360.
nrad = tointeger(radmax/drad)
nphi = tointeger(phimax/dphi)

deg2rad = acos(-1.)/180.
x = lon*cos(lat0*deg2rad)*40000./360.*1000.
y = lat*40000./360.*1000.
x0 = x(ilon0)
y0 = y(ilat0)
print("x0 y0:"+x0+" "+y0)


f=addfile("../../OUT_3D/"+filename_root+"_U.3D_atm.nc","r")
u = f->U(0,:,{latmin:latmax},{lonmin:lonmax})
print("u:"+min(u)+" "+max(u))
f=addfile("../../OUT_3D/"+filename_root+"_V.3D_atm.nc","r")
v = f->V(0,:,{latmin:latmax},{lonmin:lonmax})
print("v:"+min(v)+" "+max(v))
f=addfile("../../OUT_3D/"+filename_root+"_W.3D_atm.nc","r")
w = f->W(0,:,{latmin:latmax},{lonmin:lonmax})
print("w:"+min(w)+" "+max(w))
f=addfile("../../OUT_3D/"+filename_root+"_TABS.3D_atm.nc","r")
t = f->TABS(0,:,{latmin:latmax},{lonmin:lonmax})
print("t:"+min(t)+" "+max(t))
f=addfile("../../OUT_3D/"+filename_root+"_QV.3D_atm.nc","r")
qv = f->QV(0,:,{latmin:latmax},{lonmin:lonmax})
print("qv:"+min(qv)+" "+max(qv))
f=addfile("../../OUT_3D/"+filename_root+"_QC.3D_atm.nc","r")
qc = f->QC(0,:,{latmin:latmax},{lonmin:lonmax})
f=addfile("../../OUT_3D/"+filename_root+"_QI.3D_atm.nc","r")
qc = qc + f->QI(0,:,{latmin:latmax},{lonmin:lonmax})
print("qc:"+min(qc)+" "+max(qc))
f=addfile("../../OUT_3D/"+filename_root+"_QR.3D_atm.nc","r")
qr = f->QR(0,:,{latmin:latmax},{lonmin:lonmax})
f=addfile("../../OUT_3D/"+filename_root+"_QS.3D_atm.nc","r")
qr = qr + f->QS(0,:,{latmin:latmax},{lonmin:lonmax})
f=addfile("../../OUT_3D/"+filename_root+"_QG.3D_atm.nc","r")
qr = qr + f->QG(0,:,{latmin:latmax},{lonmin:lonmax})
print("qr:"+min(qr)+" "+max(qr))
;f=addfile("../../OUT_3D/"+filename_root+"_PP.3D_atm.nc","r")
;pp = f->PP(0,:,{latmin:latmax},{lonmin:lonmax})
;print("pp:"+min(pp)+" "+max(pp))

x = x-x0
y = y-y0

print("x:"+min(x)+" "+max(x))
print("y:"+min(y)+" "+max(y))

u!1 = "y"
u!2 = "x"
u&x = x
u&y = y
v!1 = "y"
v!2 = "x"
v&x = x
v&y = y
w!1 = "y"
w!2 = "x"
w&x = x
w&y = y
t!1 = "y"
t!2 = "x"
t&x = x
t&y = y
qv!1 = "y"
qv!2 = "x"
qv&x = x
qv&y = y
qc!1 = "y"
qc!2 = "x"
qc&x = x
qc&y = y
qr!1 = "y"
qr!2 = "x"
qr&x = x
qr&y = y
;pp!1 = "y"
;pp!2 = "x"
;pp&x = x
;pp&y = y

z = f->z
nz = dimsizes(z)
ny = dimsizes(y)
nx = dimsizes(x)

if(m.eq.0) then
 radvel = new((/nz,nphi,nrad/),float)
 azimvel = new((/nz,nphi,nrad/),float)
 vervel = new((/nz,nphi,nrad/),float)
 temp = new((/nz,nphi,nrad/),float)
 qvv = new((/nz,nphi,nrad/),float)
 qcc = new((/nz,nphi,nrad/),float)
 qrr = new((/nz,nphi,nrad/),float)
 ppr = new((/nz,nphi,nrad/),float)
 tempm = new((/nz/),float)
 ppm = new((/nz/),float)
 rads = new((/nrad/),float)
 radvel = 0.
 azimvel = 0.
 vervel = 0.
 temp = 0.
 tempm = 0.
 ppm = 0.
 qvv = 0.
 qcc = 0.
 qrr = 0.
 ppr = 0.
 radxy = new((/ny,nx/),float)
 do j=0,ny-1
  radxy(j,:) = tofloat(x(:)^2+y(j)^2)
 end do
 radxy = sqrt(radxy)
 radxy = where(radxy.gt.rad0.and.radxy.lt.rad_out,1.,0.)
 radxy = radxy/sum(radxy)
end if

do j=0,nphi-1
 phi = dphi*j
 do i=0,nrad-1
  rad = drad*i
  rads(i) = rad
  xx = rad*cos(phi*deg2rad)
  yy = rad*sin(phi*deg2rad)
;  print(rad+" "+phi+" "+xx+" "+yy)
  radvel(:,j,i) = radvel(:,j,i)+(u(:,{yy},{xx})*xx+v(:,{yy},{xx})*yy)/sqrt(xx^2+yy^2+1.e-5)
  azimvel(:,j,i) = azimvel(:,j,i)+(-u(:,{yy},{xx})*yy+v(:,{yy},{xx})*xx)/sqrt(xx^2+yy^2+1.e-5)
  vervel(:,j,i) = vervel(:,j,i)+ w(:,{yy},{xx})
  temp(:,j,i) = temp(:,j,i)+t(:,{yy},{xx})
  qvv(:,j,i) = qvv(:,j,i)+qv(:,{yy},{xx})
  qcc(:,j,i) = qcc(:,j,i)+qc(:,{yy},{xx})
  qrr(:,j,i) = qrr(:,j,i)+qr(:,{yy},{xx})
;  ppr(:,j,i) = ppr(:,j,i)+pp(:,{yy},{xx})
 end do
end do

do k=0,nz-1
 tempm(k) = tempm(k)+dim_sum_Wrap(dim_sum_Wrap(radxy*t(k,:,:)))
; ppm(k) = ppm(k)+dim_sum_Wrap(dim_sum_Wrap(radxy*pp(k,:,:)))
end do


delete(u)
delete(v)
delete(w)
delete(t)
delete(qv)
delete(qc)
delete(qr)

end do ; m

radvel = radvel/nfiles
azimvel = azimvel/nfiles
vervel = vervel/nfiles
temp = temp/nfiles
tempm = tempm/nfiles
;ppm = ppm/nfiles
do k=0,nz-1
 temp(k,:,:) = temp(k,:,:) - tempm(k)
; ppr(k,:,:) = ppr(k,:,:) - ppm(k)
end do
qvv = qvv/nfiles
qcc = qcc/nfiles
qrr = qrr/nfiles

radvel0 = dim_avg_n_Wrap(radvel,1)
azimvel0 = dim_avg_n_Wrap(azimvel,1)
vervel0 = dim_avg_n_Wrap(vervel,1)
temp0 = dim_avg_n_Wrap(temp,1)
qvv0 = dim_avg_n_Wrap(qvv,1)
qcc0 = dim_avg_n_Wrap(qcc,1)
qrr0 = dim_avg_n_Wrap(qrr,1)
qcc0 = qcc0 + qrr0
;ppr0 = dim_avg_n_Wrap(ppr,1)
print("rads: "+rads)
radvel0!0 = "z"
radvel0!1 = "x"
radvel0&z = z/1000.
radvel0&x = rads/1000.
azimvel0!0 = "z"
azimvel0!1 = "x"
azimvel0&z = z/1000.
azimvel0&x = rads/1000.
vervel0!0 = "z"
vervel0!1 = "x"
vervel0&z = z/1000.
vervel0&x = rads/1000.
temp0!0 = "z"
temp0!1 = "x"
temp0&z = z/1000.
temp0&x = rads/1000.
qvv0!0 = "z"
qvv0!1 = "x"
qvv0&z = z/1000.
qvv0&x = rads/1000.
qcc0!0 = "z"
qcc0!1 = "x"
qcc0&z = z/1000.
qcc0&x = rads/1000.
qrr0!0 = "z"
qrr0!1 = "x"
qrr0&z = z/1000.
qrr0&x = rads/1000.
;ppr0!0 = "z"
;ppr0!1 = "x"
;ppr0&z = z/1000.
;ppr0&x = rads/1000.

 wks  = gsn_open_wks("ps","ps")   ; send graphics to PNG file
  res                     = True          ; plot mods desired
  res@cnFillOn            = True          ; turn on color
  res@cnFillPalette       = "cmocean_balance" ; set color map
  res@vpWidthF            = 0.5           ; change aspect ratio of plot
  res@vpHeightF           = 0.3
  res@gsnYAxisIrregular2Linear = True
  res@cnLinesOn = True
  res@cnMonoLineColor  = False
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLineThicknessF = 3.0             ; Make lines thicker
  res@gsnContourZeroLineThicknessF = 5.


 do iter = 1,2

  if(iter.eq.1) then
   res@trYMinF=0.
   res@trYMaxF=20.
  else
   res@trYMinF=0.
   res@trYMaxF=2.
  end if
  res@trXMinF=0.
  res@trXMaxF=300.

  res@cnFillPalette       = "cmocean_balance" ; set color map
  res@cnLevelSpacingF      = 4.                  ; contour spacing
  res@cnMinLevelValF       = -28.                 ; min level
  res@cnMaxLevelValF       =  28.                 ; max level
  res@lbLabelStride         = 1
  res@tiMainString    = "Radial Wind"
  print("radvel: "+min(radvel0)+" "+max(radvel0))
  plot = gsn_csm_contour(wks,radvel0,res)       ; contour the variable

  res@cnFillPalette       = "cmocean_balance" ; set color map
  res@cnLevelSpacingF      = 5.                  ; contour spacing
  res@cnMinLevelValF       = -80.                 ; min level
  res@cnMaxLevelValF       =  80.                 ; max level
  res@lbLabelStride         = 2 
  res@tiMainString    = "Azimuthal Wind"
  print("azimvel: "+min(azimvel0)+" "+max(azimvel0))
  plot = gsn_csm_contour(wks,azimvel0,res)       ; contour the variable

  res@cnFillPalette       = "cmocean_balance" ; set color map
  res@cnLevelSpacingF      = 1.                  ; contour spacing
  res@cnMinLevelValF       = -5.                 ; min level
  res@cnMaxLevelValF       =  5.                 ; max level
  res@lbLabelStride         = 1
  res@tiMainString    = "Vertical Wind"
  print("vervel: "+min(vervel0)+" "+max(vervel0))
  plot = gsn_csm_contour(wks,vervel0,res)       ; contour the variable

  res@cnFillPalette       = "cmocean_balance" ; set color map
  res@cnLevelSpacingF      = 2.                  ; contour spacing
  res@cnMinLevelValF       = -18.                 ; min level
  res@cnMaxLevelValF       =  18.                 ; max level
  res@lbLabelStride         = 2 
  res@tiMainString    = "Temperature Perturbation"
  print("temp: "+min(temp0)+" "+max(temp0))
  plot = gsn_csm_contour(wks,temp0,res)       ; contour the variable

  res@cnFillPalette       = "WhViBlGrYeOrRe" ; set color map
  res@cnLevelSpacingF      =  1.                  ; contour spacing
  res@cnMinLevelValF       =  0.                 ; min level
  res@cnMaxLevelValF       = 25.                 ; max level
  res@lbLabelStride         = 1
  res@tiMainString    = "Water Vapor"
  print("qvv0: "+min(qvv0)+" "+max(qvv0))
  plot = gsn_csm_contour(wks,qvv0,res)       ; contour the variable

  res@cnFillPalette       = "WhViBlGrYeOrRe" ; set color map
  res@cnLevelSpacingF      =  .25                  ; contour spacing
  res@cnMinLevelValF       =  0.                 ; min level
  res@cnMaxLevelValF       = 4.                 ; max level
  res@lbLabelStride         = 2 
  res@tiMainString    = "Total Liquid/Ice Water"
  print("qcc0: "+min(qcc0)+" "+max(qcc0))
  plot = gsn_csm_contour(wks,qcc0,res)       ; contour the variable

;  res@cnFillPalette       = "NCV_bright" ; set color map
;  res@cnLevelSpacingF      = 10.                  ; contour spacing
;  res@cnMinLevelValF       = -80.                 ; min level
;  res@cnMaxLevelValF       =  0.                 ; max level
;  res@lbLabelStride         = 2 
;  plot = gsn_csm_contour(wks,ppr0,res)       ; contour the variable

end do

