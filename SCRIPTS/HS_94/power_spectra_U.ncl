
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

f=addfile("../../OUT_3D/HS94_288x144x30_20s_144_U.nc","r")

;lats = (/-60./) ; to plot
lats = (/-60.,-45.,-30.,0.,30.,45.,60./) ; to plot

fld = f->U(:,16,:,:)
print("fld:"+min(fld)+"  "+max(fld))

pi = acos(-1.)

  wks = gsn_open_wks("ps","ps")
  res=True
  res@xyXStyle ="Log"
  res@xyYStyle ="Log"
  res@tmXBMinorPerMajor = 8
  res@tmXTMinorPerMajor = 8
  res@tmYRMinorPerMajor = 8
  res@tmYLMinorPerMajor = 8
  res@tmYLLabelFont = 21    ;  21 = helvetica
  res@tmXBLabelFont = 21    ;  22 = helvetica-bold

  res@tiXAxisString          = "wavenumber (n)"
 ; res@tiYAxisString          = "Kinetic energy, Power (m~F35~:S:2:N:/~F21~s~F35~:S:2:N:)"
  res@tiYAxisString          = ""
 ; res@gsnFrame = False
  res@xyLineThicknessF = 3.

  res@trXMinF = 1
  res@trXMaxF = 10000
  res@trYMinF = 1e-6
  res@trYMaxF = 1000
do j=0,dimsizes(lats)-1
  x = fld&lon/360.*40000.*cos(lats(j)/180.*pi)
  print("lat:"+lats(j)+"  max circum:"+max(x))
  coefs = ezfftf(fld(:,{lats(j)},:))
  power = coefs(0,:,:)^2+coefs(1,:,:)^2
  power!0="time"
  power!1="wave"
  powerm = dim_avg(power(wave|:,time|:))
  dx = x(1)-x(2)
  wavenumber = powerm
  do i=0,dimsizes(powerm)-1
    wavenumber(i)=i+1
  end do
  powerm(dimsizes(powerm)-1)=powerm(dimsizes(powerm)-1)*4.
  print(powerm)


  res@tiMainString           = "latitude "+lats(j)
  res@xyDashPattern = 0
  colors = (/"blue","red","green"/)
  res@xyLineColors = colors(1)
  plot  = gsn_xy (wks,wavenumber,powerm,res)     ; day 2

end do
