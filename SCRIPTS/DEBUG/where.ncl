f=addfile("../.././OUT_3D/TROP_4608x1152x72_SAM1MOM_4608_20200120042550_0000001595.atm.3D.nc","r")
fld = f->TABS(0,2,:,:)

dims = dimsizes(fld)
nlat = dims(0)
nlon = dims(1)
print(""+dimsizes(fld))
print(min(fld)+"  "+max(fld))
do j=0,nlat-1
 fmax= max(fld(j,:))
 fmin= min(fld(j,:))
 if(fmin.lt.-900.) then
   nmin = minind(fld(j,:))
   print(j+"  "+nmin+"  "+fld(j,nmin))
 end if
 if(fmax.gt.0.3) then
   nmax = maxind(fld(j,:))
   print(nmax+"  "+j+"  "+fld(j,nmax))
 end if
end do
