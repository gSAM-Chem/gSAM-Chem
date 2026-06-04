f=addfile("../../OUT_2D/SOCRATES_9216x4608x74_SAM1MOM_4608_20180217020000_0000009360.2Dcom_1.nc","r")
g=addfile("../../OUT_INV/SOCRATES_9216x4608x74_SAM1MOM_4608_INV.2DLcom_1.nc","r")
fld = f->CFLZ(0,:,:)
fld1 = f->ZCFLZ(0,:,:)
z = g->ZSFC(0,:,:)

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
 if(fmax.gt.0.6) then
   nmax = maxind(fld(j,:))
   print(nmax+"  "+j+"  "+fld(j,nmax)+"  "+fld1(j,nmax))
   do jj=j-2,j+2
      print(""+z(jj,nmax-2:nmax+2)) 
   end do
 end if
end do
