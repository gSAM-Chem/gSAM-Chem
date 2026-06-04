f=addfile("/glade/campaign/univ/unsb0021/SOCRATES_OUT_SAM7.6/SAM_9216/SAM1MOM/OUT_3D/SOCRATES_9216x4608x74_SAM1MOM_4608_20180221020000_0000043920_TABS.nc","r")
g=addfile("/glade/campaign/univ/unsb0021/SOCRATES_OUT_SAM7.6/SAM_9216/SAM1MOM/OUT_3D/SOCRATES_9216x4608x74_SAM1MOM_4608_20180221020000_0000043920_TERRA.nc","r")
;f=addfile("../../OUT_3D/SOCRATES_9216x4608x74_M2005_4608_20180217010000_0000009000_TABS.nc","r")
;g=addfile("../../OUT_3D/SOCRATES_9216x4608x74_M2005_4608_20180217010000_0000009000_TERRA.nc","r")

z = f->z
do k=0,dimsizes(z)-1

fld = f->TABS(0,k,:,:)
terr = g->TERRA(0,k,:,:)

dims = dimsizes(fld)
nlat = dims(0)
nlon = dims(1)
print(min(fld)+"  "+max(fld))
do j=0,nlat-1
 fmax= max(fld(j,:))
 fmin= min(fld(j,:))
 if(fmin.lt.10.) then
   nmin = minind(fld(j,:))
   print(z(k)+" "+nmin+"  "+j+"  "+fld(j,nmin)+"  "+terr(j,nmin)+"  ")
 end if
 if(fmax.gt.400.) then
   nmax = maxind(fld(j,:))
   print(z(k)+" "+nmax+"  "+j+"  "+fld(j,nmax)+"  "+terr(j,nmax)+"  ")
 end if
end do

end do
