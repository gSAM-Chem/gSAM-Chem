f=addfile("../../OUT_MISC/IRMA_9216x4608x72_10s_newCFL_MISC_2017-09-05-12-00-00_0000004320.2D_atm.nc","r")
fld = f->CFL(0,:,:)

dims = dimsizes(fld)
nlat = dims(0)
nlon = dims(1)
print(""+dimsizes(fld))
print(min(fld)+"  "+max(fld))
do j=0,nlat-1
 fmax= max(fld(j,:))
 if(fmax.gt.0.) then
   nmax = maxind(fld(j,:))
   print(j+"  "+nmax+"  "+fld(j,nmax))
 end if
end do
