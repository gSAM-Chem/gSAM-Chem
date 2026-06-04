f=addfile("../../OUT_3D/ISABEL_1024x1024x72_3km_10s_NONUDGE_DELTA1200_2003-09-20-00-00-00_0000034560.3D_atm.nc","r")
fld = f->U(0,:,:,:)

dims = dimsizes(fld)
nz = dims(0)
ny = dims(1)
nx = dims(2)
z = f->z
do k=0,nz-1
print(k+"  "+min(fld(k,:,:))+"  "+max(fld(k,:,:)))
do j=0,ny-1
 fmax= max(fld(k,j,:))
 fmin= min(fld(k,j,:))
 if(fmin.lt.-40.) then
   nmin = minind(fld(k,j,:))
   print(k+"  "+z(k)+" "+j+" "+nmin+"  "+fld(k,j,nmin))
 end if
end do

end do
