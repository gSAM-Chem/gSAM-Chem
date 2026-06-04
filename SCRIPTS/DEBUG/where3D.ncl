f=addfile("../../OUT_3D/ISABEL_1024x1024x72_3km_10s_NONUDGE_DELTA1200_2003-09-20-00-00-00_0000034560.3D_atm.nc","r")
fld = f->U(0,:,:,:)

dims = dimsizes(fld)
nz = dims(0)
nlat = dims(1)
nlon = dims(2)
print(""+dims)
z = f->z
p = f->p
do k=0,nz-1
  fmax= max(fld(k,:,:))
  fmin= min(fld(k,:,:))
  print(k+"  "+z(k)+" "+fmin+" "+fmax)
end do
