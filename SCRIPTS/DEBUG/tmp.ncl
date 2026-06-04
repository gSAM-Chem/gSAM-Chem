f=addfile("../../OUT_3D/out.nc","r")
fld = f->U(0,:,:,:)

dims = dimsizes(fld)
nz = dims(0)
nlat = dims(1)
nlon = dims(2)
print(""+dims)
p = f->level
do k=0,nz-1
  fmax= max(fld(k,:,:))
  fmin= min(fld(k,:,:))
  print(k+"  "+p(k)+" "+fmin+" "+fmax)
end do
