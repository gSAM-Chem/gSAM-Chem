f=addfile("../../OUT_3D/JW06_512x256x30_20s_BARO_r8_nadv2_nosgs_tmp_256_0000000001.atm.3D.nc","w")
f1=addfile("../../OUT_3D/JW06_512x256x30_20s_BARO_r8_nadv2_nosgs_tmp_256_0000000000.atm.3D.nc","w")
nx=dimsizes(f->lon)
u = f->U
u1 = f1->U
v = f->V
do i=0,nx-1
 u(:,:,:,i) = u(:,:,:,i) - u1(:,:,:,nx-1)
end do
f->U1 = u
f->V1 = v
