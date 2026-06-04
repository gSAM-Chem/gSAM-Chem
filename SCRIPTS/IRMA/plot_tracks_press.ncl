files = systemfunc("ls -1 ../../OUT_2D/IRMA_2304x1152x74_15s_mom23_0.8_ens0*_IRMA.2D_atm.nc")
print(""+files)
system("rm -f ../../OUT_2D/out.nc")
g=addfile("../../OUT_2D/out.nc","c")

nfiles = dimsizes(files)
do k=0,nfiles-1
 print("k="+k)
 f = addfile(files(k),"r")
 ntime = dimsizes(f->time)
 fld0 = dim_avg_Wrap(f->PSFC(lat|:,lon|:,time|:))
 if(k.eq.0) then
  fld = fld0
  fld = 0.
 end if
 fff = f->PSFC(:,:,:)
 do m=0,ntime-1
  fld = fld + fld0 - fff(m,:,:)
 end do
 print(min(fld)+" "+max(fld))
end do
g->PSFC = fld
