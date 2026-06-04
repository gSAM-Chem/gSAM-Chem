use gmg

implicit none

integer, parameter ::  ny_gl=4608, nzm=74
real(8) yy(ny_gl+1),zi8(nzm+1),pp(nzm+1)
character(len=1) mu_discr,rho_discr
logical unit
real(8):: earth_factor
real(8) pfy(nzm,ny_gl)  
integer  ny_gl1, nzm1
real(8) alpha
real(8) ff(nzm,ny_gl)
integer v_cyc, niter
logical debug
real(8) eps,resnorm,error
integer k,j

open(20,file='../../OUT_3D/out1.bin',form='unformatted')
read(20) ny_gl1,nzm1,yy,zi8,pp,mu_discr,rho_discr,unit,earth_factor
close(20)
print*,'ny_gl, nzm:',ny_gl1,nzm1
print*,'yy:',yy
print*
print*,'zi8:',zi8
print*
print*,'pp:',pp
print*
print*,'mu_discr,rho_discr,unit,earth_factor:',mu_discr,rho_discr,unit,earth_factor

call pressure_gmg_init(ny_gl,nzm,yy,zi8,pp,mu_discr,rho_discr,unit,earth_factor)

open(20,file='../../OUT_3D/out2.bin',form='unformatted')
read(20) alpha,ff,eps,v_cyc,pfy,resnorm,niter,debug
print*,'alpha=',alpha
do j=1,ny_gl
  do k=1,nzm
    if(isnan(ff(k,j))) print*,'ff:',j,k,ff(k,j)
    if(isnan(pfy(k,j))) print*,'pfy',j,k,pfy(k,j)
  end do
end do
print*,'ff(min/max)=',minval(ff),maxval(ff)
print*,'eps=',eps
print*,'v_cyc=',v_cyc
print*,'pfy(min/max)=',minval(pfy),maxval(pfy)
print*,'resnorm=',resnorm
print*,'niter=',niter
print*,'debug=',debug

call pressure_gmg_solv(alpha,ff,eps,v_cyc,pfy(:,:),resnorm,niter,debug)

print*,'solution:'
print*,'pfy(min/max)=',minval(pfy),maxval(pfy)
print*,'resnorm=',resnorm
print*,'niter=',niter


end
