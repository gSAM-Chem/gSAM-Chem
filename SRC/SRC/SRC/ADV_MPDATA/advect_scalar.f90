subroutine advect_scalar (f,fadv,flux,f2leadv,f2legrad,fwleadv,doit,dofix,dosubtr,dodble,fsubtr)
 	
!     positively definite monotonic advection with non-oscillatory option

use grid
use vars, only: u1, v1, w1, rho, rhow, misc, misc1
use params, only: docolumn, doterrain

implicit none

real, intent(inout) :: f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real, intent(out) :: flux(nz), fadv(nz)
real, intent(out) :: f2leadv(nzm),f2legrad(nzm),fwleadv(nzm)
logical, intent(in) ::  doit ! do some extended statistics
logical, intent(in) :: dofix ! enforce no scalar flux through walls
logical, intent(in) :: dosubtr ! subtract iminimum for better precision for real4
logical, intent(in) :: dodble ! advect with double precision
real, intent(in) :: fsubtr

real, allocatable ::  df(:,:,:)
real f0(nzm),df0(nzm),fff(nz),factor
real coef
integer i,j,k

if(docolumn) then
  fadv = 0.
  flux = 0.
  f2leadv = 0.
  f2legrad = 0.
  fwleadv = 0.
  return
end if

if(dostatis) then
	
 allocate (df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))
 df(:,:,:) = f(:,:,:)

endif

if(RUN3D) then
 if(dofix.and.doterrain) then
  call advect_scalar3D_tracer(f, u1, v1, w1, rho, rhow, flux)
 else
  if(dodble) then
    call advect_scalar3D_dble(f, u1, v1, w1, rho, rhow, flux)
  else
    call advect_scalar3D(f, u1, v1, w1, rho, rhow, flux, dosubtr, fsubtr)
  end if
 end if
else
  call advect_scalar2D(f, u1, w1, rho, rhow, flux)	  
endif

if(dostatis) then

  do k=1,nzm
    fadv(k)=0.
    do j=1,ny
     do i=1,nx
      fadv(k)=fadv(k)+f(i,j,k)-df(i,j,k)
     end do
    end do
  end do

end if

if(dostatis.and.doit) then
	
  misc1(:,:,:) = misc(:,:,:) 

  call stat_varscalar(f,df,f0,df0,f2leadv)
  call stat_sw2(f,df,fwleadv)


!  Compute advection flux of variance
 

  do k=1,nzm
    do j=dimy1_s,dimy2_s
     do i=dimx1_s,dimx2_s
      df(i,j,k) = (df(i,j,k)-df0(k))**2
     end do
    end do
  end do

  coef = max(1.e-10,maxval(df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, 1:nzm)))
  df(:,:,:) = df(:,:,:) / coef
  if(RUN3D) then
   call advect_scalar3D(df, u1, v1, w1, rho, rhow, fff, dosubtr, fsubtr)
  else
   call advect_scalar2D(df, u1, w1, rho, rhow, fff)	  
  endif
  df(:,:,:) = df(:,:,:) * coef

  factor=dz/(nx*ny*dtn)
  do k = 1,nzm
    fff(k)=fff(k) * factor
  end do
  fff(nz)=0.
  do k = 1,nzm
    f2legrad(k) = f2leadv(k)
    f2leadv(k)=-(fff(k+1)-fff(k))/(dz*adz(k)*rho(k))	 
    f2legrad(k)=f2legrad(k)-f2leadv(k)
  end do

  misc(:,:,:) = misc1(:,:,:)

endif

if(dostatis) deallocate(df)

end subroutine advect_scalar

