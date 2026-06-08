subroutine pressure

! call a pressure solver 

use vars
use terrain, only: terra, terrau, terrav, terraw, k_terra
use params, only: dolatlon, doflat, donodynamics, npressure_iter
use params, only: dowallx, dowally, doregion, docolumn
use check_energy
implicit none
real(8) p0(nzm),p01(nzm),coef,pmin,pmax
integer i,j,k

if(donodynamics.or.docolumn) return

call t_startf ('pressure')

call sumenergy(1,'pressure:')

! do iteration to reduce residual velocity inside obstacles or terrain
! if required

do k = 0,npressure_iter

  if(k.gt.0) then
   if(k.eq.1) p(:,:,:,nc) = p(:,:,:,na)
   u(1:nx,1:ny,1:nzm) = u(1:nx,1:ny,1:nzm)*terrau(1:nx,1:ny,1:nzm)
   v(1:nx,1:ny,1:nzm) = v(1:nx,1:ny,1:nzm)*terrav(1:nx,1:ny,1:nzm)
   w(1:nx,1:ny,1:nzm) = w(1:nx,1:ny,1:nzm)*terraw(1:nx,1:ny,1:nzm)
  end if

  if(RUN3D) then
  ! if(mod(nx_gl,nsubdomains).ne.0.or.mod(ny_gl,nsubdomains).ne.0) then
   if(RUN2D.or.mod(nx_gl,nsubdomains).ne.0) then
    call pressure_orig()
   else
    if(dolatlon.or.doyvar) then
     if(doflat) then
       call pressure_big(.false.)
     else
       call pressure_gmg(.false.)
     end if
    else
     call pressure_big(.false.)
    end if
   end if
  else
    call pressure_orig()
   ! if(masterproc) print*,'pressure_orig is missing'
   ! call task_abort()
  end if

  if(k.gt.0) p(:,:,:,nc) = p(:,:,:,nc)+p(:,:,:,na)

end do ! k

if(npressure_iter.gt.0) p(:,:,:,na) = p(:,:,:,nc)

! compute global pressure field by subtracting global mean pressure perturbation
! (as presure from solver is defined within addition of an arbitrary constant
! at each level) and then adding to the reference profile:

! compute global pressure perturbation
p0(:) = 0.
do k=1,nzm
 do j=1,ny
  do i=1,nx
   p0(k)=p0(k)+p(i,j,k,na)*rho(k)*wgt(j,k)*terra(i,j,k)
!   p0(k)=p0(k)+p(i,j,k,na)*rho(k)*wgts(j,k)
  end do
 end do
end do

if(dompi) then
  call task_sum_real8(p0,p01,nzm)
  p0(:) = p01(:)
end if
coef = 1./dble(nsubdomains*nx*ny)
p0(:) = p0(:)*coef
! average in the vertical now:

do k=1,nzm
 pmin = 0.85*pres(k) ! put some guard against ocasional "flare-up" 
 pmax = 1.15*pres(k) ! of GMG pressure solver 
 do j=1-YES3D,ny
  do i=0,nx
     pp(i,j,k) = min(pmax,max(pmin,(p(i,j,k,na)*rho(k)-p0(k))*0.01+pres(k)))
  end do
 end do
end do

! compute surface pressure:

do j=1,ny
  do i=1,nx
   k = k_terra(i,j)
   psfc_xy(i,j) = psfc_xy(i,j) + pp(i,j,k)*dtfactor
  end do
end do

if(.not.doregion.and.dowallx.and.mod(rank,nsubdomains_x).eq.0) then

    do k=1,nzm
     do j=1,ny
      dudt(1,j,k,na) = 0.
      u(1,j,k) = 0.
     end do
    end do

end if

if(.not.doregion.and.dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then

    do k=1,nzm
     do i=1,nx
      dvdt(i,1,k,na) = 0.
      v(i,1,k) = 0.
     end do
    end do

end if

call sumenergy(-1,'pressure:')

call t_stopf ('pressure')

end
