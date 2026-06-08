
! ---------------------------------------------------------
! Fill-in the scalar inside the terrain
! assuming zero laplacian

subroutine terrain_fill(niter)
use grid
use terrain, only: terra, kmax
use vars, only: t
use microphysics, only: micro_field, index_water_vapor, nmicro_fields, flag_wmass
use sgs, only: sgs_field
implicit none
integer, intent(in) :: niter
real fld(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm)

real(8) a(ny),b(ny),e(ny)
integer i,j,k,kb,kc,iter,m
real tmp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), coef, coef1

return

call t_startf('terrain_fill')

m = index_water_vapor
where(terra(:,:,:).lt.1.) micro_field(:,:,:,m) = 0.

!do k=1,nzm
! do j=1,ny
!  do i=1,nx
!   if(terra(i,j,k).lt.1.) then
!     micro_field(i,j,k,m) = (terra(i-1,j,k)*micro_field(i-1,j,k,m)+terra(i+1,j,k)*micro_field(i+1,j,k,m)+ &
!                             terra(i,j-1,k)*micro_field(i,j-1,k,m)+terra(i,j+1,k)*micro_field(i,j+1,k,m))/ &
!                             (terra(i-1,j,k)+terra(i+1,j,k)+terra(i,j-1,k)+terra(i,j+1,k)+1.e-5) 
!   end if
!  end do
! end do
!end do
!
!tmp(:,:,:) = t(:,:,:)
!coef = 1./100.
!coef1 = 1./(1.+coef)
!do k=1,nzm
! do j=1,ny
!  do i=1,nx
!   if(terra(i,j,k).lt.1.) then
!     t(i,j,k) = (terra(i-1,j,k)*tmp(i-1,j,k)+terra(i+1,j,k)*tmp(i+1,j,k)+ &
!                 terra(i,j-1,k)*tmp(i,j-1,k)+terra(i,j+1,k)*tmp(i,j+1,k))/ &
!                 (terra(i-1,j,k)+terra(i+1,j,k)+terra(i,j-1,k)+terra(i,j+1,k)+1.e-5)  
!     t(i,j,k) = (t(i,j,k)+coef*tmp(i,j,k))*coef1
!   end if
!  end do
! end do
!end do


!do m=1,nmicro_fields
!   if(m.ne.index_water_vapor) micro_field(:,:,:,index_water_vapor) = &
!         micro_field(:,:,:,index_water_vapor) * terra(:,:,:)
!!   if(m.ne.index_water_vapor) then
!!    do k=1,nzm
!!      do j=1,ny
!!        do i=1,nx
!!          if(terra(i,j,k).lt.1.) then
!!           if(flag_wmass(m).eq.1) micro_field(i,j,k,index_water_vapor) = &
!!             micro_field(i,j,k,index_water_vapor) + micro_field(i,j,k,m)
!!           micro_field(i,j,k,m) = 0.
!!          end if
!!        end do
!!      end do
!!    end do
!!   end if
! end do

!do j=1,ny
!  a(j) = mu(j)*muv(j+1)*(dx/dy)**2/(ady(j)*adyv(j+1))
!  b(j) = mu(j)*muv(j)*(dx/dy)**2/(ady(j)*adyv(j))
!  e(j) = 1./(2.+a(j)+b(j))
!end do
!do iter=1,niter
! if(niter.gt.1) call boundaries(3)
! do k=1,kmax
!  kc = min(k+1,nzm)
!  kb = max(1,k-1)
!  do j=1,ny
!   do i=1,nx
!     if(terra(i,j,k).lt.1.) then
!      t(i,j,k) = (t(i+1,j,k)+t(i-1,j,k)+a(j)*t(i,j+1,k)+b(j)*t(i,j-1,k))*e(j)
!      micro_field(i,j,k,index_water_vapor) =  &
!       (micro_field(i+1,j,k,index_water_vapor)+micro_field(i-1,j,k,index_water_vapor)+ &
!        a(j)*micro_field(i,j+1,k,index_water_vapor)+b(j)*micro_field(i,j-1,k,index_water_vapor))*e(j)
!      sgs_field(i,j,k,:) = 0.
!     end if
!   end do
!  end do
! end do
!end do

call t_stopf('terrain_fill')
end subroutine terrain_fill

