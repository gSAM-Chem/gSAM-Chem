
! A version with non-blocking receive and blocking send. Doesn't need 
! EAGER_LIMIT set for communication.

subroutine pressure_gmg(doflag)
	
!   parallel iterative pressure-solver for 3D large domains.
!   with multigrid solver in y-z
!   (C) 2015, Marat Khairoutdinov

! CHANGE: Jan 2025: implemented non-blocking send for 4D bufy1 instead of 3D 
!                   modified subroutines: transpose_y and transpose_y_inv - MK

use vars
use mpi_stuff
use params, only: dowallx, dowally, dolatlon, ggr, &
                  earth_factor, gmg_precision, doflat, gamma_RAVE
use gmg, only: gmg_init, gmg_solve
implicit none

logical, intent(in) :: doflag ! switch-flag; if .true. then not call certain subroutines and
                              ! handle special cases for doregion=.true.

integer, parameter :: nslab = max(1,nsubdomains/ny_gl) !number of slabs in the vertical
integer, parameter :: nx_s=nx_gl/nsubdomains ! width of the x-slabs
integer, parameter :: ny_s=nslab*ny_gl/nsubdomains ! width of the y-slabs
integer, parameter :: nzm2 = nzm/nslab

! Slabs:
real fx(nx_gl, ny_s, nzm2)! slab for x-pass Fourier coefs
real fy(nzm, ny_gl, nx_s) ! slab for GMG sover (y-z slabs)
real(8) gx(nx_gl+2, ny_s) ! array to perform FFT in x

real bufx1(nx, ny_s, nzm2, max(1,nsubdomains_x))
real bufx2(nx, ny_s, nzm2, max(1,nsubdomains_x))
real bufy1(nx_s, ny, nzm, max(1,nsubdomains_y))
real bufy2(nx_s, ny, nzm, max(1,nsubdomains_y))

! save memory by collocating some of the large arrays
equivalence (bufx1(1,1,1,1),bufy1(1,1,1,1))
equivalence (bufx2(1,1,1,1),bufy2(1,1,1,1))

! rhs and solution (without extra boundaries)
real ppp(nx,ny,nzm)

! FFT stuff:
real(8) work(max((nx_gl+3)*(ny_s+1),(nx_s+1)*(ny_gl+2)))
real(8) trigxi(3*nx_gl/2+1)
integer ifaxi(100)

! iterative solver stuff:
real(8) alpha
real(8) yy(ny_gl+1)
real(8) pz(nz)
real(8) zi8(nz)
real(8) rho8(nzm)
real(8) rad_reduc_factor
real(8), allocatable ::  pfy0(:,:,:)
integer niter,v_cyc,ierr
logical unit
character(len=1) mu_discr,rho_discr

! Misc
real(8) ff(nzm,ny_gl)
real(8) gg(ny_gl+1,nzm+1)
real(8) xi,xnx,ddx,pii,factx
real(8) eps,resnorm,error
integer reqs_in(nsubdomains), ranks(nsubdomains), tags(nsubdomains)
integer reqs_out(nsubdomains)
integer i, j, k, id, m, n, it, jt, tag, ii,jj,jb,jc,kb,kc
integer irank, rnk
integer n_in, n_out, counter
logical flag(nsubdomains)
integer jwall

! for wrapping p
real buff_ew1(ny,nzm), buff_ns1(nx,nzm)
real buff_ew2(ny,nzm), buff_ns2(nx,nzm)
integer rf, tagrf
logical waitflag
character(4) rankchar

call t_startf('pressure_gmg')

! Make sure that the grid is suitable for the solver:

if(.not.dowally) then
  if(masterproc) print*,'pressure_gmg: iterative solver requires a wall in y (dowall=.true.). STOP'
  call task_abort
end if
if(mod(nx_gl,nsubdomains).ne.0) then
  if(masterproc) print*,'pressure_gmg: nx_gl/nsubdomains is not round number. STOP'
  call task_abort
endif
if(mod(nzm,nslab).ne.0) then
  if(masterproc) print*,'pressure_gmg: nzm is not divisible by nslab. STOP'
  call task_abort
endif
!if(mod(ny_gl,nsubdomains).ne.0) then
!  if(masterproc) print*,'pressure_gmg: ny_gl/nsubdomains is not round number. STOP'
!  call task_abort
!endif

!-----------------------------------------------------------------

if(RUN2D) then
  print*,'pressure3D() cannot be called for 2D domains. Quitting...'
  call task_abort()
endif

!-----------------------------------------------------------------
if(.not.allocated(pfy)) then
  allocate(pfy(nzm,ny_gl,nx_gl/nsubdomains))
  pfy = 0.
end if
if(doflag) then
  allocate(pfy0(nzm,ny_gl,nx_gl/nsubdomains))
  pfy0 = pfy
  pfy = 0.
end if
!==========================================================================
!  Compute the r.h.s. of the Poisson equation for pressure

call t_startf ('press_rhs')

if(doflag) then
  ppp(1:nx,1:ny,1:nzm) = p(1:nx,1:ny,1:nzm,na)
else
  call press_rhs(ppp)
end if

call t_stopf ('press_rhs')
! variable p will also be used as grid decomposition placeholder between transposes
!   for the fourier coefficients

!==========================================================================
!   Form the vertical slabs (x-z) of right-hand-sides of Poisson equation 
!   for the FFT - one slab per a processor.


call t_startf ('transpose_x')
call transpose_x(ppp)
call t_stopf ('transpose_x')

!==========================================================================
! Perform Fourier transformation n x-direction for a slab:
! This routing produces  N+2 coefficients for input array of langth N.
! However, two of those (2, and N+2) are identical zeros, so there is no point to
! send them around during transpose. That's why the temporary array gx is used to
! conform to requirements of fft991_crm subroutine.

call t_startf ('press_fft')
call fftfax_crm(nx_gl,ifaxi,trigxi)
if(dowallx) then
  do k=1,nzm2
    gx(1:nx_gl,1:ny_s) = fx(1:nx_gl,1:ny_s,k)
    call cosft_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,-1)
    fx(1:nx_gl,1:ny_s,k) = gx(1:nx_gl,1:ny_s)
  end do
else
  do k=1,nzm2
    gx(1:nx_gl,1:ny_s) = fx(1:nx_gl,1:ny_s,k)
    call fft991_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,-1)
    fx(1,1:ny_s,k) = gx(1,1:ny_s)
    fx(2:nx_gl,1:ny_s,k) = gx(3:nx_gl+1,1:ny_s)
  end do
end if
call task_barrier()
call t_stopf ('press_fft')

!=========t=================================================================
!   Form blocks again from x-z slabs   

call t_startf ('transpose_x_inv')
call transpose_x_inv(ppp)
call t_stopf ('transpose_x_inv')

!==========================================================================
!   Form the vertical slabs (y-z) of Fourier coefs  
!   for the FFT - in y, one slab per a processor.

call t_startf ('transpose_y')
call transpose_y()
call t_stopf ('transpose_y')

call t_startf ('press_solver')
!==========================================================================
! solve 2D elliptic problem by iterations:

if(dolatlon) then
 if(doflat) then
   unit = .true.
 else
   unit = .false.
 end if
else
 unit = .true.
end if
v_cyc=100
mu_discr='B'
rho_discr='B'
error_max = 0.
niter_max = 0
pii = acos(-1._8)
xnx=nx_gl
it=rank*nx_s
pz(1)=presi(1)*100._8/ggr
zi8(1) = zi(1)*gamma_RAVE
do k=1,nzm
 rho8(k) = rho(k)
 zi8(k) = zi(k)*gamma_RAVE
end do
zi8(nz) = zi(nz)*gamma_RAVE
do k=2,nz
 pz(k) = pz(k-1)-rho8(k-1)*(zi8(k)-zi8(k-1))
end do
do j=1,ny_gl+1
 yy(j) = yv_gl(j)
! yy(j) = (j-ny_gl/2-1.)*dy
end do
ddx=1._8/(dx*dx)
if(doflag.or.nstep.le.3) then
  eps = 1.e-10
else
 eps = gmg_precision
end if
if(.not.gmg_initialized) then   
  rad_reduc_factor = dble(earth_factor)
!  if(masterproc) write(1,*) ny_gl,nzm,yy(:),zi8(:),pz(:),mu_discr,rho_discr,unit,rad_reduc_factor
  call gmg_init(ny_gl,nzm,yy,zi8,pz,mu_discr,rho_discr,unit,rad_reduc_factor)
  gmg_initialized = .true.
end if
do i=1,nx_s
    if(dowallx) then
      id=i+it-1
      factx = 1._8
    else
      id=(i+1+it-0.1)/2.
      factx = 2._8
    end if
    xi=id
    do k=1,nzm
     ff(k,1:ny_gl)=fy(k,1:ny_gl,i)*mu_gl(1:ny_gl)**2
    end do
    alpha = (2._8-2._8*cos(factx*pii/xnx*xi))*ddx
    call gmg_solve(alpha,ff,eps,v_cyc,pfy(:,:,i),resnorm,niter,ierr)
!    if(alpha.eq.0.) write(1,*)alpha,ff(:,:),eps,v_cyc,pfy(:,:,i)
!    write(*,'(3I5,2G20.5)')rank,niter,ierr,alpha,resnorm
    niter_max = max(niter_max, niter)
    error_max = max(error_max, real(resnorm))
end do

call task_barrier()
call t_stopf ('press_solver')
!==========================================================================
!   Form blocks again from y-z slabs   

call t_startf ('transpose_y_inv')
call transpose_y_inv()
call t_stopf ('transpose_y_inv')

!==========================================================================
!   Form the vertical slabs (x-z) of Fourier coefs
!   for the inverse FFT - in x, one slab per a processor.

 call t_startf ('transpose_x')
 call transpose_x(ppp)
 call t_stopf ('transpose_x')

! Perform inverse Fourier transform n x-direction for a slab:

 call t_startf ('press_fft')
 call fftfax_crm(nx_gl,ifaxi,trigxi)

if(dowallx) then
 do k=1,nzm2
  gx(1:nx_gl,1:ny_s) = fx(1:nx_gl,1:ny_s,k)
  gx(nx_gl+1:nx_gl+2,:) = 0.
  call cosft_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,1)
  fx(1:nx_gl,1:ny_s,k) = gx(1:nx_gl,1:ny_s)
 end do
else
 do k=1,nzm2
  gx(1,1:ny_s) = fx(1,1:ny_s,k)
  gx(2,:) = 0.
  gx(3:nx_gl+1,1:ny_s) = fx(2:nx_gl,1:ny_s,k)
  gx(nx_gl+2,:) = 0.
  call fft991_crm(gx,work,trigxi,ifaxi,1,nx_gl+2,nx_gl,ny_s,1)
  fx(1:nx_gl,1:ny_s,k) = gx(1:nx_gl,1:ny_s)
 end do
end if

call task_barrier()
call t_stopf ('press_fft')

call t_startf ('transpose_x_inv')
call transpose_x_inv(ppp)
call t_stopf ('transpose_x_inv')

call t_startf ('press_grad')
!==========================================================================
!  Update the pressure fields in the subdomains
!
p(1:nx,1:ny,1:nzm,na) = ppp(:,:,:)

! when we cut back on the ffts wrap p here - look to sib for the model
!DD temporary measure for dompi
  if(dompi) then
!DD custom build a wrap that sends from the north and east edges, and receives at the
!DD    south and west edges

      if(rank==rankee) then
        p(0,1:ny,1:nzm,na) = p(nx,1:ny,1:nzm,na)
      else
        call task_receive_float(buff_ew2(:,:),ny*nzm,reqs_in(1))
          buff_ew1(1:ny,1:nzm) = p(nx,1:ny,1:nzm,na)
          call task_bsend_float(rankee,buff_ew1(:,:),ny*nzm,1)
          waitflag = .false.
        do while (.not.waitflag)
            call task_test(reqs_in(1),waitflag,rf,tagrf)
        end do
        call task_barrier()
        p(0,1:ny,1:nzm,na) = buff_ew2(1:ny,1:nzm)
      endif

      if(rank==ranknn) then
         p(:,1-YES3D,1:nzm,na) = p(:,ny,1:nzm,na)
      else
             call task_receive_float(buff_ns2(:,:),nx*nzm,reqs_in(1))
             buff_ns1(1:nx,1:nzm) = p(1:nx,ny,1:nzm,na)
             call task_bsend_float(ranknn,buff_ns1(:,:),nx*nzm,1)
             waitflag = .false.
         do while (.not.waitflag)
               call task_test(reqs_in(1),waitflag,rf,tagrf)
             end do
         call task_barrier()
         p(1:nx,1-YES3D,1:nzm,na) = buff_ns2(1:nx,1:nzm)
      endif

  else
    p(0,:,1:nzm,na) = p(nx,:,1:nzm,na)
    p(:,1-YES3D,1:nzm,na) = p(:,ny,1:nzm,na)
  endif
!DD end ugly wrap code.

  ! overwrite for the case of walls in x 
  call task_rank_to_index(rank,it,jt)
  if(dowallx.and.it.eq.0) then
    p(0,:,1:nzm,na) = p(1,:,1:nzm,na)
  end if
  if(dowally.and.jt.eq.0) then
    p(:,1-YES3D,1:nzm,na) = p(:,1,1:nzm,na)
  end if


!==========================================================================
!  Ad2d pressure gradient term to the rhs of the momentum equation:

if(.not.doflag) call press_grad()

call t_stopf ('press_grad')

if(doflag) then
  pfy = pfy0
  deallocate(pfy0)
end if

call t_stopf('pressure_gmg')
!==========================================================================
!==========================================================================
!==========================================================================

contains

!==========================================================================
   subroutine transpose_x(pm)

! transpose from blocks to x-z slabs

      REAL, INTENT(in) :: pm(nx, nslab*ny, nzm2)
      
      irank = rank-mod(rank,nsubdomains_x)  

      n_in = 0
      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufx2(1:nx,1:ny_s,1:nzm2,n_in),nx*ny_s*nzm2,reqs_in(n_in))
          flag(n_in) = .false.

        end if

      end do ! m

      n_out = 0
      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          n = m-irank

          n_out = n_out + 1
          bufx1(1:nx,1:ny_s,1:nzm2,n_out) = pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2)
          call task_send_float(m,bufx1(1:nx,1:ny_s,1:nzm2,n_out),nx*ny_s*nzm2, 33, reqs_out(n_out)) 

        endif

      end do ! m


! don't sent a buffer to itself, just fill directly.

      n = rank-irank
      call task_rank_to_index(rank,it,jt)
      fx(1+it:nx+it,1:ny_s,1:nzm2) = pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2)


      ! Fill slabs when receive buffers are full:

      counter = n_in

      do while (counter .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
      	    call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then 
            	 counter=counter-1
                 call task_rank_to_index(rnk,it,jt)	  
                 fx(1+it:nx+it,1:ny_s,1:nzm2) = bufx2(1:nx,1:ny_s,1:nzm2,m)
              endif   
          endif
         end do
      end do

      call task_waitall(n_out,reqs_out,ranks,tags)

!      call task_waitall(counter,reqs_in,ranks,tags)
!      do m = 1,n_in
!         call task_rank_to_index(ranks(m),it,jt)	  
!         fx(1+it:nx+it,1:ny_s,1:nzm2) = bufx2(1:nx,1:ny_s,1:nzm2,m)
!      end do

      call task_barrier()

   end subroutine transpose_x
   
!==========================================================================
   subroutine transpose_x_inv(pm)

! transpose from x-z slabs to blocks

      REAL, INTENT(out) :: pm(nx, nslab*ny, nzm2)
      
      irank = rank-mod(rank,nsubdomains_x)
      n_in = 0
      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufx2(1:nx,1:ny_s,1:nzm2,n_in),nx*ny_s*nzm2,reqs_in(n_in))
          flag(n_in) = .false.

        endif

      end do ! m

      n_out = 0
      do m = irank, irank+nsubdomains_x-1

        if(m.ne.rank) then

          call task_rank_to_index(m,it,jt)
          n_out = n_out + 1
          bufx1(1:nx,1:ny_s,1:nzm2,n_out) = fx(1+it:it+nx,1:ny_s,1:nzm2)
          call task_send_float(m,bufx1(1:nx,1:ny_s,1:nzm2,n_out),nx*ny_s*nzm2, 33, reqs_out(n_out))

        endif

      end do ! m

! don't sent a buffer to itself, just fill directly.

      n = rank-irank
      call task_rank_to_index(rank,it,jt)
      pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2) = fx(1+it:nx+it,1:ny_s,1:nzm2)  

! Fill slabs when receive buffers are full:

      counter = n_in

      do while (counter .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
              call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then
                 counter=counter-1
                 n = rnk-irank
                 pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2) = bufx2(1:nx,1:ny_s,1:nzm2,m)
              endif
         endif
        end do
      end do

      call task_waitall(n_out,reqs_out,ranks,tags)

!      call task_waitall(counter,reqs_in,ranks,tags)
!      do m = 1,n_in
!         n = ranks(m)-irank
!         pm(1:nx,n*ny_s+1:n*ny_s+ny_s,1:nzm2) = bufx2(1:nx,1:ny_s,1:nzm2,m)
!      end do

      call task_barrier()
   end subroutine transpose_x_inv

!==========================================================================
   subroutine transpose_y()

! transpose from blocks to y-z slabs

      irank = rank / nsubdomains_y

      call task_rank_to_index(rank,it,jt)
      n_in = 0
      do m = irank, nsubdomains-1, nsubdomains_x

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufy2(1:nx_s,1:ny,1:nzm,n_in),ny*nx_s*nzm,reqs_in(n_in))
          flag(n_in) = .false.

        else
! don't sent a buffer to itself, just fill directly.

          n = mod(rank,nsubdomains_y)
          do i = 1,nx_s
            fy(1:nzm,1+jt:ny+jt,i) = transpose(ppp(n*nx_s+i,1:ny,1:nzm))
          enddo

        end if

      end do ! m

      irank = nsubdomains_y*mod(rank,nsubdomains_x)
      n_out = 0
      do m = irank, irank+nsubdomains_y-1

        if(m.ne.rank) then

          n = m-irank

          n_out = n_out + 1
          bufy1(1:nx_s,1:ny,1:nzm,n_out) = ppp(n*nx_s+1:n*nx_s+nx_s,1:ny,1:nzm)
          call task_send_float(m,bufy1(1:nx_s,1:ny,1:nzm,n_out),ny*nx_s*nzm, 33, reqs_out(n_out))

        endif

      end do ! m


      ! Fill slabs when receive buffers are full:

      counter = n_in

      do while (counter .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
            call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then
                 counter=counter-1
                 call task_rank_to_index(rnk,it,jt)
                 do i = 1,nx_s
                   fy(1:nzm,1+jt:ny+jt,i) = transpose(bufy2(i,1:ny,1:nzm,m))
                 enddo
              endif
          endif
         end do
      end do

      call task_waitall(n_out,reqs_out,ranks,tags)


!      call task_waitall(counter,reqs_in,ranks,tags)
!      do m = 1,n_in
!         call task_rank_to_index(ranks(m),it,jt)
!         do i = 1,nx_s
!           fy(1:nzm,1+jt:ny+jt,i) = transpose(bufy2(i,1:ny,1:nzm,m))
!         enddo
!      end do


      call task_barrier()

   end subroutine transpose_y

!==========================================================================
   subroutine transpose_y_inv()

! transpose from y-z slabs to blocks

      n_in = 0
      irank = nsubdomains_y*mod(rank,nsubdomains_x)
      call task_rank_to_index(rank,it,jt)
      do m = irank, irank+nsubdomains_y-1

        if(m.ne.rank) then

          n_in = n_in + 1
          call task_receive_float(bufy2(1:nx_s,1:ny,1:nzm,n_in),ny*nx_s*nzm,reqs_in(n_in))
          flag(n_in) = .false.

        else

! don't sent a buffer to itself, just fill directly.

          n = rank-irank
          do i = 1,nx_s
            ppp(n*nx_s+i,1:ny,1:nzm) = transpose(pfy(1:nzm,1+jt:ny+jt,i))
          enddo

        endif

      end do ! m

      irank = rank / nsubdomains_y
      n_out = 0
      do m = irank, nsubdomains-1, nsubdomains_x

        if(m.ne.rank) then

          call task_rank_to_index(m,it,jt)
          n_out = n_out + 1
          do i = 1,nx_s
            bufy1(i,1:ny,1:nzm,n_out) = transpose(pfy(1:nzm,1+jt:jt+ny,i))
          enddo
          call task_send_float(m,bufy1(1:nx_s,1:ny,1:nzm,n_out),ny*nx_s*nzm, 33, reqs_out(n_out))

        endif

      end do ! m

! Fill slabs when receive buffers are full:

      irank = nsubdomains_y*mod(rank,nsubdomains_x)
      counter = n_in

      do while (counter .gt. 0)
        do m = 1,n_in
         if(.not.flag(m)) then
              call task_test(reqs_in(m), flag(m), rnk, tag)
              if(flag(m)) then
                 counter=counter-1
                 n = rnk-irank
                 ppp(n*nx_s+1:n*nx_s+nx_s,1:ny,1:nzm) = bufy2(1:nx_s,1:ny,1:nzm,m)
              endif
         endif
        end do
      end do

      call task_waitall(n_out,reqs_out,ranks,tags)

!      call task_waitall(counter,reqs_in,ranks,tags)
!      do m = 1,counter
!         n = ranks(m)-irank
!         ppp(n*nx_s+1:n*nx_s+nx_s,1:ny,1:nzm) = bufy2(1:nx_s,1:ny,1:nzm,m)
!      end do


      call task_barrier()

   end subroutine transpose_y_inv

end 

