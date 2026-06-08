module sgs

! module for original SAM subgrid-scale SGS closure (Smagorinsky or 1st-order TKE)
! Marat Khairoutdinov, 2012

use grid, only: nx,nxp1,ny,nyp1,YES3D,nzm,nz,dimx1_s,dimx2_s,dimy1_s,dimy2_s 
use params, only: dosgs
implicit none

!----------------------------------------------------------------------
! Required definitions:

!!! prognostic scalar (need to be advected arround the grid):

integer, parameter :: nsgs_fields = 1   ! total number of prognostic sgs vars

real sgs_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nsgs_fields)

!!! sgs diagnostic variables that need to exchange boundary information (via MPI):

integer, parameter :: nsgs_fields_diag = 2   ! total number of diagnostic sgs vars

! diagnostic fields' boundaries:
integer, parameter :: dimx1_d=0, dimx2_d=nxp1, dimy1_d=1-YES3D, dimy2_d=nyp1

real sgs_field_diag(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm, nsgs_fields_diag)

logical:: advect_sgs = .false. ! advect prognostics or not, default - not (Smagorinsky)
logical, parameter:: do_sgsdiag_bound = .true.  ! exchange boundaries for diagnostics fields

! SGS fields that output by default (if =1).
integer, parameter :: flag_sgs3Dout(nsgs_fields) = (/0/)
integer, parameter :: flag_sgsdiag3Dout(nsgs_fields_diag) = (/0,0/)

real fluxbsgs (nx, ny, 1:nsgs_fields) ! surface fluxes 
real fluxtsgs (nx, ny, 1:nsgs_fields) ! top boundary fluxes 

!!! these arrays may be needed for output statistics:

real sgswle(nz,1:nsgs_fields)  ! resolved vertical flux
real sgswsb(nz,1:nsgs_fields)  ! SGS vertical flux
real sgsadv(nz,1:nsgs_fields)  ! tendency due to vertical advection
real sgslsadv(nz,1:nsgs_fields)  ! tendency due to large-scale vertical advection
real sgsdiff(nz,1:nsgs_fields)  ! tendency due to vertical diffusion

real tkmax(0:ny,nzm)  ! maximum allowed coef for SGS diffusion (scalar)
real tkmin(0:ny,nzm)  ! minimum coef for SGS diffusion (scalar)
real tkmaxu(0:ny,nzm)  ! maximum allowed coef for SGS diffusion for w
real tkminu(0:ny,nzm)  ! minimum coef for SGS diffusion for u
real tkmaxv(0:ny,nzm)  ! maximum allowed coef for SGS diffusion for w
real tkminv(0:ny,nzm)  ! minimum coef for SGS diffusion for v
real tkmaxw(0:ny,nzm)  ! maximum allowed coef for SGS diffusion for w
real tkminw(0:ny,nzm)  ! minimum coef for SGS diffusion for w

real, parameter :: cfl_diff_max = 0.46 ! maximum CFL allowed to run for diffusion of momentum
real, parameter :: cfl_diffsc_max = 0.46  ! maximum CFL allowed to run for diffusion of scalars

!------------------------------------------------------------------
! internal (optional) definitions:

! make aliases for prognostic variables:

real tke(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! SGS TKE
equivalence (tke(dimx1_s,dimy1_s,1),sgs_field(dimx1_s,dimy1_s,1,1))

! make aliases for diagnostic variables:

real tk  (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy viscosity
real tkh (dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
equivalence (tk(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,1))
equivalence (tkh(dimx1_d,dimy1_d,1), sgs_field_diag(dimx1_d, dimy1_d,1,2))

!-------------------------------------------------------
! namelist variables:
logical:: dosmagor = .true.  ! if true, then use Smagorinsky closure
real:: Cs = 0.19 ! Smagorinsky constant; for CRM Cs=0.22 and for LES Cs=0.15 are recommended - MK 2026
logical:: dodns = .false.   ! if true, run SGS with constant eddy coefficients
real:: RHO_DNS = 0.         ! rederence density for DNS
real:: PRES_DNS = 0.        ! rederence surface pressure for DNS
real:: PR_DNS = 0.71 ! Prandtle number for DNS (ratio of viscosoty to diffusiveity coefs)
real:: DIFF_DNS = 1.5e-5 ! Kinematic diffusivity for FNS (m2/s) 
real:: tk_factor = 0. ! factor to compute minimum tk from maximum allowed for momentum (in tke_full.f90)
logical:: dochamber = .false. ! do a chamber simulation (solid wallss and top and bottom
                              ! kept a prescribed temperatures)
real:: t_top =0.              ! chamber top temperature
real:: t_bot =0.              ! chamber bottom temperature
real:: t_wall =0.             ! chamber wall temperature
logical:: dosidewalls = .false. ! do sidewalls at diff temperature: two at t_bot and two and t_top 

real:: delta_max = 1000. ! max hor. grid spacing in length-scale computation in tke_full
                         ! affects eddy diffusivities.

logical:: doaddTKEdiss = .false. ! add SGS tke dissipation as heating to t() - MK 2025

! Local diagnostics:

real:: tkesbbuoy(nz)=0., tkesbshear(nz)=0., tkesbdiss(nz)=0. , tkesbdiff(nz)=0.

CONTAINS

subroutine sgs_setparm()

  use grid, only: case, caseid, masterproc
  use params, only: dosurface, dowallx, dowally
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder

  !======================================================================
  NAMELIST /SGS_TKE/ dosmagor, tk_factor, dodns, PR_DNS, DIFF_DNS, RHO_DNS, PRES_DNS, &
                     dochamber, t_top, t_bot, t_wall, dosidewalls, delta_max, doaddTKEdiss 

  NAMELIST /BNCUIODSBJCB/ place_holder

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./CASES/'//trim(case)//'/prm', status='old',form='formatted')

  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  read (55,SGS_TKE,IOSTAT=ios)

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in SGS_TKE namelist'
        rewind(55)
        read (55,SGS_TKE)
        call task_abort()
     end if
  end if
  close(55)
   ! write namelist values out to file for documentation
   if(masterproc) then
      open(unit=55,file='./OUT_STAT/'//trim(case)//'_'//trim(caseid)//'.nml', form='formatted', position='append')
      write (unit=55,nml=SGS_TKE,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
  end if

  if(dochamber) then
      dodns = .true.
      dowallx = .true.
      dowally = .true.
      if(t_top.eq.0..or.t_bot.eq.0..or.t_wall.eq.0.) then
       if(masterproc) print*,'for dochamber=true, t_wall, t_bot, and t_top should be set. exit...'
       call task_abort()
      end if
  end if
  if(dodns) then
      dosmagor = .true. ! make sure that no SGS TKE equation is used
      dosurface = .false. ! don't use Monin-Obukhov for surface fluxes. 
      if(t_top.eq.0..or.t_bot.eq.0.) then
       if(masterproc) print*,'for dodns=true, t_bot, and t_top should be set. exit...'
       call task_abort()
      end if
      if(RHO_DNS.eq.0..or.PRES_DNS.eq.0.) then
       if(masterproc) print*,'for dodns=true, RHO_DNS and PRES_DNS should be set. exit...'
       call task_abort()
      end if
  end if
  advect_sgs = .not.dosmagor

  !======================================================================

end subroutine sgs_setparm

!----------------------------------------------------------------------
!!! Initialize sgs:


subroutine sgs_init()

   use grid, only: nrestart,nx_gl, ny_gl,dt,dx,dy,dz, &
                   adz,ady,adyv,adzw,mu,muv,masterproc
   use mpi_stuff, only: rank
   use terrain, only: terra, terrau, terrav, terraw
   integer k, j, it, jt
   real tmp

  if(nrestart.eq.0) then

     sgs_field = 0.
     sgs_field_diag = 0.

     fluxbsgs = 0.
     fluxtsgs = 0.

     if(dodns) then
      tkh(:,:,:) = DIFF_DNS 
      tk(:,:,:) = PR_DNS*DIFF_DNS 
      ! for chamber, make interior of lateral walls as if inside terrain
      ! which is used to correct molecular fluxes at them as only half grid step
      ! are used to compute molecular flux at the walls - MK
      if(dochamber) then
       call task_rank_to_index(rank,it,jt)
       if(it.eq.0) then
         terra(:0,:,:) = 0.
         terrav(:0,:,:) = 0.
         terraw(:0,:,:) = 0.
       end if
       if(it+nx.eq.nx_gl) then
         terra(nx+1:,:,:) = 0.
         terrav(nx+1:,:,:) = 0.
         terraw(nx+1:,:,:) = 0.
       end if
       if(jt.eq.0) then
         terra(:,:0,:) = 0.
         terrau(:,:0,:) = 0.
         terraw(:,:0,:) = 0.
       end if
       if(jt+ny.eq.ny_gl) then
         terra(:,ny+1:,:) = 0.
         terrau(:,ny+1:,:) = 0.
         terraw(:,ny+1:,:) = 0.
       end if
      end if
     end if
 
  end if

  if(masterproc) then
     if(dodns) then
        write(*,*) 'Direct Numerical Simulation (DNS)'
        write(*,*) 'PR_DNS = ', PR_DNS
        write(*,*) 'DIFF_DNS = ', DIFF_DNS,'m2/s'
        write(*,*) 'T_top = ', t_top
        write(*,*) 'T_bot = ', t_bot
       write(*,*) 'Air Density RHO_DNS=',RHO_DNS 
       write(*,*) 'Pressure PRES_DNS=',PRES_DNS 
        tmp = cfl_diffsc_max/(1./dx**2+1./dy**2+1./dz**2)/DIFF_DNS ! maximum timestep for stability
        write(*,*) 'maximum allowed time step for DNS:',tmp
        if(dt.gt.tmp) then
         if(masterproc) print*,'dt is greater than maximum allowed time step for DNS. Stop..'
         call task_abort()
        end if
      if(dochamber) then
       write(*,*) 'Simulation of A chamber'
       write(*,*) 'wall temperature t_wall=', t_wall
       if(.not.dosidewalls) write(*,*) 'T_wall = ', t_wall
       if(dosidewalls) write(*,*) 'dosidewalls = ', dosidewalls
      end if
     else 
      print*,'dosmagor=',dosmagor
      print*,'tk_factor=',tk_factor
      print*,'delta_max=',delta_max
      print*,'doaddTKEdiss=',doaddTKEdiss
      if(dosmagor) then
        write(*,*) 'Smagorinsky SGS Closure'
      else
        write(*,*) 'Prognostic TKE 1.5-order SGS Closure'
      end if
     end if
     print*,'cfl_diff_max=',cfl_diff_max
     print*,'cfl_diffsc_max=',cfl_diffsc_max
  end if

  sgswle = 0.
  sgswsb = 0.
  sgsadv = 0.
  sgsdiff = 0.
  sgslsadv = 0.

  tkmax = 0.
  tkmin = 0.
  tkmaxu = 0. 
  tkminu = 0.
  tkmaxv = 0.
  tkminv = 0.
  tkmaxw = 0.
  tkminw = 0.

end subroutine sgs_init

!----------------------------------------------------------------------
!!! make some initial noise in sgs:
!
subroutine setperturb_sgs(ptype)

use vars, only: q0, z
use terrain
integer, intent(in) :: ptype
integer i,j,k

select case (ptype)

  case(-2)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(k.lt.10) then
           tke(i,j,k)=0.1
         end if
       end do
      end do
     end do


  case(0)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(k.le.4.and..not.dosmagor) then
            tke(i,j,k)=0.04*(5-k)
         endif
       end do
      end do
     end do

  case(1)

     do k=1,nzm
      do j=1,ny
       do i=1,nx
        ! if(q0(k).gt.6.e-3.and..not.dosmagor) then
         if(k.ge.k_terra(i,j).and.k.le.k_terra(i,j)+5) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do

  case(2)

  case(3)   ! gcss wg1 smoke-cloud case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.0.5e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case(4)  ! gcss wg1 arm case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.150..and..not.dosmagor) then
            tke(i,j,k)=0.15*(1.-z(k)/150.)
         endif
       end do
      end do
     end do


  case(5)  ! gcss wg1 BOMEX case

     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(z(k).le.3000..and..not.dosmagor) then
            tke(i,j,k)=1.-z(k)/3000.
         endif
       end do
      end do
     end do

  case(6)  ! GCSS Lagragngian ASTEX


     do k=1,nzm
      do j=1,ny
       do i=1,nx
         if(q0(k).gt.6.e-3.and..not.dosmagor) then
            tke(i,j,k)=1.
         endif
       end do
      end do
     end do


  case default

end select

end subroutine setperturb_sgs

!----------------------------------------------------------------------
!!! Estimate Courant number limit for SGS
!

subroutine kurant_sgs(cfl)

use grid, only: dtn, dx, dy, dz, adz, adzw, imu, ady, mu 
use params, only: doimplicitdiff
implicit none
real, intent(out) :: cfl
real cx, cy, cz, coef
integer j,k,i

if(doimplicitdiff) then
 coef = 0.
else
 coef = 1.
end if
cfl = 0.
do k = 1,nzm
 cz = coef/(dz*adz(k))**2
 do j = 1,ny
  cx = 1./(dx*mu(j))**2
  cy = 1./(dy*ady(j))**2
  do i = 1,nx
    cfl = max(cfl,dtn*min(tkmax(j,k),tkh(i,j,k))*(cx+cy+cz))
  end do
 end do
end do

end subroutine kurant_sgs


!----------------------------------------------------------------------
!!! compute sgs diffusion of momentum:
!
subroutine sgs_mom()

   use vars
   use terrain, only: k_terrau, k_terrav

   real rdz, rdz1
   integer i,j,k

   if(dodns) then

     rdz = 2./dz
     rdz1 = 2./(dz*adzw(nz))
     do j=1,ny
      do i=1,nx
       k = k_terrau(i,j)
       fluxbu(i,j) = -rdz*DIFF_DNS*PR_DNS*u(i,j,k)
       fluxtu(i,j) =  rdz1*DIFF_DNS*PR_DNS*u(i,j,nzm)
      end do
     end do
     do j=1,ny
      do i=1,nx
       k = k_terrav(i,j)
       fluxbv(i,j) = -rdz*DIFF_DNS*PR_DNS*v(i,j,k)
       fluxtv(i,j) =  rdz1*DIFF_DNS*PR_DNS*v(i,j,nzm)
      end do
     end do
     taux_xy(:,:) = taux_xy(:,:) + fluxbu(:,:) * dtfactor
     tauy_xy(:,:) = tauy_xy(:,:) + fluxbv(:,:) * dtfactor
     taux_top_xy(:,:) = taux_top_xy(:,:) + fluxtu(:,:) * dtfactor
     tauy_top_xy(:,:) = tauy_top_xy(:,:) + fluxtv(:,:) * dtfactor

   end if ! dodns

   call diffuse_mom()

end subroutine sgs_mom

!----------------------------------------------------------------------
!!! compute sgs diffusion of scalars:
!
subroutine sgs_scalars()

  use vars
  use microphysics, only: index_water_vapor, fluxbmk, fluxtmk, mkdiff, mkwsb, &
           flag_precip, flag_advect, micro_field, q, micro_flux, total_water
  use chemistry
  use tracers
  use params, only: dotracers, doterrain, docloud, doprecip, dosmoke, dochem
  use terrain, only: k_terra, terra, terraw, kmax
  use consts, only: cp, lcond
  use stat_coars

  implicit none

   real dummy(nz),dummy1(nzm)
   real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
   real, allocatable :: tkk(:,:,:)
   real q_bot, q_top, qsatw, rdz, rdz1, coef, rdx, rdy
   integer i,j,k,it,jt
   real buf(nz)
 
   call t_startf ('sgs_scalars')

!=====================================================================
   if(dodns) then

     if(dosmoke) then
       coef = 0.
     else
       coef = 1.
     end if

     q_bot = qsatw(t_bot,pres(1))*coef
     q_top = qsatw(t_top,pres(nzm))*coef

! fill temperature and water inside terrain for lateral fluxes:
     if(doterrain) then
      do k=1,kmax
       do j=0,ny+1
        do i=0,nx+1
         t(i,j,k) = terra(i,j,k)*t(i,j,k)+(1.-terra(i,j,k))*(t_bot + gamaz(k))
         q(i,j,k) = terra(i,j,k)*q(i,j,k)+(1.-terra(i,j,k))*q_bot
        end do
       end do
      end do
     end if

! compute fluxes at the bottom and top:

     rdz = 2./dz
     rdz1 = 2./(dz*adzw(nz))
     do j=1,ny
      do i=1,nx
       k = k_terra(i,j)
       fluxbt(i,j) = -rdz*DIFF_DNS*(t(i,j,k)-t_bot-gamaz(k))
       fluxtt(i,j) = -rdz1*DIFF_DNS*(t_top+gamaz(nzm)-t(i,j,nzm))
       shf_all(i,j) = cp*rhow(k)*fluxbt(i,j)
       fluxbq(i,j) = -rdz*DIFF_DNS*(q(i,j,k)-q_bot)
       fluxtq(i,j) = -rdz1*DIFF_DNS*(q_top-q(i,j,nzm))
       lhf_all(i,j) = lcond*rhow(k)*fluxbq(i,j)
      end do
     end do
     shf_xy(:,:) = shf_xy(:,:) + fluxbt(:,:) * dtfactor
     lhf_xy(:,:) = lhf_xy(:,:) + fluxbq(:,:) * dtfactor
     shf_top(:,:) = cp*rhow(nz)*fluxtt(:,:)
     lhf_top(:,:) = lcond*rhow(nz)*fluxtq(:,:)
     shf_top_xy(:,:) = shf_top_xy(:,:) + fluxtt(:,:) * dtfactor
     lhf_top_xy(:,:) = lhf_top_xy(:,:) + fluxtq(:,:) * dtfactor


     if(dochamber) then
       if(dostatis) then
        fluxwallt(:,:) = 0.
        fluxwallq(:,:) = 0.
       end if
       call task_rank_to_index(rank,it,jt)
       if(it.eq.0) then
        if(dosidewalls) t_wall = t_bot
        do k=1,nzm 
         t(:0,:,k) = t_wall+gamaz(k) 
         q(:0,:,k) = qsatw(t_wall,pres(k))*coef 
        end do
        if(dostatis) then
         do k=1,nzm
          do j=1,ny
            fluxwallt(k,1) = fluxwallt(k,1) - rho(k)*(t(1,j,k)-t(0,j,k))
            fluxwallq(k,1) = fluxwallq(k,1) - rho(k)*(q(1,j,k)-q(0,j,k))
          end do
         end do
        end if
       end if
       if(it+nx.eq.nx_gl) then
        if(dosidewalls) t_wall = t_top
        do k=1,nzm
         t(nx+1:,:,k) = t_wall+gamaz(k)
         q(nx+1:,:,k) = qsatw(t_wall,pres(k))*coef
        end do
        if(dostatis) then
         do k=1,nzm
          do j=1,ny
            fluxwallt(k,3) = fluxwallt(k,3) - rho(k)*(t(nx+1,j,k)-t(nx,j,k))
            fluxwallq(k,3) = fluxwallq(k,3) - rho(k)*(q(nx+1,j,k)-q(nx,j,k))
          end do
         end do
        end if
       end if
       if(jt.eq.0) then
        if(dosidewalls) t_wall = t_bot
        do k=1,nzm
         t(:,:0,k) = t_wall+gamaz(k)
         q(:,:0,k) = qsatw(t_wall,pres(k))*coef
        end do
        if(dostatis) then
         do k=1,nzm
          do i=1,nx
            fluxwallt(k,4) = fluxwallt(k,4) - rho(k)*(t(i,1,k)-t(i,0,k))
            fluxwallq(k,4) = fluxwallq(k,4) - rho(k)*(q(i,1,k)-q(i,0,k))
          end do
         end do
        end if
       end if
       if(jt+ny.eq.ny_gl) then
        if(dosidewalls) t_wall = t_top
        do k=1,nzm
         t(:,ny+1:,k) = t_wall+gamaz(k)
         q(:,ny+1:,k) = qsatw(t_wall,pres(k))*coef
        end do
        if(dostatis) then
         do k=1,nzm
          do i=1,nx
            fluxwallt(k,2) = fluxwallt(k,2) - rho(k)*(t(i,ny+1,k)-t(i,ny,k))
            fluxwallq(k,2) = fluxwallq(k,2) - rho(k)*(q(i,ny+1,k)-q(i,ny,k))
          end do
         end do
        end if
       end if
       if(dostatis) then
        if(dompi) then
         call task_sum_real(fluxwallt(:,1),buf,nzm)
         fluxwallt(:,1) = buf(:)
         call task_sum_real(fluxwallt(:,2),buf,nzm)
         fluxwallt(:,2) = buf(:)
         call task_sum_real(fluxwallt(:,3),buf,nzm)
         fluxwallt(:,3) = buf(:)
         call task_sum_real(fluxwallt(:,4),buf,nzm)
         fluxwallt(:,4) = buf(:)
         call task_sum_real(fluxwallq(:,1),buf,nzm)
         fluxwallq(:,1) = buf(:)
         call task_sum_real(fluxwallq(:,2),buf,nzm)
         fluxwallq(:,2) = buf(:)
         call task_sum_real(fluxwallq(:,3),buf,nzm)
         fluxwallq(:,3) = buf(:)
         call task_sum_real(fluxwallq(:,4),buf,nzm)
         fluxwallq(:,4) = buf(:)
        end if
        coef = DIFF_DNS*2./(dy*ny_gl)
        fluxwallt(:,1) = fluxwallt(:,1)*coef*cp 
        fluxwallt(:,3) = fluxwallt(:,3)*coef*cp 
        fluxwallq(:,1) = fluxwallq(:,1)*coef*lcond 
        fluxwallq(:,3) = fluxwallq(:,3)*coef*lcond 
        coef = DIFF_DNS*2./(dx*nx_gl)
        fluxwallt(:,2) = fluxwallt(:,2)*coef*cp 
        fluxwallt(:,4) = fluxwallt(:,4)*coef*cp 
        fluxwallq(:,2) = fluxwallq(:,2)*coef*lcond 
        fluxwallq(:,4) = fluxwallq(:,4)*coef*lcond 
       end if
       
     end if ! dochamber

   end if ! dodns
!=====================================================================

    call diffuse_scalar(t,fluxbt,fluxtt,tkh,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.,.true.)
    if(collect_coars) then
        call coars_fld(misc,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4))
        call coars_fld(misc1,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4+nmicro_fields+2))
    end if
    
    if(advect_sgs) then
         allocate (tkk(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm),stat=i) 
         if(i.gt.0) call task_abort_msg("diffuse scalars: alloc tkk failed!") 
         do k = 1,nzm
          do j = dimy1_d,dimy2_d
           do i = dimx1_d,dimx2_d 
            tkk(i,j,k) = 3.*tkh(i,j,k)
           end do
          end do
         end do
         call diffuse_scalar(tke,fzero,fzero,tkk,dummy,sgswsb, &
                                    dummy1,dummy1,dummy1,.false.,.false.)
         deallocate (tkk)
    end if


!
!    diffusion of microphysics prognostics:
!
    call micro_flux()

    total_water_evap = total_water_evap - total_water()

    do k = 1,nmicro_fields
        if(collect_coars) misc = 0.
        if(   k.eq.index_water_vapor ) then    ! transport water-vapor variable no metter what
           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k) 
           call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, &
                tkh,mkdiff(:,k),mkwsb(:,k), q2lediff,q2lediss,qwlediff,.true.,.false.)
        else if(docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k) 
           if(flag_advect(k).eq.1) call diffuse_scalar(micro_field(:,:,:,k), &
              fluxbtmp,fluxttmp, tkh,mkdiff(:,k),mkwsb(:,k), dummy1,dummy1,dummy1,.false.,.false.)
           if(collect_coars) call coars_fld(misc,mu(1:ny),ady(1:ny), &
                                         terraw(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4+k))
        end if
        if(collect_coars) call coars_fld(misc,mu(1:ny),ady(1:ny), &
                                         terraw(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4+k))
    end do

    total_water_evap = total_water_evap + total_water()

 ! diffusion of tracers:

    if(dotracers) then

      call tracers_flux()

      do k = 1,ntracers

        fluxbtmp = fluxbtr(:,:,k)
        fluxttmp = fluxttr(:,:,k)
        call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, tkh, &
             trdiff(:,k),trwsb(:,k),dummy1,dummy1,dummy1,.false.,.false.)

      end do
      do j=1,ny
       do i = 1,nx
         tr_ac_xy(i,j,1:ntracers) = tr_ac_xy(i,j,1:ntracers) - dtn*fluxbtr(i,j,1:ntracers)*rhow(k_terra(i,j))
       end do
      end do


    end if

    if ( dochem ) then
      call chem_flux()

      do k = 1,ngchem_fields
        
          fluxbtmp = fluxbch(:,:,k)
          fluxttmp = fluxtch(:,:,k)

          call diffuse_scalar(gchem_field(:,:,:,k),fluxbtmp,fluxttmp, &
                tkh, gchdiff(:,k),gchwsb(:,k), dummy, dummy, dummy, .false., .false.) 
      end do
    end if

    call t_stopf ('sgs_scalars')


end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc()

   use grid, only: dt,dtn,dx,dy,dz,adz,ady,adyv,adzw,mu,muv,zi
   use params, only : nub, dodamping_w, doimplicitdiff
   use consts, only: pi
   use vars, only: taudamp
   real cx, cy, cz
   integer k, j, kb
   real tk_factorz, nu, coef

!    SGS TKE equation:

if(dosgs) then

  if(.not.dodns) call tke_full()

  if(doimplicitdiff) then
    coef = 0.
  else
    coef = 1.
  end if

  do k=1,nzm
            ! make high horizontal eddy diffusivity in the region of spange layer
            ! to fight the noise there MK 01/2023
            nu = (zi(k)-zi(1))/(zi(nzm)-zi(1))
            if(dodamping_w.and.nu.gt.nub) then
              tk_factorz = max(tk_factor,0.01*taudamp(k)*dt/dtn)
            else
              tk_factorz = tk_factor
            end if
            cz = coef*dtn/((dz*min(adzw(k),adzw(k+1)))**2)
            do j=1-YES3D,ny
             cx = dtn/(dx*mu(j))**2
             cy = dtn/(dy*min(adyv(j),adyv(j+YES3D)))**2
             tkmax(j,k) = cfl_diffsc_max/(cx+cy+cz)  ! maximum value of eddy visc/cond
             tkmin(j,k) = tk_factorz*tkmax(j,k)  ! minimum value of eddy visc/cond
             tkmaxu(j,k) = tkmax(j,k)*cfl_diff_max/cfl_diffsc_max
             tkminu(j,k) = tk_factorz*tkmaxu(j,k)
             cx = dtn/(dx*muv(j))**2
             cy = dtn/(dy*min(ady(j),ady(j+YES3D)))**2
             tkmaxv(j,k) = cfl_diff_max/(cx+cy+cz)  ! maximum value of eddy visc/cond
             tkminv(j,k) = tk_factorz*tkmaxv(j,k)  ! minimum value of eddy visc/cond
            end do
            kb=max(1,k-1)
            cz = coef*dtn/(dz*min(adz(k),adz(kb)))**2
            do j=1-YES3D,ny
             cx = dtn/(dx*mu(j))**2
             cy = dtn/(dy*min(adyv(j),adyv(j+YES3D)))**2
             tkmaxw(j,k) = cfl_diff_max/(cx+cy+cz)  ! maximum value of eddy visc/cond
             tkminw(j,k) = tk_factorz*tkmaxw(j,k)  ! minimum value of eddy visc/cond
            end do
  end do

end if

end subroutine sgs_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine sgs_diagnose()
! None 

end subroutine sgs_diagnose


!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!
subroutine sgs_statistics()
  
  use vars
  use hbuffer, only: hbuf_put, hbuf_avg_put
  use params, only : lcond

  real factor_xy 
  real tkz(nzm), tkhz(nzm)
  integer i,j,k,n
  character(LEN=6) :: statname  !bloss: for conditional averages

  if(.not.dosgs) return

  call t_startf ('sgs_statistics')

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
    tkz(k) = 0.
    tkhz(k) = 0.
    do j=1,ny
    do i=1,nx
      tkz(k)=tkz(k)+tk(i,j,k)
      tkhz(k)=tkhz(k)+tkh(i,j,k)
    end do
    end do
  end do

  call hbuf_avg_put('TKES',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)

  call hbuf_put('TK',tkz,factor_xy)
  call hbuf_put('TKH',tkhz,factor_xy)

!---------------------------------------------------------
! SGS TKE Budget:

         call hbuf_put('ADVTRS',sgswle(:,1),factor_xy)
         call hbuf_put('BUOYAS',tkesbbuoy,factor_xy)
         call hbuf_put('SHEARS',tkesbshear,factor_xy)
         call hbuf_put('DISSIPS',tkesbdiss,factor_xy)

  if(dochamber) then
         call hbuf_put('FLXWLWT',fluxwallt(:,1),1.)
         call hbuf_put('FLXWLNT',fluxwallt(:,2),1.)
         call hbuf_put('FLXWLET',fluxwallt(:,3),1.)
         call hbuf_put('FLXWLST',fluxwallt(:,4),1.)
         call hbuf_put('FLXWLWQ',fluxwallq(:,1),1.)
         call hbuf_put('FLXWLNQ',fluxwallq(:,2),1.)
         call hbuf_put('FLXWLEQ',fluxwallq(:,3),1.)
         call hbuf_put('FLXWLSQ',fluxwallq(:,4),1.)
  end if

  call t_stopf ('sgs_statistics')

end subroutine sgs_statistics

!----------------------------------------------------------------------
! called when stepout() called

subroutine sgs_print()

 call fminmax_print('tke:',tke,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
 call fminmax_print('tk:',tk,0,nxp1,1-YES3D,nyp1,nzm)
 call fminmax_print('tkh:',tkh,0,nxp1,1-YES3D,nyp1,nzm)
 

end subroutine sgs_print

!----------------------------------------------------------------------
!!! Initialize the list of sgs statistics 
!
subroutine sgs_hbuf_init(namelist,deflist,unitlist,status,average_type,count,sgscount)
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,sgscount

end subroutine sgs_hbuf_init


end module sgs



