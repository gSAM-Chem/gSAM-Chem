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

real grdf_x(0:ny,nzm)! grid factor for eddy diffusion in x
real grdf_y(0:ny,nzm)! grid factor for eddy diffusion in y
real grdfu_x(0:ny,nzm)! grid factor for eddy diffusion in x
real grdfu_y(0:ny,nzm)! grid factor for eddy diffusion in y
real grdfv_x(0:ny,nzm)! grid factor for eddy diffusion in x
real grdfv_y(0:ny,nzm)! grid factor for eddy diffusion in y
real grdfw_x(0:ny,nzm)! grid factor for eddy diffusion in x
real grdfw_y(0:ny,nzm)! grid factor for eddy diffusion in y

logical:: dosmagor   ! if true, then use Smagorinsky closure

! Local diagnostics:

real tkesbbuoy(nz), tkesbshear(nz),tkesbdiss(nz), tkesbdiff(nz)


CONTAINS

! required microphysics subroutines and function:
!----------------------------------------------------------------------
!!! Read microphysics options from prm (namelist) file

subroutine sgs_setparm()

  use grid, only: case, caseid, masterproc
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder

  !======================================================================
  NAMELIST /SGS_TKE/ &
       dosmagor ! Diagnostic Smagorinsky closure

  NAMELIST /BNCUIODSBJCB/ place_holder

  dosmagor = .true. ! default 

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./CASES/'//trim(case)//'/prm', status='old',form='formatted')

  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  read (55,SGS_TKE,IOSTAT=ios)

  advect_sgs = .not.dosmagor

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

  ! END UW ADDITION
  !======================================================================

end subroutine sgs_setparm

!----------------------------------------------------------------------
!!! Initialize sgs:


subroutine sgs_init()

   use grid, only: nrestart,dx,dy,dz,adz,ady,adyv,adzw,mu,muv,masterproc
   integer k, j

  if(nrestart.eq.0) then

     sgs_field = 0.
     sgs_field_diag = 0.

     fluxbsgs = 0.
     fluxtsgs = 0.

  end if

  if(masterproc) then
     if(dosmagor) then
        write(*,*) 'Smagorinsky SGS Closure'
     else
        write(*,*) 'Prognostic TKE 1.5-order SGS Closure'
     end if
  end if

  do k=1,nzm
       do j=0,ny
        grdf_x(j,k) = min(16.,(dx*mu(j)/(adz(k)*dz))**2)
        grdf_y(j,k) = min(16.,(dy*adyv(j+1)/(adz(k)*dz))**2)
        grdfu_x(j,k) = min(16.,(dx*mu(j)/(adz(k)*dz))**2)
        grdfu_y(j,k) = min(16.,(dy*adyv(j+1)/(adz(k)*dz))**2)
        grdfv_x(j,k) = min(16.,(dx*muv(j)/(adz(k)*dz))**2)
        grdfv_y(j,k) = min(16.,(dy*ady(j)/(adz(k)*dz))**2)
        grdfw_x(j,k) = min(16.,(dx*mu(j)/(adzw(k)*dz))**2)
        grdfw_y(j,k) = min(16.,(dy*adyv(j+1)/(adzw(k)*dz))**2)
       end do
  end do

  sgswle = 0.
  sgswsb = 0.
  sgsadv = 0.
  sgsdiff = 0.
  sgslsadv = 0.

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

use grid, only: dtn, dx, dy, dz, adz, adzw, imu, ady 
use params, only: doimplicitdiff
implicit none

real, intent(out) :: cfl

integer j,k, i
real tkhxmax(nzm),tkhymax(nzm),tkhmax(nzm)

tkhxmax(:) = 0.
tkhymax(:) = 0.
do k = 1,nzm
 tkhmax(k) = maxval(tkh(1:nx,1:ny,k))
 do j = 1,ny
  tkhxmax(k) = max(tkhxmax(k),maxval(tkh(:,j,k)*grdf_x(j,k)*imu(j)**2))
  tkhymax(k) = max(tkhymax(k),maxval(tkh(:,j,k)*grdf_y(j,k)/ady(j)**2))
 end do
end do

cfl = 0.
if(doimplicitdiff) then
  do k=1,nzm
   cfl = max(cfl,        &
     tkhxmax(k)*dtn/dx**2 + &
     YES3D*tkhymax(k)*dtn/dy**2)
  end do
else
  do k=1,nzm
    cfl = max(cfl,        &
     tkhmax(k)*dtn/(dz*adzw(k))**2 + &
     tkhxmax(k)*dtn/dx**2 + &
     YES3D*tkhymax(k)*dtn/dy**2)
  end do
end if

end subroutine kurant_sgs


!----------------------------------------------------------------------
!!! compute sgs diffusion of momentum:
!
subroutine sgs_mom()

   call diffuse_mom()

end subroutine sgs_mom

!----------------------------------------------------------------------
!!! compute sgs diffusion of scalars:
!
subroutine sgs_scalars()

  use vars
  use microphysics
  use tracers
  use params, only: dotracers
  use terrain, only: k_terra, terraw
  use stat_coars

  implicit none

   real dummy(nz)
   real fluxbtmp(nx,ny), fluxttmp(nx,ny) !bloss
   real, allocatable :: tkk(:,:,:)
   integer i,j,k
 
    call t_startf ('diffuse_scalars')

    call diffuse_scalar(t,fluxbt,fluxtt,tkh,tdiff,twsb, &
                           t2lediff,t2lediss,twlediff,.true.)
    if(collect_coars) then
        call coars_fld(misc,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4))
        call coars_fld(misc1,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4+nmicro_fields+2))
    end if
    
    if(advect_sgs) then
         allocate (tkk(dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm))
         do k = 1,nzm
          do j = dimy1_d,dimy2_d
           do i = dimx1_d,dimx2_d 
            tkk(i,j,k) = min(tkmax(j,k),3.*tkh(i,j,k))
           end do
          end do
         end do
         call diffuse_scalar(tke,fzero,fzero,tkk,dummy,sgswsb, &
                                    dummy,dummy,dummy,.false.)
         deallocate (tkk)
    end if


!
!    diffusion of microphysics prognostics:
!
    call micro_flux()

    total_water_evap = total_water_evap - total_water()

    do k = 1,nmicro_fields
        if(collect_coars) misc = 0.
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
           fluxbtmp(1:nx,1:ny) = fluxbmk(1:nx,1:ny,k)
           fluxttmp(1:nx,1:ny) = fluxtmk(1:nx,1:ny,k) 
              if(flag_advect(k).eq.1) call diffuse_scalar(micro_field(:,:,:,k),fluxbtmp,fluxttmp, tkh, &
                mkdiff(:,k),mkwsb(:,k), dummy,dummy,dummy,.false.)
         end if
         if(collect_coars) then
             call coars_fld(misc,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_flux_sgs(:,:,:,4+k))
         end if
    end do

    total_water_evap = total_water_evap + total_water()

 ! diffusion of tracers:

    if(dotracers) then

      call tracers_flux()

      do k = 1,ntracers

        fluxbtmp = fluxbtr(:,:,k)
        fluxttmp = fluxttr(:,:,k)
        call diffuse_scalar(tracer(:,:,:,k),fluxbtmp,fluxttmp, tkh, &
             trdiff(:,k),trwsb(:,k), &
             dummy,dummy,dummy,.false.)
!!$          call diffuse_scalar(tracer(:,:,:,k),fluxbtr(:,:,k),fluxttr(:,:,k), tkh, &
!!$                              trdiff(:,k),trwsb(:,k), &
!!$                           dummy,dummy,dummy,.false.)

      end do
      do j=1,ny
       do i = 1,nx
         tr_ac_xy(i,j,1:ntracers) = tr_ac_xy(i,j,1:ntracers) - dtn*fluxbtr(i,j,1:ntracers)*rhow(k_terra(i,j))
       end do
      end do


    end if

    call t_stopf ('diffuse_scalars')


end subroutine sgs_scalars

!----------------------------------------------------------------------
!!! compute sgs processes (beyond advection):
!
subroutine sgs_proc()

   use grid, only: dtn,dx,dy,dz,adz,ady,adyv,adzw,mu,muv
   use params, only : doimplicitdiff, tk_factor, cfl_diff_max
   real cx, cy, cz, cz1, coef
   integer k, j, kb


!    SGS TKE equation:

if(dosgs) then

  call tke_full()

  if(doimplicitdiff) then
    coef=1.e10
  else
    coef=1.
  end if
  do k=1,nzm
            cz1 = (dz*min(adzw(k),adzw(k+1)))**2/dtn
            cz = coef*cz1
            do j=1-YES3D,ny
             cx = (dx*mu(j))**2/dtn/grdf_x(j,k)
             cy = (dy*min(adyv(j),adyv(j+YES3D)))**2/dtn/grdf_y(j,k)
             tkmax(j,k) = cfl_diff_max/(1./cx+1./cy+1./cz)  ! maximum value of eddy visc/cond
             tkmin(j,k) = tk_factor*tkmax(j,k)  ! minimum value of eddy visc/cond
             tkmaxu(j,k) = tkmax(j,k)
             tkminu(j,k) = tkmin(j,k)
             cx = (dx*muv(j))**2/dtn/grdfv_x(j,k)
             cy = (dy*min(ady(j),ady(j-YES3D)))**2/dtn/grdfv_y(j,k)
             tkmaxv(j,k) = cfl_diff_max/(1./cx+1./cy+1./cz)  ! maximum value of eddy visc/cond
             tkminv(j,k) = tk_factor*tkmaxv(j,k)  ! minimum value of eddy visc/cond
            end do
  end do
  do k=1,nzm
            kb=max(1,k-1)
            cz1 = (dz*min(adz(k),adz(kb)))**2/dtn
            cz = coef*cz1
            do j=1-YES3D,ny
             cx = (dx*mu(j))**2/dtn/grdfw_x(j,k)
             cy = (dy*min(adyv(j),adyv(j+YES3D)))**2/dtn/grdfw_y(j,k)
             tkmaxw(j,k) = cfl_diff_max/(1./cx+1./cy+1./cz)  ! maximum value of eddy visc/cond
             tkminw(j,k) = tk_factor*tkmaxw(j,k)  ! minimum value of eddy visc/cond
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



