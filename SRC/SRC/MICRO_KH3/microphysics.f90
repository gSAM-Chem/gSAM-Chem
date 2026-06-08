module microphysics

! module for Khairoutdinov & Heymsfield (2015) bulk microphysics
! Marat Khairoutdinov, 2015

use grid, only: nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s ! subdomain grid information 
use params, only: doprecip, docloud
use micro_params
implicit none

!----------------------------------------------------------------------
!!! required definitions:

integer, parameter :: nmicro_fields = 3   ! total number of prognostic water vars

!!! microphysics prognostic variables are storred in this array:

real micro_field(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, nmicro_fields)

integer, parameter :: flag_wmass(nmicro_fields) = (/1,1,1/)
integer, parameter :: index_water_vapor = 1 ! index for variable that has water vapor
integer, parameter :: index_cloud_ice = -1   ! index for cloud ice (sedimentation)
integer, parameter :: flag_precip(nmicro_fields) = (/0,0,1/)

! all variables correspond to mass, not number
integer, parameter :: flag_number(nmicro_fields) = (/0,0,0/)
integer, parameter :: flag_advect(nmicro_fields) = (/1,1,1/) ! dummy parameter (don't change -MK)

!  3D and 2DZ microphysical fields are output by default.
integer :: flag_micro3Dout(nmicro_fields) = (/0,0,0/)
integer, parameter :: flag_micro2DZout(nmicro_fields) = (/0,0,0/)
integer total_micro_3Dout

real fluxbmk (nx, ny, 1:nmicro_fields) ! surface flux of tracers
real fluxtmk (nx, ny, 1:nmicro_fields) ! top boundary flux of tracers

!!! these arrays are needed for output statistics:

real mkwle(nz,1:nmicro_fields)  ! resolved vertical flux
real mkwsb(nz,1:nmicro_fields)  ! SGS vertical flux
real mkadv(nz,1:nmicro_fields)  ! tendency due to vertical advection
real mklsadv(nz,1:nmicro_fields)  ! tendency due to large-scale vertical advection
real mkdiff(nz,1:nmicro_fields)  ! tendency due to vertical diffusion

!======================================================================

!bloss: arrays with names/units for microphysical outputs in statistics.
character*3, dimension(nmicro_fields) :: mkname
character*80, dimension(nmicro_fields) :: mklongname
character*10, dimension(nmicro_fields) :: mkunits
real, dimension(nmicro_fields) :: mkoutputscale

!bloss: dummy arrays for effective radius and other variables
!    useful in computing cloud radiative properties, mainly
!    used by M2005 microphysics but included here so that
!    RRTM will compile against both schemes
real, allocatable, dimension(:,:,:) :: reffc, reffi, reffs, &
     CloudLiquidMassMixingRatio, CloudLiquidGammaExponent, CloudLiquidLambda, &
     CloudIceMassMixingRatio, SnowMassMixingRatio

! Flags related to M2005 cloud optics routines, included here
!   so that RRTM will compile against this microphysics as well.
logical :: dosnow_radiatively_active = .false.
logical :: dorrtm_cloud_optics_from_effrad_LegacyOption = .true.
logical :: doallice_radiatively_active = .false. ! treat all ice as radiative active

!======================================================================

!------------------------------------------------------------------
! Optional (internal) definitions)

! make aliases for prognostic variables:
! note that the aliases should be local to microphysics

real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)   ! water vapor + liquid cloud water
real qi(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! cloud ice/snow
real qp(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)  ! graupel + rain
equivalence (q(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,1))
equivalence (qi(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,2))
equivalence (qp(dimx1_s,dimy1_s,1),micro_field(dimx1_s,dimy1_s,1,3))

real qc(nx,ny,nzm)  ! cloud liquid water
real qn(nx,ny,nzm)  ! cloud liquid water
equivalence (qc(1,1,1),qn(1,1,1))

real qpsrc(nz)  ! source of precipitation (including ice) microphysical processes
real qpfall(nz) ! source of precipitating water due to fall out in a given level
real qpevp(nz)  ! sink of precipitating water due to evaporation

real vice(nzm),vrain,vgrau,cice(nzm),crain, cgrau  ! precomputed coefs for precip terminal velocity

CONTAINS

! required microphysics subroutines and function:

!copied from Bloss
!----------------------------------------------------------------------
function micro_scheme_name()
  character(len=32) :: micro_scheme_name
  ! Return the scheme name, normally the same as the directory name with leading
  ! "MICRO_" removed
  micro_scheme_name = "kh3"
end function   micro_scheme_name

!----------------------------------------------------------------------
!!! Read microphysics options from prm file

subroutine micro_setparm()

  use grid, only: case, caseid, masterproc
  implicit none

  integer ierr

  ! no user-definable options in this  microphysics.
   NAMELIST /MICRO_KH3/ qcw0, qci0, icefall_fudge
  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  if(masterproc) print*,'Microphysics: KH3'
  open(55,file='./CASES/'//trim(case)//'/prm', status='old',form='formatted')
  read (55,MICRO_KH3,IOSTAT=ierr,end=111)
  if (ierr.ne.0) then
     !namelist error checking
        write(*,*) '****** ERROR: bad specification in MICRO_SAM1MOM namelist'
        rewind(55)
        read (55,MICRO_KH3)
        call task_abort()
  end if
111  close(55)
! write namelist values out to file for documentation
  if(masterproc) then
      open(unit=55,file='./OUT_STAT/'//trim(case)//'_'//trim(caseid)//'.nml',&
            form='formatted',position='append')
      write (55,nml=MICRO_KH3)
      write(55,*)
      close(55)
  end if


end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:


subroutine micro_init()

  use grid, only: nrestart
  use vars, only: qv, q0, qci, qcl, qpi, qpl
  use params, only: dosmoke, readinit
  integer k

  a_pr = 1./(tprmax-tprmin)

  call precip_init() 

  if(nrestart.eq.0) then

     micro_field = 0.
     if(readinit) then
        q(1:nx,1:ny,1:nzm) = qv(1:nx,1:ny,1:nzm)+qcl(1:nx,1:ny,1:nzm)
        qi(1:nx,1:ny,1:nzm) = qci(1:nx,1:ny,1:nzm)
        qp(1:nx,1:ny,1:nzm) = qpl(1:nx,1:ny,1:nzm)+qpi(1:nx,1:ny,1:nzm)
     else
       do k=1,nzm
        q(:,:,k) = q0(k)
       end do
     end if
     qc = 0.
     fluxbmk = 0.
     fluxtmk = 0.

     if(docloud) then
       call cloud()
       call micro_diagnose()
     end if
     if(dosmoke) then
       call micro_diagnose()
     end if

  end if

  mkwle = 0.
  mkwsb = 0.
  mkadv = 0.
  mkdiff = 0.
  mklsadv = 0.

  qpsrc = 0.
  qpevp = 0.

  mkname(1) = 'QT'
  mklongname(1) = 'TOTAL WATER (VAPOR + LIQUID CONDENSATE)'
  mkunits(1) = 'g/kg'
  mkoutputscale(1) = 1.e3

  mkname(2) = 'QI'
  mklongname(2) = 'PRISTINE ICE + SNOW'
  mkunits(2) = 'g/kg'
  mkoutputscale(2) = 1.e3

  mkname(3) = 'QP'
  mklongname(3) = 'RAIN+GRAUPEL WATER'
  mkunits(3) = 'g/kg'
  mkoutputscale(3) = 1.e3

  total_micro_3Dout = 0
  if(docloud) then
      flag_micro3Dout(1) = 1
      total_micro_3Dout = 2
  end if
  if(doprecip) then
      flag_micro3Dout(3) = 0
      total_micro_3Dout = total_micro_3Dout + 2
  end if
  total_micro_3Dout = total_micro_3Dout + sum(flag_micro3Dout)

end subroutine micro_init

!----------------------------------------------------------------------
! given temperature and bulk cloud/precip fields, make them internally
! consistent
! with current microphysics scheme and keeping temperature the same

subroutine micro_set(tabs, t, qv, qcl, qci, qpl, qpi, qpg)
  use vars, only: gamaz, pres
  use params, only: fac_cond, fac_sub
  real, intent(in) ::  tabs(nx,ny,nzm)
  real, intent(inout) ::  qv(nx,ny,nzm), qcl(nx,ny,nzm), qci(nx,ny,nzm)
  real, intent(inout) ::  qpl(nx,ny,nzm), qpi(nx,ny,nzm), qpg(nx,ny,nzm), t(nx,ny,nzm)
  real om, qq, dq, qsatt
  integer i,j,k
  real, external :: qsatw,qsati
      a_pr = 1./(tprmax-tprmin)
      do k=1,nzm
       do j=1,ny
        do i=1,nx
          qv(i,j,k) = max(0.,qv(i,j,k))
          if(pres(k).lt.50.) then  ! no cloud above 50 mb level
           qcl(i,j,k) = 0.
           qci(i,j,k) = 0.
          else
           qcl(i,j,k) = max(0.,qcl(i,j,k))
           qci(i,j,k) = max(0.,qci(i,j,k))
          end if
          om = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
          qq = qpl(i,j,k) + qpi(i,j,k)
          qpl(i,j,k) = max(0.,qq)*om
          qpi(i,j,k) = max(0.,qq)*(1.-om)
          qpg(i,j,k) = qpi(i,j,k)
          t(i,j,k)=tabs(i,j,k)+gamaz(k)-fac_cond*(qcl(i,j,k)+qpl(i,j,k)) &
                                        -fac_sub*(qci(i,j,k)+qpi(i,j,k))
        end do
       end do
      end do
end subroutine micro_set


!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
subroutine micro_flux()

  use vars, only: fluxbq, fluxtq

  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:)
  fluxtmk(:,:,index_water_vapor) = fluxtq(:,:)

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (bayond advection and SGS diffusion):
!
subroutine micro_proc()

   use grid, only: nstep,dt,icycle
   use params, only: dosmoke

   call t_startf ('micro_proc')


   ! Update bulk coefficient
   if(icycle.eq.1) call precip_init() 

   if(docloud) then
     call cloud()
     if(doprecip) call precip_proc()
     call micro_diagnose()
   end if
   if(dosmoke) then
     call micro_diagnose()
   end if

   call t_stopf ('micro_proc')

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and statistics:
!
subroutine micro_diagnose()
 
   use vars

   real omn, omp
   integer i,j,k

   do k=1,nzm
    do j=1,ny
     do i=1,nx
       qv(i,j,k) = q(i,j,k) - qc(i,j,k)
       qcl(i,j,k) = qc(i,j,k)
       qci(i,j,k) = qi(i,j,k)
       omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
       qpl(i,j,k) = qp(i,j,k)*omp
       qpi(i,j,k) = qp(i,j,k)*(1.-omp)
     end do
    end do
   end do
   where(qcl(1:nx,1:ny,1:nzm)+qci(1:nx,1:ny,1:nzm).gt.0)
      cld(1:nx,1:ny,1:nzm) = 1.
   elsewhere
      cld(1:nx,1:ny,1:nzm) = 0.
   end where

end subroutine micro_diagnose

!----------------------------------------------------------------------
!!! function to compute terminal velocity for precipitating variables:
! In this particular case there is only one precipitating variable.

real function term_vel_qp(i,j,k,ind)
  
  use vars
  integer, intent(in) :: i,j,k,ind
  real omp, qrr, qgg

  if(ind.eq.1) then
    term_vel_qp = vice(k)*qi(i,j,k)**cice(k)
  else
    term_vel_qp = 0.
    if(qp(i,j,k).gt.qp_threshold) then
      omp = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
      if(omp.eq.1.) then
         term_vel_qp = vrain*(rho(k)*qp(i,j,k))**crain
      elseif(omp.eq.0.) then
         term_vel_qp = vgrau*(rho(k)*qp(i,j,k))**cgrau 
      else
         qrr=omp*qp(i,j,k)
         qgg=qp(i,j,k)-qrr
         term_vel_qp = omp*vrain*(rho(k)*qrr)**crain+(1.-omp)*vgrau*(rho(k)*qgg)**cgrau 
      endif
    end if  
  end if
end function term_vel_qp

!----------------------------------------------------------------------
!!! compute sedimentation 
!
subroutine micro_precip_fall()
  
  use vars
  use params, only : pi

  real omega(nx,ny,nzm)
  real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real f0(nzm),df0(nzm)
  real dummy(1)
  integer ind
  integer i,j,k

  crain = b_rain / 4.
  cgrau = b_grau / 4.
  vrain = a_rain * gamma(4.+b_rain) / 6. / (pi * rhor * nzeror) ** crain
  vgrau = a_grau * gamma(4.+b_grau) / 6. / (pi * rhog * nzerog) ** cgrau
  do k=1,nzm
   cice(k) = bv_ice(k)/b_ice
   vice(k) = av_ice(k)*gamma(b_ice+bv_ice(k)+gami(k)+1)/gamma(b_ice+gami(k)+1) &
           *(gamma(gami(k)+1)*rho(k)/(a_ice*niz(k)*gamma(b_ice+gami(k)+1)))**(bv_ice(k)/b_ice)
  end do

! Initialize arrays that accumulate surface precipitation flux

 if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
    do i=1,nx
     precsfc(i,j)=0.
    end do
   end do
   do k=1,nzm
    precflux(k) = 0.
   end do
 end if

 do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
    qpfall(k)=0.
    tlat(k) = 0.
    qifall(k)=0.
 end do
   
 do k=1,nzm
  do j=1,ny
   do i=1,nx
       omega(i,j,k) = 0.
   end do
  end do
 end do
 call precip_fall(qi, term_vel_qp, 1, omega, 1)

 do k=1,nzm
  do j=1,ny
   do i=1,nx
       omega(i,j,k) = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
   end do
  end do
 end do
 call precip_fall(qp, term_vel_qp, 2, omega, 0)

 if(dostatis) then
   do k=1,nzm
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,dimx2_s
          df(i,j,k) = t(i,j,k)
       end do
     end do
   end do
   call stat_varscalar(t,df,f0,df0,t2leprec)
   call setvalue(twleprec,nzm,0.)
   call stat_sw2(t,df,twleprec)
 endif

end subroutine micro_precip_fall

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!
subroutine micro_statistics()
  
  use vars
  use hbuffer, only: hbuf_put
  use params, only : lcond

  real tmp(2), factor_xy 
  real qcz(nzm), qiz(nzm), qrz(nzm), qsz(nzm), qgz(nzm)
  integer i,j,k,n
  character(LEN=6) :: statname  !bloss: for conditional averages

  call t_startf ('micro_statistics')

  factor_xy = 1./float(nx*ny)

  do k=1,nzm
      tmp(1) = dz/rhow(k)
      tmp(2) = tmp(1) / dtn
      mkwsb(k,1) = mkwsb(k,1) * tmp(1) * rhow(k) * lcond
      mkwle(k,1) = mkwle(k,1)*tmp(2)*rhow(k)*lcond + mkwsb(k,1)
      if(docloud.and.doprecip) then
        mkwsb(k,2) = mkwsb(k,2) * tmp(1) * rhow(k) * lcond
        mkwle(k,2) = mkwle(k,2)*tmp(2)*rhow(k)*lcond + mkwsb(k,2)
      endif
  end do

  call hbuf_put('QTFLUX',mkwle(:,1),factor_xy)
  call hbuf_put('QTFLUXS',mkwsb(:,1),factor_xy)
  call hbuf_put('QPFLUX',mkwle(:,2),factor_xy)
  call hbuf_put('QPFLUXS',mkwsb(:,2),factor_xy)

  do k=1,nzm
    qcz(k) = 0.
    qiz(k) = 0.
    qrz(k) = 0.
    qsz(k) = 0.
    qgz(k) = 0.
    do j=1,ny
    do i=1,nx
      qcz(k)=qcz(k)+qcl(i,j,k)
      qiz(k)=qiz(k)+qci(i,j,k)
      qrz(k)=qrz(k)+qpl(i,j,k)
      qgz(k)=qgz(k)+qpi(i,j,k)
    end do
    end do
  end do

  call hbuf_put('QC',qcz,1.e3*factor_xy)
  call hbuf_put('QI',qiz,1.e3*factor_xy)
  call hbuf_put('QR',qrz,1.e3*factor_xy)
  call hbuf_put('QG',qgz,1.e3*factor_xy)
  call hbuf_put('QS',qsz,1.e3*factor_xy)

  call hbuf_put('QTADV',mkadv(:,1)+qifall,factor_xy*86400000./dtn)
  call hbuf_put('QTDIFF',mkdiff(:,1),factor_xy*86400000./dtn)
  call hbuf_put('QTSINK',qpsrc,-factor_xy*86400000./dtn)
  call hbuf_put('QTSRC',qpevp,-factor_xy*86400000./dtn)
  call hbuf_put('QPADV',mkadv(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPDIFF',mkdiff(:,2),factor_xy*86400000./dtn)
  call hbuf_put('QPFALL',qpfall,factor_xy*86400000./dtn)
  call hbuf_put('QPSRC',qpsrc,factor_xy*86400000./dtn)
  call hbuf_put('QPEVP',qpevp,factor_xy*86400000./dtn)

  do n = 1,nmicro_fields
     call hbuf_put(trim(mkname(n))//'LSADV', &
          mklsadv(:,n),mkoutputscale(n)*factor_xy*86400.)
  end do

  do n = 1,ncondavg

     do k=1,nzm
        qcz(k) = 0.
        qiz(k) = 0.
        qrz(k) = 0.
        qsz(k) = 0.
        qgz(k) = 0.
        do j=1,ny
           do i=1,nx
              qcz(k)=qcz(k)+qcl(i,j,k)*condavg_mask(i,j,k,n)
              qiz(k)=qiz(k)+qci(i,j,k)*condavg_mask(i,j,k,n)
              qrz(k)=qrz(k)+qpl(i,j,k)*condavg_mask(i,j,k,n)
              qgz(k)=qgz(k)+qpi(i,j,k)*condavg_mask(i,j,k,n)
           end do
        end do
     end do

     call hbuf_put('QC' // TRIM(condavgname(n)),qcz,1.e3)
     call hbuf_put('QI' // TRIM(condavgname(n)),qiz,1.e3)
     if(doprecip) then
        call hbuf_put('QR' // TRIM(condavgname(n)),qrz,1.e3)
        call hbuf_put('QG' // TRIM(condavgname(n)),qgz,1.e3)
        call hbuf_put('QS' // TRIM(condavgname(n)),qrz,1.e3)
     end if
  end do

  ncmn = 0.
  nrmn = 0.

  call t_stopf ('micro_statistics')

end subroutine micro_statistics

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
end subroutine micro_print

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics 
!
subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,trcount)

  use vars


   character(*) namelist(*), deflist(*), unitlist(*)
   integer status(*),average_type(*),count,trcount
   integer ntr, n


   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUX'
   deflist(count) = 'Nonprecipitating water flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QTFLUXS'
   deflist(count) = 'Nonprecipitating-water flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUX'
   deflist(count) = 'Precipitating-water turbulent flux (Total)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QPFLUXS'
   deflist(count) = 'Precipitating-water turbulent flux (SGS)'
   unitlist(count) = 'W/m2'
   status(count) = 1    
   average_type(count) = 0
   
   do n = 1,nmicro_fields
      count = count + 1
      trcount = trcount + 1
      namelist(count) = TRIM(mkname(n))//'LSADV'
      deflist(count) = 'Source of '//TRIM(mklongname(n))//' due to large-scale vertical advection'
      unitlist(count) = TRIM(mkunits(n))//'day'
      status(count) = 1    
      average_type(count) = 0
   end do

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QC'
   deflist(count) = 'Liquid cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QI'
   deflist(count) = 'Icy cloud water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QR'
   deflist(count) = 'Rain water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QS'
   deflist(count) = 'Snow water'
   unitlist(count) = 'g/kg'
   status(count) = 1
   average_type(count) = 0

   count = count + 1
   trcount = trcount + 1
   namelist(count) = 'QG'
   deflist(count) = 'Graupel water'
   unitlist(count) = 'g/kg'
   status(count) = 1    
   average_type(count) = 0

  !bloss: setup to add an arbitrary number of conditional statistics
   do n = 1,ncondavg

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QC' // TRIM(condavgname(n))
      deflist(count) = 'Mean Liquid cloud water in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      count = count + 1
      trcount = trcount + 1
      namelist(count) = 'QI' // TRIM(condavgname(n))
      deflist(count) = 'Mean Icy cloud water in ' // TRIM(condavglongname(n))
      unitlist(count) = 'g/kg'
      status(count) = 1    
      average_type(count) = n

      if(doprecip) then
         count = count + 1
         trcount = trcount + 1
         namelist(count) = 'QR' // TRIM(condavgname(n))
         deflist(count) = 'Mean Rain water in ' // TRIM(condavglongname(n))
         unitlist(count) = 'g/kg'
         status(count) = 1    
         average_type(count) = n

         count = count + 1
         trcount = trcount + 1
         namelist(count) = 'QS' // TRIM(condavgname(n))
         deflist(count) = 'Mean Snow water in ' // TRIM(condavglongname(n))
         unitlist(count) = 'g/kg'
         status(count) = 1
         average_type(count) = n

         count = count + 1
         trcount = trcount + 1
         namelist(count) = 'QG' // TRIM(condavgname(n))
         deflist(count) = 'Mean Graupel water in ' // TRIM(condavglongname(n))
         unitlist(count) = 'g/kg'
         status(count) = 1    
         average_type(count) = n
      end if

   end do

end subroutine micro_hbuf_init


subroutine micro_write_fields3D(nfields1)

  use grid, only: nsubdomains,save3Dbin,nfiles3D
  use vars, only: qv, qcl, qci, qpl, qpi, qpg
  use mpi_stuff, only: dompi, rank
  integer, intent(inout) :: nfields1 ! used by writefields3D as a field counter
  real(4) tmp(nx,ny,nzm)
  character *80 long_name
  character *8 name
  character *10 units

  integer i,j,k

  if(docloud) then
    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=qcl(i,j,k)*1.e3
      end do
     end do
    end do
    name='QC'
    long_name='Cloud Water'
    units='g/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                   save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=qci(i,j,k)*1.e3
      end do
     end do
    end do
    name='QI'
    long_name='Cloud ice/Snow'
    units='g/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                   save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)

  end if

  if(doprecip) then

    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=qpl(i,j,k)*1.e3
      end do
     end do
    end do
    name='QR'
    long_name='Rain Water'
    units='g/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
    nfields1=nfields1+1
    do k=1,nzm
     do j=1,ny
      do i=1,nx
        tmp(i,j,k)=qpg(i,j,k)*1.e3
      end do
     end do
    end do
    name='QG'
    long_name='Graupel Water'
    units='g/kg'
    call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains,nfiles3D,1)
  end if

end subroutine micro_write_fields3D


!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use grid, only: wgty,nx_gl,ny_gl
  use vars, only : nstep,nprint,adz,dz,rho
  use terrain, only: terra
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)*wgty(j)*terra(i,j,k)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do
  total_water = total_water/dble(nx_gl*ny_gl)

end function total_water

! -------------------------------------------------------------------------------
! dummy effective radius functions:

function Get_reffc() ! liquid water
  real, pointer, dimension(:,:,:) :: Get_reffc
end function Get_reffc

function Get_reffi() ! ice
  real, pointer, dimension(:,:,:) :: Get_reffi
end function Get_reffi


end module microphysics



