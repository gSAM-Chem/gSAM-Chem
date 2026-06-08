#include "fppmacros"

! Module to compute coarsened fields

module stat_coars

use grid, only: nx_gl, ny_gl, nx, ny, nzm, ncoars_x, ncoars_y, &
                nstep,nstatcoars,nstatcoarsstart,nstatcoarsend,collect_coars,icycle
use microphysics, only: nmicro_fields
use params, only: dostatcoars, donudge3D

implicit none

integer nxm, nym ! subdomain's coarse grid sizes (maybe different for different subdomains)
integer nym_max ! maximum nym among all subdomains (for output)
integer j_start(ny), j_end(ny) ! start abd end indexes for coarsening (ny is just maximum possible index)

integer, parameter :: ntends = 2  ! number of tendency fields
integer, parameter :: nsfc = 5  ! maximum number of sfc fields
integer, parameter :: nflds = 11  ! maximum number of 3D prognostic fields
integer nprog  ! maximum number of 3D prognostic fields
real, allocatable :: fld_init(:,:,:,:)  ! fields at the begining of time step
real, allocatable :: fld_tend(:,:,:,:)  ! tendency due to process
real, allocatable :: fld_flux(:,:,:,:)  ! resolved flux (including precip related)
real, allocatable :: fld_flux_sgs(:,:,:,:)  ! SGS flux
real, allocatable :: fld_tend_nudge(:,:,:,:)  ! tendency due to 3D nudging

real(4), allocatable :: lon_c(:)  ! global coarse grid's lon if dolatlon = .true.
real(4), allocatable :: lat_c(:) ! global coarse grid's lat if dolatlon = .true.
real(4), allocatable :: y_c(:) ! global coarse grid's y if dolatlon = .true.

CONTAINS

!----------------------------------------------------------------------

subroutine coars_init()
  use grid
  use terrain, only: terra,terrau, terrav, terraw
  use vars
  use microphysics
  use params, only: docloud
  real size
  integer i,j,k,it,jt,ibuf(2)
  if(.not.dostatcoars) return
  call t_startf('coars')
  nprog = 4+nmicro_fields
  if(.not.allocated(fld_init)) then
      nxm = max(1,nx/ncoars_x)
      nym = 0
      if(abs(latitude(1,1)).lt.abs(latitude(1,ny))) then
        j = 1
        do while(j.le.ny)
          nym = nym+1
          j_start(nym) = j
          size = 0.
          do i=1,ncoars_y
            j_end(nym) = j
            size = size+ady(j) 
            j=j+1
            if(j.gt.ny.or.size.gt.ncoars_y) exit
          end do
        end do
      else
        j = ny
        do while(j.ge.1)
          nym = nym+1
          j_end(nym) = j
          size = 0.
          do i=1,ncoars_y
            j_start(nym) = j
            size = size+ady(j)
            j=j-1
            if(j.lt.1.or.size.gt.ncoars_y) exit
          end do
        end do
        j_start(1:nym) = j_start(nym:1:-1)
        j_end(1:nym) = j_end(nym:1:-1)
      end if
      ibuf(1) = nym
      call task_max_integer(ibuf(1),ibuf(2),1)
      nym_max = ibuf(2)
      allocate(lon_c(1:nxm),lat_c(1:nym_max),y_c(1:nym_max))
      allocate(fld_init(nxm,nym,nzm,nflds))
      allocate(fld_tend(nxm,nym,nzm,ntends))
      allocate(fld_flux(nxm,nym,nzm,nprog+4))
      allocate(fld_flux_sgs(nxm,nym,nzm,nprog+2))
      if(donudge3D) allocate(fld_tend_nudge(nxm,nym,nzm,3))
      ! lat/lon cooedinates of the coarse grid:
      call task_rank_to_index(rank,it,jt)
      do j=1,nym
        lat_c(j) = 0.5*(latv_gl(jt+j_end(j)+1)+latv_gl(jt+j_start(j)))
        y_c(j) = 0.5*(yv_gl(jt+j_end(j)+1)+yv_gl(jt+j_start(j)))
      end do
      do i=1,nxm
       lon_c(i) = 0.5*(lon_gl((i+it-1)*ncoars_x+1)+lon_gl((i+it)*ncoars_x))
      end do
  end if
  collect_coars = icycle.eq.1.and.mod(nstep,nstatcoars).eq.0.and.nstep.ge.nstatcoarsstart &
                                   .and.nstep.le.nstatcoarsend
  if(icycle.eq.1.and.collect_coars) then
    call coars_fld(u(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terrau(1:nx,1:ny,1:nzm),fld_init(:,:,:,1))
    call coars_fld(v(1:nx,1:ny,1:nzm),muv(1:ny),adyv(1:ny),terrav(1:nx,1:ny,1:nzm),fld_init(:,:,:,2))
    call coars_fld(w(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_init(:,:,:,3))
    call coars_fld(t(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,4))
    call coars_fld(tabs(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,5))
    call coars_fld(qv(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,6))
    call coars_fld(qcl(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,7))
    call coars_fld(qci(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,8))
    call coars_fld(qpl(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,9))
    call coars_fld(qpg(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,10))
    call coars_fld(qpi(1:nx,1:ny,1:nzm)-qpg(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),&
                   fld_init(:,:,:,11))
  end if
  fld_tend = 0.
  fld_flux = 0.
  fld_flux_sgs = 0.
  if(donudge3D) fld_tend_nudge = 0.

  call t_stopf('coars')
end subroutine coars_init

!----------------------------------------------------------------------

subroutine coars_fld(fld,mu,ady,terra,fld_coars)

! store the field before the process to compute tendency due to that process

real, intent(in) ::  fld(nx,ny,nzm) ! hi-res field
real, intent(in) ::  terra(nx,ny,nzm) ! terrain indicator
real(8), intent(in) ::  mu(ny) ! cos(lat)
real(8), intent(in) ::  ady(ny) ! grid-cell width in y
real, intent(out) ::  fld_coars(nxm,nym,nzm) ! coarsend field

integer i, j, k, ii, jj, starti, startj
real(8) factor, www
real(8) tmp

if(.not.dostatcoars) return

call t_startf ('coars')

fld_coars(:,:,:) = 0.

do k=1,nzm
 startj=1
 do jj=1,nym
  starti=1
  do ii=1,nxm
    factor=0.
    tmp = 0.
    do j=j_start(jj),j_end(jj)
        do i=starti,starti+(ncoars_x-1)
          www = mu(j)*ady(j)*terra(i,j,k)  
          factor = factor + www
          tmp = tmp+fld(i,j,k)*www  
        enddo ! i
    enddo ! j
    fld_coars(ii,jj,k)=tmp/(factor+1.e-10) 
    starti=starti+ncoars_x
  enddo  ! ii
  startj=startj+ncoars_y
 enddo ! jj
end do ! k

call t_stopf ('coars')

end subroutine coars_fld

!--------------------------------------------------------------------

subroutine coars_sfc(fld,mu,ady,mask,fld_coars)

! store the field before the process to compute tendency due to that process

real, intent(in) ::  fld(nx,ny) ! hi-res field
real(8), intent(in) ::  mu(ny) ! cos(lat)
real(8), intent(in) ::  ady(ny) ! grid-cell width in y
real, intent(in) ::  mask(nx,ny) ! conditional mask (1/0)
real, intent(out) ::  fld_coars(nxm,nym) ! coarsened field

integer i, j, ii, jj, starti, startj
real(8) factor, www
real(8) tmp

if(.not.dostatcoars) return

call t_startf ('coars')

fld_coars(:,:) = 0.

startj=1
do jj=1,nym
 starti=1
 do ii=1,nxm
   factor=0.
   tmp = 0.
   do j=j_start(jj),j_end(jj)
       do i=starti,starti+(ncoars_x-1)
         www = mu(j)*ady(j)*mask(i,j)
         factor = factor + www
         tmp = tmp+fld(i,j)*www
       enddo ! i
   enddo ! j
   fld_coars(ii,jj)=tmp/(factor+1.e-10)
   starti=starti+ncoars_x
 enddo  ! ii
 startj=startj+ncoars_y
enddo ! jj

call t_stopf ('coars')

end subroutine coars_sfc
!----------------------------------------------------------------------

subroutine coars_tend(fld,mu,ady,terra,nfld)

! compute tendency due to a process

real, intent(in) ::  fld(nx,ny,nzm) ! initial field
real, intent(in) ::  terra(nx,ny,nzm) ! terrain indicator
real(8), intent(in) ::  mu(ny) ! cos(lat)
real(8), intent(in) ::  ady(ny) ! grid-cell width in y
integer, intent(in) :: nfld ! field identifier

integer i, j, k, ii, jj, starti, startj
real(8) factor, www
real(8) tmp

if(.not.dostatcoars) return

call t_startf ('coars')

do k=1,nzm
 startj=1
 do jj=1,nym
  starti=1
  do ii=1,nxm
    factor=0.
    tmp = 0.
    do j=j_start(jj),j_end(jj)
        do i=starti,starti+(ncoars_x-1)
          www = mu(j)*ady(j)*terra(i,j,k)
          factor = factor + www
          tmp = tmp+fld(i,j,k)*www
        enddo ! i
    enddo ! j
    fld_tend(ii,jj,k,nfld) = tmp/(factor+1.e-10)-fld_tend(ii,jj,k,nfld)
    starti=starti+ncoars_x
  enddo  ! ii
  startj=startj+ncoars_y
 enddo ! jj
end do ! k

call t_stopf ('coars')

end subroutine coars_tend

!----------------------------------------------------------------------


subroutine write_coars3D()

use grid
use vars, only: u,v,w,t,rho,rhow
use params, only: docloud, dolatlon
use microphysics
use terrain, only: terra,terrau,terrav,terraw
implicit none
integer, external :: lenstr

character *120 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *10 filetype
character *10 units
character *5 sepchar
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real(4) tmp(nxm,nym_max,nzm),lats(nym_max*nsubdomains),lons(nxm*nsubdomains)
real(4) ys(nym_max*nsubdomains)
integer nyms(nsubdomains),ibuf(1)

if(.not.dostatcoars) return

call t_startf('write_coars')

nfields=21
if(docloud) nfields=nfields+nmicro_fields*2+8
if(donudge3D) nfields=nfields+3

ibuf(1) = nym
call task_gatherv_integer(ibuf(1),nyms,1)
call task_gatherv_float4(lat_c,lats,nym_max)
call task_gatherv_float4(y_c,ys,nym_max)
call task_gatherv_float4(lon_c,lons,nxm)

if(rank.eq.rank-mod(rank,nsubdomains/nfilescoars)) then

  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do
  if(nfilescoars.eq.1) then
     sepchar=""
  else
     write(rankchar,'(i4)') rank/(nsubdomains/nfilescoars)
     sepchar="_"//rankchar(5-lenstr(rankchar):4)
  end if

    filetype = '.coarse.3D'
    filename='./OUT_MOMENTS/'//trim(case)//'_coars_'//trim(caseid)//'_'// &
        trim(date_pr)//'_'//timechar(1:10)//filetype//sepchar
    open(46,file=filename,status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')

  if(masterproc) then

     write(46) nstep, time
     write(46) nxm,nym_max,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields,nfilescoars
   !  print*,nxm,nym,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields,nfilescoars
     write(46) real(dx*float(ncoars_x),4), real(dy*float(ncoars_y),4)
     write(46) real(float(nstep)*dt/(3600.*24.)+day0,4),datechar
     write(46) dolatlon
     write(46) 'atm'
     write(46) real(z(1:nzm),4),real(pres(1:nzm),4), &
                real(lats,4),real(lons,4),&
                real(ys,4),real(dz*adz(1:nzm),4),real(rho(1:nzm),4),real(rhow(1:nzm),4)

   end if ! msterproc 
end if 

nfields1=0

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_init(i,j,k,1)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='U'
long_name='Zonal wind (before step)'
units='m/s'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_init(i,j,k,2)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='V'
long_name='Meridional wind (before step)'
units='m/s'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,1)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_init(i,j,k,3)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='W'
long_name='Vertical wind (before step)'
units='m/s'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,2)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_init(i,j,k,4)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='T'
long_name='Liquid/ice static energy (before step)'
units='K'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_init(i,j,k,5)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='TABS'
long_name='Absolute Temperature(before step)'
units='K'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux_sgs(i,j,k,4+nmicro_fields+1)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='KEDDY'
long_name='Eddy Diffusion Coefficient for momentum'
units='m2/s'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux_sgs(i,j,k,4+nmicro_fields+2)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='KEDDYSC'
long_name='Eddy Diffusion Coefficient for scalars'
units='m2/s'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)


nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux(i,j,k,1)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoUW'
long_name='Verical Flux of Zonal wind (Resolved)'
units='N/m2'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,2)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux_sgs(i,j,k,1)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoUWS'
long_name='Verical Flux of Zonal wind (SGS)'
units='N/m2'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,2)


nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux(i,j,k,2)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoVW'
long_name='Vertical Flux of Meridional wind (Resolved)'
units='N/m2'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux_sgs(i,j,k,2)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoVWS'
long_name='Vertical Flux of Meridional wind (SGS)'
units='N/m2'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)


nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux(i,j,k,3)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoWW'
long_name='Vertical Flux of Vertical wind (Resolved)'
units='N/m2'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux_sgs(i,j,k,3)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoWWS'
long_name='Vertical Flux of Vertical wind (SGS)'
units='N/m2'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux(i,j,k,4)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoTW'
long_name='Vertical Flux of Liquid/ice static energy (Resolved)'
units='K kg/m2/s'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,2)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_flux_sgs(i,j,k,4)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='rhoTWS'
long_name='Vertical Flux of Liquid/ice static energy (SGS)'
units='K kg/m2/s'
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,2)

nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
     do i=1,nxm
        tmp(i,j,k)=fld_tend(i,j,k,2)/dtn
     end do
  end do
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
end do
name="QRAD"
long_name='Tendency of T due to radiation'
units="K/s"
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi,&
               rank,nsubdomains,nfilescoars,1)


if(docloud) then

      nfields1=nfields1+1
      do k=1,nzm
       do j=1,nym
         do i=1,nxm
            tmp(i,j,k)=fld_init(i,j,k,6)
         enddo
       enddo
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      enddo
      name='QV'
      long_name='Water Vapor'
      units='kg/kg'
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

      nfields1=nfields1+1
      do k=1,nzm
       do j=1,nym
         do i=1,nxm
            tmp(i,j,k)=fld_init(i,j,k,7)
         enddo
       enddo
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      enddo
      name='QC'
      long_name='Cloud Water'
      units='kg/kg'
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

      nfields1=nfields1+1
      do k=1,nzm
       do j=1,nym
         do i=1,nxm
            tmp(i,j,k)=fld_init(i,j,k,8)
         enddo
       enddo
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      enddo
      name='QI'
      long_name='Cloud Ice'
      units='kg/kg'
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

      nfields1=nfields1+1
      do k=1,nzm
       do j=1,nym
         do i=1,nxm
            tmp(i,j,k)=fld_init(i,j,k,9)
         enddo
       enddo
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      enddo
      name='QR'
      long_name='Rain'
      units='kg/kg'
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

      nfields1=nfields1+1
      do k=1,nzm
       do j=1,nym
         do i=1,nxm
            tmp(i,j,k)=fld_init(i,j,k,10)
         enddo
       enddo
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      enddo
      name='QG'
      long_name='Graupel'
      units='kg/kg'
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

      nfields1=nfields1+1
      do k=1,nzm
       do j=1,nym
         do i=1,nxm
            tmp(i,j,k)=fld_init(i,j,k,11)
         enddo
       enddo
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      enddo
      name='QS'
      long_name='Snow'
      units='kg/kg'
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                rank,nsubdomains,nfilescoars,1)

 do n = 1,nmicro_fields

      nfields1=nfields1+1
      do k=1,nzm
         do j=1,nym
            do i=1,nxm
               tmp(i,j,k)=fld_flux(i,j,k,4+n)
            end do
         end do
         tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      end do
      name="rho"//TRIM(mkname(n))//"W"
      long_name="Vertical Flux of "//TRIM(mklongname(n))//' (Resolved)'
      units="kg/m2/s "//TRIM(mkunits(n))
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)

      nfields1=nfields1+1
      do k=1,nzm
         do j=1,nym
            do i=1,nxm
               tmp(i,j,k)=fld_flux_sgs(i,j,k,4+n)
            end do
         end do
         tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      end do
      name="rho"//TRIM(mkname(n))//"WS"
      long_name="Vertical Flux of "//TRIM(mklongname(n))//' (SGS)'
      units="kg/m2/s "//TRIM(mkunits(n))
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)


 end do

      nfields1=nfields1+1
      do k=1,nzm
         do j=1,nym
            do i=1,nxm
               tmp(i,j,k)=fld_flux(i,j,k,4+nmicro_fields+1)
            end do
         end do
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      end do
      name="PREC"
      long_name="Precipitation Flux"
      units="kg/m2/s "
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)
      nfields1=nfields1+1
      do k=1,nzm
         do j=1,nym
            do i=1,nxm
               tmp(i,j,k)=fld_flux(i,j,k,4+nmicro_fields+2)
            end do
         end do
       tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      end do
      name="LPREC"
      long_name="Latent Heat Flux due to Precipitation"
      units="W/m2"
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)
      nfields1=nfields1+1
      do k=1,nzm
         do j=1,nym
            do i=1,nxm
               tmp(i,j,k)=fld_flux(i,j,k,4+nmicro_fields+3)
            end do
         end do
         tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      end do
      name="SED"
      long_name="Sedimentation Flux"
      units="kg/m2/s "
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)
      nfields1=nfields1+1
      do k=1,nzm
         do j=1,nym
            do i=1,nxm
               tmp(i,j,k)=fld_flux(i,j,k,4+nmicro_fields+4)
            end do
         end do
         tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      end do
      name="LSED"
      long_name="Latent Heat Flux due to Sedimentation"
      units="W/m2"
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
               rank,nsubdomains,nfilescoars,2)


      nfields1=nfields1+1
      do k=1,nzm
         do j=1,nym
            do i=1,nxm
               tmp(i,j,k)=fld_tend(i,j,k,1)/dtn
            end do
         end do
         tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
      end do
      name="QP_MICRO"
      long_name='Tendency of QP due to microphysics'
      units="kg/kg/s"
      call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi,&
                      rank,nsubdomains,nfilescoars,1)

end if


if(donudge3D) then

  nfields1=nfields1+1
  do k=1,nzm
    do j=1,nym
      do i=1,nxm
        tmp(i,j,k)=fld_tend_nudge(i,j,k,1)*dtfactor*86400.
      enddo
    enddo
    tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
  enddo
  name='U_NUDG'
  long_name='Nudging tendency of zonal wind'
  units='m/s/day'
  call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                  rank,nsubdomains,nfilescoars,1)

  nfields1=nfields1+1
  do k=1,nzm
    do j=1,nym
      do i=1,nxm
        tmp(i,j,k)=fld_tend_nudge(i,j,k,2)*dtfactor*86400.
      enddo
    enddo
    tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
  enddo
  name='V_NUDG'
  long_name='Nudging tendency of meridional wind'
  units='m/s/day'
  call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi, &
                  rank,nsubdomains,nfilescoars,1)

  nfields1=nfields1+1
  do k=1,nzm
    do j=1,nym
      do i=1,nxm
        tmp(i,j,k)=fld_tend_nudge(i,j,k,3)*dtfactor*86400.
      enddo
    enddo
    tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
  enddo
  name='T_NUDG'
  long_name='Nudging tendency of liquid/ice static energy'
  units='K/day'
  call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi,&
                  rank,nsubdomains,nfilescoars,1)

end if

call coars_fld(terra(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terra(1:nx,1:ny,1:nzm),fld_init(:,:,:,1))
nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_init(i,j,k,1)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='TERRA'
long_name='Terrain Mask for Scalars'
units=''
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi,&
                rank,nsubdomains,nfilescoars,1)

call coars_fld(terraw(1:nx,1:ny,1:nzm),mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_init(:,:,:,1))
nfields1=nfields1+1
do k=1,nzm
  do j=1,nym
    do i=1,nxm
      tmp(i,j,k)=fld_init(i,j,k,1)
    enddo
  enddo
  tmp(:,nym+1:nym_max,k) = tmp(1,1,k)
enddo
name='TERRAW'
long_name='Terrain Mask for Fuxes and Vert. Vel'
units=''
call compress3D(tmp,nxm,nym_max,nzm,name,long_name,units,savecoarsbin,dompi,&
                rank,nsubdomains,nfilescoars,2)

call task_barrier()

if (nfields .ne. nfields1) then
  if(masterproc) print*,'write_coars error: nfields=',nfields,' nfields1=',nfields1
  call task_abort()
endif
if (masterproc) then
  close(46)
endif

if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_coars error: nfields'
    call task_abort()
end if
if(masterproc) then
     print*, 'Writting tendency data. file:'//filename
endif

call t_stopf('write_coars')

end subroutine write_coars3D


!=================================================================
!=================================================================
!=================================================================

subroutine write_coars2D

use vars
use params
use terrain, only: k_terra,elevationg
use slm_vars, only: landtype,landicemask,seaicemask,vege_YES,soilt,soilw,t_canop
implicit none
character *120 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *10 filetype
character *10 units
integer i,j,k,nfields,nfields1
real tmph(nx,ny),unity(nx,ny),tmp1(nxm,nym)
real(4) tmp(nxm,nym_max),lats(nym_max*nsubdomains),lons(nxm*nsubdomains)
real(4) ys(nym_max*nsubdomains)
integer nyms(nsubdomains)
real coef, coefa
integer, external :: lenstr
integer ibuf(1)

nfields= 13
if(OCEAN.or..not.SLM)  nfields = nfields - 8
ibuf(1) = nym
call task_gatherv_integer(ibuf(1),nyms,1)
call task_gatherv_float4(lat_c,lats,nym_max)
call task_gatherv_float4(y_c,ys,nym_max)
call task_gatherv_float4(lon_c,lons,nxm)

if(masterproc) then

  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do i=1,11-lenstr(timechar)-1
    timechar(i:i)='0'
  end do

! Make sure that the new run doesn't overwrite the file from the old run

    filetype = '.coarse.2D'

    filename='./OUT_MOMENTS/'//trim(case)//'_coars_'//trim(caseid)//'_'// &
          rankchar(5-lenstr(rankchar):4)//trim(date_pr)//'_'//timechar(1:10)//filetype
    open(46,file=filename,status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')

    write(46) nstep, dolatlon
    write(46) 'atm'
    write(46) nxm,nym_max,nzm,nsubdomains, nsubdomains_x,nsubdomains_y,nfields
    write(46) real(dx*float(ncoars_x),4)
    write(46) real(dy*float(ncoars_y),4)
    write(46) real(nstep*dt/(3600.*24.)+day0,4)

end if! masterproc
   unity = 1.
   nfields1 = 0

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      tmph(i,j)=fluxbt(i,j)*rhow(k)
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='rhoWT'
   long_name='Surface Flux of Liquid/ice static energy '
   units='K kg/m2/s'
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      tmph(i,j)=fluxbq(i,j)*rhow(k)
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='rhoWQ'
   long_name='Surface Flux of Total Water'
   units='kg/m2/s'
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      tmph(i,j)=fluxbu(i,j)*rhow(k)
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='rhoWU'
   long_name='Surface Zonal Stress'
   units='kg/m/s'
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      tmph(i,j)=fluxbv(i,j)*rhow(k)
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='rhoWV'
   long_name='Surface Meridional Stress'
   units='kg/m/s'
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      tmph(i,j)=sstxy(i,j)+t00
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='SKT'
   long_name='Surface Skin Temperature'
   units='K'              
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   if(.not.OCEAN.and.SLM) then

   nfields1=nfields1+1
   call coars_sfc(real(landmask),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='LANDMASK'
   long_name='Land Fraction'
   units=''
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   call coars_sfc(real(1.-landmask),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='SEAMASK'
   long_name='Sea Fraction'
   units=''
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)


   nfields1=nfields1+1
   call coars_sfc(real(landicemask),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='LANDICEMASK'
   long_name='Land Ice Fraction'
   units=''
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   call coars_sfc(real(seaicemask),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='SEAICEMASK'
   long_name='Land Ice Fraction'
   units=''
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   call coars_sfc(real(vege_YES),mu(1:ny),ady(1:ny),unity(1:nx,1:ny),tmp1(:,:))
   name='VEG'
   long_name='Vegetated Land Fraction'
   units=''
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmph(i,j)=t_canop(i,j)
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),real(vege_YES),tmp1(:,:))
   name='TCANOP'
   long_name='Canopy Temperatire'
   units='K'
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmph(i,j)=soilt(1,i,j)
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),real(landmask),tmp1(:,:))
   name='SOILT'
   long_name='Soil Temperature'
   units='K'
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmph(i,j)=soilw(1,i,j)
    end do
   end do
   call coars_sfc(tmph(1:nx,1:ny),mu(1:ny),ady(1:ny),real(landmask),tmp1(:,:))
   name='SOILW'
   long_name='Soil Water Content'
   units='m3/m3'
   tmp(1:nxm,1:nym) = tmp1(1:nxm,1:nym)
   tmp(1:nxm,nym+1:nym_max) = tmp(1,1)
   call compress2DC(tmp,nxm,nym_max,name,long_name,units, &
                               savecoarsbin,dompi,rank,nsubdomains)

end if ! .not.OCEAN

end subroutine write_coars2D

!==========================================================================================

subroutine compress2DC (f,nx,ny,name, long_name, units, &
                       savebin, dompi, rank, nsubdomains)


use grid, only: masterproc,nsubdomains_x,nsubdomains_y
implicit none
! Input:

integer nx,ny
real(4) f(nx,ny)
character(*) name,long_name,units
integer rank,irank,nsubdomains
logical savebin, dompi

! Local:

integer(2), allocatable :: byte(:)
real(4), allocatable :: byte4(:)
integer(2), allocatable :: byte_2D(:)
real(4), allocatable :: byte4_2D(:)

character(7) form
integer integer_max, integer_min
parameter (integer_min=-32000, integer_max=32000)
!       parameter (integer_min=-127, integer_max=127)
real(4) f_max(1),f_min(1), f_max1(1), f_min1(1), scale
integer i,j,m,it,jt,rrr,ttt,nx_glm, ny_glm
integer n_in, count, reqs_in(max(1,nsubdomains-1)), tag
logical flag(nsubdomains-1)

if(masterproc) then
    write(46) name,' ',long_name,' ',units
    write(46) savebin
end if

nx_glm = nx*nsubdomains_x
ny_glm = ny*nsubdomains_y

if(savebin) then

  allocate (byte4(nx*ny))
  allocate(byte4_2D(nx_glm*ny_glm))

   count = 0
   do j=1,ny
    do i=1,nx
     count = count+1
     byte4(count) = f(i,j)
    end do
   end do

  if(.not.dompi) then
    write(46) byte4(:)
  else
    call task_gatherv_float4(byte4,byte4_2D,nx*ny)
    if(masterproc) write(46) byte4_2D(:)
  end if

  deallocate(byte4)
  deallocate(byte4_2D)


else

   allocate (byte(nx*ny))
   allocate(byte_2D(nx_glm*ny_glm))

   f_max=-1.e30
   f_min= 1.e30
   do j=1,ny
    do i=1,nx
     f_max(1) = max(f_max(1),f(i,j))
     f_min(1) = min(f_min(1),f(i,j))
    end do
   end do
   if(dompi) then
     f_max1=f_max
     f_min1=f_min
     call task_max_real4(f_max1,f_max,1)
     call task_min_real4(f_min1,f_min,1)
   endif

   if(abs(f_max(1)).lt.10..and.abs(f_min(1)).lt.10.) then
          form='(f10.7)'
   else if(abs(f_max(1)).lt.100..and.abs(f_min(1)).lt.100.) then
          form='(f10.6)'
   else if(abs(f_max(1)).lt.1000..and.abs(f_min(1)).lt.1000.) then
          form='(f10.5)'
   else if(abs(f_max(1)).lt.10000..and.abs(f_min(1)).lt.10000.) then
          form='(f10.4)'
   else if(abs(f_max(1)).lt.100000..and.abs(f_min(1)).lt.100000.) then
          form='(f10.3)'
   else if(abs(f_max(1)).lt.1000000..and.abs(f_min(1)).lt.1000000.) then
          form='(f10.2)'
   else if(abs(f_max(1)).lt.10000000..and.abs(f_min(1)).lt.10000000.) then
          form='(f10.1)'
   else if(abs(f_max(1)).lt.100000000..and.abs(f_min(1)).lt.100000000.) then
          form='(f10.0)'
   else
          form='(f10.0)'
!          f_min=-999.
!          f_max= 999.
   end if


   if(f_max(1)-f_min(1).lt.1.e-10) then
      scale = 0.
   else
      scale = float(integer_max-integer_min)/(f_max(1)-f_min(1))
   end if

   count = 0
   do j=1,ny
    do i=1,nx
      count = count+1
      byte(count)= integer_min+scale*(f(i,j)-f_min(1))
    end do
   end do

  if(masterproc) then
     write(46) real(f_max(1),4),real(f_min(1),4)
    ! print*,'writing ',trim(long_name)
  end if

  if(.not.dompi) then
    write(46) byte(:)
  else
    call task_gatherv_integer2(byte,byte_2D,nx*ny)
    if(masterproc) write(46) byte_2D(:)
  end if

  deallocate(byte)
  deallocate(byte_2D)

end if ! savebin


call task_barrier()

end subroutine compress2DC



end module stat_coars






