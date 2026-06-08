#include "./../fppmacros"
!=====================================================================
! Based on write_fields2D.f90 of SAM
! Write 2D variables into the output file
!	
!=====================================================================
subroutine slm_write2D

use grid
use slm_vars
use params, only: dolatlon
use pnetcdf_stuff, only: open_file2D_pnetcdf, close_file_pnetcdf, dopnetcdf

implicit none

!!!!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!
!
! MAKE SURE THE NUMBER OF DIAGNOSTICS HERE MATCH THOSE IN
! THE CALLING ROUTINE, write_fields2D.f90
!
!!!!!!!!!!!!!!!!!!!!!!!CAUTION!!!!!!!!!!!!!!!!!!!!

! argument list
character(200)filename
character(80) long_name
character(8)  name
character(10) timechar
character(7)  filetype
character(10) units
character(4)  numchar
integer i,j,k,nfields,nfields1,nsteplast
real(4) tmp(nx,ny)
character*7 filestatus
real(DBL) coef_dble
integer, external :: lenstr
logical flag

dopnetcdf = save2DLnetcdf

nfields = (42+(nsoil*2))
if(dosoiltnudging) nfields = nfields + nsoil
if(dosoilwnudging) nfields = nfields + nsoil
if(dooutthermal) nfields = nfields + 2*nsoil

nfields1 = 0 

if(.not.dopnetcdf) then

 if(masterproc) then

  write(timechar,'(i10)') nstep
  do i=1,11-lenstr(timechar)-1
    timechar(i:i)='0'
  end do

! Make sure that the new run doesn't overwrite the file from the old run

    if(.not.save2DLsep.and.notopened2DL.and.(nrestart.eq.0.or.nrestart.eq.2)) then
      filestatus='new'
    else
      filestatus='unknown'
    end if

    filetype = '.2D'

    if(save2DLsep) then
       filename='./OUT_2DL/'//trim(case)//'_'//trim(caseid)//'_'// &
          trim(date_pr)//'_'//timechar(1:10)//filetype
          open(46,file=filename,status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')
       print*, 'Writting to file: '//trim(filename)
    else
       filename='./OUT_2DL/'//trim(case)//'_'//trim(caseid)//filetype
       open(46,file=filename,status=filestatus,form='unformatted',BUFFEREDYES ACTION='READWRITE')
       print*, 'Writting to file: '//trim(filename)
       do while(.true.)
         read(46,end=222)  nsteplast
         if(nsteplast.ge.nstep) then
           backspace(46)
           goto 222   ! yeh, I know, it's bad ....
         end if
         read(46)
         read(46)
         read(46)
         read(46)
         read(46)
         do i=1,nfields
           read(46)
           read(46) flag
           if(.not.flag) read(46)
           read(46)
         end do
       end do
 222    continue
       print*,'nsteplast=',nsteplast
       notopened2DL=.false.
    end if

    write(46) nstep, dolatlon, time, datechar, timeUTsec
    write(46) 'lnd'
    write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    write(46) real(dx,4),lat_gl,latv_gl,real(y_gl,4),real(yv_gl,4)
    write(46) real(dy,4),lon_gl,lonu_gl
    write(46) dble(nstep)*dt/(3600._8*24._8)+day0

 end if! masterproc

else

! open netcdf file:

 write(timechar,'(i10)') nstep
 do k=1,11-lenstr(timechar)-1
   timechar(k:k)='0'
 end do

 filename ='./OUT_2DL/'//trim(case)//'_'//trim(caseid)//'_'// &
       trim(date_pr)//'_'//timechar(1:10)//'.2D_lnd.nc'

 call open_file2D_pnetcdf(filename)

 if(masterproc) print*, 'Writting to file: '//trim(filename)

end if ! .not.dopnetcdf



if(.not.nstep.eq.1.and.save2DLavg) then
   coef_dble = 1./float(nsave2DL)
else
   coef_dble = 1.
end if



k = 0   
!====================================================================
!! Surface momentum flux in x
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_taux_sfc(i,j)*coef_dble,4)
     end do
   end do
  name='TAUX'
  long_name='SFC Momentum Flux X'
  units='m2/s2'
 call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
!! Surface momentum flux in y
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_tauy_sfc(i,j)*coef_dble,4)
     end do
   end do
  name='TAUY'
  long_name='SFC Momentum Flux Y'
  units='m2/s2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! precipitation interception rate by canopy 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_precip(i,j)*86400._DBL*coef_dble,4)
     end do
   end do
  name='PCNP'
  long_name='Precip. Interception Rate'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! precipitation drainage rate from canopy
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_drain(i,j)*86400._DBL*coef_dble,4)
     end do
   end do
  name='DCNP'
  long_name='Precip. Drainage Rate'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! precip rate at surface 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_precip_sfc(i,j)*86400._DBL*coef_dble,4)
     end do
   end do
  name='PSFC'
  long_name='Precip. Rate at soil surface'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! net  rad  on canopy 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_net_rad(1,i,j)*coef_dble,4)
     end do
   end do
  name='NetRC'
  long_name='Net Rad on canopy'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! net  rad  on soil 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_net_rad(2,i,j)*coef_dble,4)
     end do
   end do
  name='NetRS'
  long_name='Net Rad on soil'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!====================================================================
! Friction Velocity 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_ustar(i,j)*coef_dble,4)
     end do
   end do
  name='USTAR'
  long_name='Friction Velocity'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! LAI
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_LAI(i,j)*coef_dble,4)
     end do
   end do
  name='LAI'
  long_name='Leaf Area Index'
  units='m2/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!====================================================================
! ref. resistance, Ra 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_r_a(i,j)*coef_dble,4)
     end do
   end do
  name='Ra'
  long_name='Aerodynamic Resistance'
  units='s/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! leaf boundary resistance 
!========================
  nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_r_b(i,j)*coef_dble,4)
     end do
   end do
  where(tmp(:,:).gt.5000.) tmp(:,:) = 5000.
  name='Rb'
  long_name='Leaf boundary Resistance'
  units='s/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! stomatal resistance 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_r_c(i,j)*coef_dble,4)
     end do
   end do
  where(tmp(:,:).gt.5000.) tmp(:,:) = 5000.
  name='Rc'
  long_name='Stomatal Resistance'
  units='s/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! undercanopy resistance 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_r_d(i,j)*coef_dble,4)
     end do
   end do
  where(tmp(:,:).gt.5000.) tmp(:,:) = 5000.
  name='Rd'
  long_name='Undercanopy Resistance'
  units='s/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil surface resistance 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_r_soil(i,j)*coef_dble,4)
     end do
   end do
  where(tmp(:,:).gt.5000.) tmp(:,:) = 5000.
  name='Rsoil'
  long_name='Soil Resistance'
  units='s/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! canopy sensible heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_shf_canop(i,j)*coef_dble,4)
     end do
   end do
  name='SHFC'
  long_name='canopy sensible heat flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil sensible heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_shf_soil(i,j)*coef_dble,4)
     end do
   end do
  name='SHFS'
  long_name='soil sensible heat flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! canopy latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_lhf_canop(i,j)*coef_dble,4)
     end do
   end do
  name='LHFC'
  long_name='canopy latent heat flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_lhf_soil(i,j)*coef_dble,4)
     end do
   end do
  name='LHFS'
  long_name='soil latent heat flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! wet latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_lhf_wet(i,j)*coef_dble,4)
     end do
   end do
  name='LHFCW'
  long_name='canopy wet latent heat flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! canopy evaporation flux
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_evp_canop(i,j)*coef_dble,4)*86400.
     end do
   end do
  name='EVPC'
  long_name='canopy evaporation flux'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil evaporation flux
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_evp_soil(i,j)*coef_dble,4)*86400.
     end do
   end do
  name='EVPS'
  long_name='soil evaporation flux'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! wet evaporation flux
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_evp_wet(i,j)*coef_dble,4)*86400.
     end do
   end do
  name='EVPCW'
  long_name='canopy wet evaporation flux'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! soil surface wetness
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_fh(i,j)*coef_dble,4)
     end do
   end do
  name='FH'
  long_name='soil surface wetness (0-1)'
  units=' '
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!====================================================================
! Snow Cover
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       if(snow_mass(i,j).gt.0.) then
         tmp(i,j) = 1.
       else
         tmp(i,j) = 0.
       end if
     end do
   end do
  name='SNOW'
  long_name='Snow Cover'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!====================================================================
! Snow Depth
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=snow_mass(i,j)/rho_snow
     end do
   end do
  name='SNOWD'
  long_name='Snow_Depth (Inst)'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! Snow Temperature
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=snowt(i,j)
     end do
   end do
  name='SNOWT'
  long_name='Snow Temperature (Inst)'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! dry latent heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_lhf_tr(i,j)*coef_dble,4)
     end do
   end do
  name='LHFCD'
  long_name='canopy dry latent heat flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! dry evaporation flux
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_evp_tr(i,j)*coef_dble,4)*86400.
     end do
   end do
  name='EVPCD'
  long_name='canopy dry evaporation flux'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! ground heat flux 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_grflux0(i,j)*coef_dble,4)
     end do
   end do
  name='GRFLUX0'
  long_name='ground heat flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! Precipitation infiltrate rate 
!========================

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_precip_in(i,j)*86400._DBL*coef_dble,4)
     end do
   end do
  name='PIN'
  long_name='Precipitation infiltrate rate'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! Precipitation at reference level
!========================

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_precip_ref(i,j)*86400._DBL*coef_dble,4)
     end do
   end do
  name='PREF'
  long_name='Precipitation at reference level'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! surface run off rate 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_run_off_sfc(i,j)*86400._DBL*coef_dble,4)
     end do
   end do
  name='ROFFS'
  long_name='surface runoff rate'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! canopy temperature 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_tcnp(i,j)*coef_dble,4)
     end do
   end do
  name='TCNP'
  long_name='canopy temperature'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! Vegetation Mask
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(vege_YES(i,j),4)
     end do
   end do
  name='VEGMASK'
  long_name='Vegetation Mask'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! CAS temperature 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_tcas(i,j)*coef_dble,4)
     end do
   end do
  name='TCAS'
  long_name='canopy air space temperature'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! CAS q 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_qcas(i,j)*1000._DBL*coef_dble,4)
     end do
   end do
  name='QCAS'
  long_name='canopy air space q'
  units='g/kg'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! canopy water storage 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_mw(i,j)*coef_dble,4)
     end do
   end do
  name='MW'
  long_name='canopy water storage'
  units='kg/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! drainage rate 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_drainage(i,j)*86400._DBL*coef_dble,4)
     end do
   end do

  name='WDRAIN'
  long_name='water drainage rate from bottom soil layer'
  units='mm/d'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! puddle water / flood  
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_mws(i,j)*coef_dble,4)
     end do
   end do
  name='MWS'
  long_name='puddle water storage'
  units='mm'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! net upward SW to ref.level 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_net_swup(1,i,j)*coef_dble,4)
     end do
   end do
  name='SUP'
  long_name='SW upward to ref.level'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! net upward LW to ref. level
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_net_lwup(1,i,j)*coef_dble)
     end do
   end do
  name='LUP'
  long_name='LW upward to ref. level'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! soil total water
!========================
  nfields1=nfields1+1
  do j=1,ny
   do i=1,nx
     tmp(i,j) = 0.
     do k = 1,nsoil
       tmp(i,j)=tmp(i,j)+real(s_soilw(k,i,j)*s_depth(k,i,j)*coef_dble,4)
     end do
     if(landmask(i,j).eq.1) then
      tmp(i,j) = tmp(i,j)*rho_water
     else
      tmp(i,j) = tmp(i,j)*rho_ice
     end if
   end do
  end do
  name='TSW'
  long_name='Total soil/ice water content'
  units='kg'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

do k = 1,nsoil
write(numchar,'(i0)') k
!========================
! soil wetness 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_soilw(k,i,j)*poro_soil(k,i,j)*coef_dble,4)
     end do
   end do
  name='WS'//numchar
  long_name='soil water content'//numchar
  units='m3/m3 '
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil temperature
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_soilt(k,i,j)*coef_dble,4)
     end do
   end do
  name='TS'//numchar
  long_name='Soil temperature'//numchar
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

if(dooutthermal) then
!========================
! soil thermal conductivity
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_scond(k,i,j)*coef_dble,4)
     end do
   end do
  name='SCOND'//numchar
  long_name='Soil thermal conductivity'//numchar
  units='W/m/K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil heat capacity
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_scapa(k,i,j)*coef_dble,4)
     end do
   end do
  name='SCAPA'//numchar
  long_name='Soil heat capacity'//numchar
  units='J/kg/K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
end if

if(dosoilwnudging) then
!========================
! Soil moisture nudging
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_soilw_nudge(k,i,j)*coef_dble,4)
     end do
   end do
  name='SWNDG'//numchar
  long_name='SW nudging'//numchar
  units='/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
end if
if(dosoiltnudging) then
!========================
! Soil temperature nudging
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_soilt_nudge(k,i,j)*coef_dble,4)
     end do
   end do
  name='TSNDG'//numchar
  long_name='TS nudging'//numchar
  units='/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
end if
end do ! k


call task_barrier()

if(nfields.ne.nfields1) then
  if(masterproc) print*,'slm_fields2D: error in nfields!!',' nfields=',nfields,'nfields1=',nfields1
  call task_abort()
end if

if(dopnetcdf) then
     call close_file_pnetcdf()
     dopnetcdf = .false.
else
 if(masterproc) then
     close(46)
     if(save2DLsep.and.dogzip2D) call systemf('gzip -f '//filename)
 endif
end if
if(masterproc) print*, 'Done.', nfields1,'fields'


end subroutine
