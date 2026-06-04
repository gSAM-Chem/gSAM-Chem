#include "fppmacros"
     
subroutine write_fields2D
	
use vars
use params
use terrain, only: k_terra,elevationg,elevation
use slm_vars, only: landtype,landicemask,seaicemask,icemask
use tracers, only: ntracers, tr_xy, tr_ac_xy, tr_acs_xy, tracerunits, tracername
use buildings, only: doshadows, doirgroundfromwalls, &
                     doirwallsfromground, doirwallsfromwalls, &
                     shadow_mask, shadow_height, bld_mask, facemask, cell_sunny, &
                     ind_wall_sfc, flx_wall_sfc, sangle_wall_sfc, bld_index, doAC, doHEAT, &
                     t_bld, building_cooling_energy, building_heating_energy
use pnetcdf_stuff, only: open_file2D_pnetcdf, close_file_pnetcdf, dopnetcdf
use sgs, only: dodns

implicit none

character(200)filename
character(80) long_name
character(8)  name
character(10) timechar
character(3)  filetype
character(10) units
integer i,j,k,nfields,nfields1,nsteplast,itr
real(4) tmp(nx,ny)
real cape(nx,ny), cin(nx,ny)
real coef, coefa, factor
integer, external :: lenstr
real, external :: heat_index, qsatw
character(9) fileaction
logical flag

if(docheck) return

dopnetcdf = save2Dnetcdf

nfields= 20
if(docloud) then
    nfields = nfields+16
    if(save2Davg) nfields=nfields+1
end if
if(doprecip) then
    nfields = nfields+7
    if(snow2Dout) nfields = nfields+3
end if
if(dolongwave) then
    nfields = nfields+6
    if(save2Dradac) nfields = nfields+3
    if(.not.save2Davg.and.save2Drada) nfields = nfields+6
end if
if(doshortwave) then
    nfields = nfields+7
    if(save2Dradac) nfields = nfields+3
    if(.not.save2Davg.and.save2Drada) nfields = nfields+7
end if
if(doslabocean) nfields=nfields+1
if((ocean_type.ne.0.or.doslabocean.or.doequilocean.or.dodynamicocean).and..not.dossthomo.or. &
    (ISLAND.or.LAND).and.SLM.or.readsst) nfields=nfields+1
if(ISLAND.or.doseaice) nfields=nfields+1
if(doseaice.and.doseaiceevol) nfields = nfields+1
if(.not.SFC_FLX_FXD) Then
  nfields = nfields+13
  if(dodns) nfields = nfields+4 
end if
if(dotracers) nfields=nfields+2*ntracers
if(dotracers.and.dotrsfcflux) nfields=nfields+ntracers
if(doterrain) nfields=nfields+1
if(doterrain.and.dobuildings) nfields=nfields+3
if(doterrain.and.dobuildings.and..not.(doAC.and.doHEAT)) nfields=nfields+1
if(doterrain.and.dobuildings.and.(doAC.or.doHEAT)) nfields=nfields+1
if(doterrain.and.dobuildings.and.doshadows) nfields=nfields+3
if(doterrain.and.dobuildings.and.doirgroundfromwalls) nfields=nfields+3
if(donudge3D.or.doregion) nfields=nfields+8
!===================================================================

nfields1=0

if(.not.dopnetcdf) then

 if(masterproc) then

  write(timechar,'(i10)') nstep
  do i=1,11-lenstr(timechar)-1
    timechar(i:i)='0'
  end do

! Make sure that the new run doesn't overwrite the file from the old run 

    filetype = '.2D'
    
    if(save2Dsep) then
       fileaction = 'WRITE'
       filename='./OUT_2D/'//trim(case)//'_'//trim(caseid)//'_'// &
          trim(date_pr)//'_'//timechar(1:10)//filetype 
    else
       fileaction = 'READWRITE'
       filename='./OUT_2D/'//trim(case)//'_'//trim(caseid)//filetype 
    end if
    open(46,file=filename,status='unknown', &
                  form='unformatted',BUFFEREDYES ACTION=trim(fileaction))
    print*, 'Writting to file: '//trim(filename)
    if(.not.save2Dsep) then 
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
 222   continue
       print*,'nsteplast=',nsteplast
    end if

    write(46) nstep, dolatlon, time, datechar, timeUTsec
    write(46) 'atm'
    write(46) nx,ny,nzm,nsubdomains, nsubdomains_x,nsubdomains_y,nfields
    write(46) real(dx,4),lat_gl,latv_gl,real(y_gl,4),real(yv_gl,4),wgty
    write(46) real(dy,4),lon_gl,lonu_gl
    write(46) dble(nstep)*dt/(3600._8*24._8)+day0

 end if! masterproc

else

! open netcdf file:

 write(timechar,'(i10)') nstep
 do k=1,11-lenstr(timechar)-1
   timechar(k:k)='0'
 end do

 filename ='./OUT_2D/'//trim(case)//'_'//trim(caseid)//'_'// &
       trim(date_pr)//'_'//timechar(1:10)//'.2D_atm.nc'

 call open_file2D_pnetcdf(filename)

 if(masterproc) print*, 'Writting to file: '//trim(filename)

end if ! .not.dopnetcdf


if(.not.nstep.eq.1.and.save2Davg) then
   coef = 1./float(nsave2D)
else
   coef = 1.
end if
coefa = 1./float(nsave2D)


! 2D fields:


if(doprecip) then

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=prec_xy(i,j)/dtp*86400.
    end do
   end do
  name='Prec'
  long_name='Surface Precip. Rate'
  units='mm/day'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=preca_xy(i,j)/dtpa*86400.
    end do
   end do
  name='Preca'
  long_name='Surface Precip. Rate (aver)'
  units='mm/day'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      precac_xy(i,j) = precac_xy(i,j)+preca_xy(i,j)
      tmp(i,j)=precac_xy(i,j)
    end do
   end do
  name='Precac'
  long_name='Surface Accum Precip.'
  units='mm'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)

  if(snow2Dout) then

    nfields1=nfields1+1
    do j=1,ny
     do i=1,nx
       tmp(i,j)=precs_xy(i,j)/dtp*86400.
     end do
    end do
   name='Precs'
   long_name='Surface LWE Snow Rate'
   units='mm/day'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                                save2Dbin,dompi,rank,nsubdomains)
    nfields1=nfields1+1
    do j=1,ny
     do i=1,nx
       tmp(i,j)=precsa_xy(i,j)/dtpa*86400.
     end do
    end do
   name='Precsa'
   long_name='Surface LWE Snow Rate (aver)'
   units='mm/day'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                                save2Dbin,dompi,rank,nsubdomains)
 
    nfields1=nfields1+1
    do j=1,ny
     do i=1,nx
       precsac_xy(i,j) = precsac_xy(i,j)+precsa_xy(i,j)
       tmp(i,j)=precsac_xy(i,j)
     end do
    end do
   name='Precsac'
   long_name='Surface Accum iLWE Snow'
   units='mm'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                                .true.,dompi,rank,nsubdomains)

  end if ! snow2Dout
end if

if(doterrain) then
  nfields1=nfields1+1
  do j=1,ny
    do i=1,nx
      tmp(i,j)=elevation(i,j)
!      tmp(i,j)=elevationg(i,j)
    end do
  end do
  name='ZSFC'
  long_name='Surface Elevation (on grid)'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
  if(dobuildings) then

    nfields1=nfields1+1
    do j=1,ny
      do i=1,nx
        tmp(i,j)=bld_mask(i,j)
    end do
    end do
    name='BLDMASK'
    long_name='Building Mask'
    units=' '
    call compress2D(tmp,nx,ny,name,long_name,units, &
                             save2Dbin,dompi,rank,nsubdomains)

    nfields1=nfields1+1
    do j=1,ny
      do i=1,nx
        tmp(i,j)=bld_index(i,j)
    end do
    end do
    name='BLDNUMBER'
    long_name='Building Number Mask'
    units=' '
    call compress2D(tmp,nx,ny,name,long_name,units, &
                             save2Dbin,dompi,rank,nsubdomains)

    if(.not.(doAC.and.doHEAT)) then
      nfields1=nfields1+1
      do j=1,ny
        do i=1,nx
          if(bld_index(i,j).gt.0) then
            tmp(i,j)=t_bld(bld_index(i,j))
          else
            tmp(i,j) = 0.
          end if
        end do
      end do
      name='TBLD'
      long_name='Building Internal Temperature'
      units='C'
      call compress2D(tmp,nx,ny,name,long_name,units, &
                             save2Dbin,dompi,rank,nsubdomains)
    end if

    if(doAC.or.doHEAT) then
      nfields1=nfields1+1
      do j=1,ny
        do i=1,nx
          if(bld_index(i,j).gt.0) then
            tmp(i,j) = (building_cooling_energy(bld_index(i,j))+ &
                     building_heating_energy(bld_index(i,j)))*0.000001
          else
            tmp(i,j) = 0.
          end if
        end do
      end do
      name='HVAC Energy'
      long_name='Building HVAC Energy'
      units='MJ'
      call compress2D(tmp,nx,ny,name,long_name,units, &
                             save2Dbin,dompi,rank,nsubdomains)
    end if

    nfields1=nfields1+1
    do j=1,ny
      do i=1,nx
        tmp(i,j)=facemask(i,j)
      end do
    end do
    name='FACEMASK'
    long_name='Building Face Mask'
    units=' '
    call compress2D(tmp,nx,ny,name,long_name,units, &
                             save2Dbin,dompi,rank,nsubdomains)

    if(doshadows) then
      nfields1=nfields1+1
      do j=1,ny
        do i=1,nx
          tmp(i,j)=shadow_height(i,j)
        end do
      end do
      name='SHADOW_Z'
      long_name='Shadow Height'
      units='m'
      call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)

      nfields1=nfields1+1
      do j=1,ny
        do i=1,nx
          tmp(i,j)=shadow_mask(i,j)
        end do
      end do
      name='SHADOW'
      long_name='Shadow Mask'
      units=' '
      call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)

      nfields1=nfields1+1
      do j=1,ny
        do i=1,nx
          tmp(i,j)=sum(cell_sunny(i,j,:))
        end do
      end do
      name='FACESUN'
      long_name='Number of sunny faces'
      units=' '
      call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)

    end if

    if(doirgroundfromwalls) then

      nfields1=nfields1+1
       do j=1,ny
         do i=1,nx
           tmp(i,j)=ind_wall_sfc(i,j)
         end do
       end do
       name='INDVIS'
       long_name='nimber of visible faces on the ground'
       units=' '
       call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)

      nfields1=nfields1+1
       do j=1,ny
         do i=1,nx
           tmp(i,j)=flx_wall_sfc(i,j)
         end do
       end do
       name='FLXSB'
       long_name='IR flux on surface from walls'
       units='W/m2'
       call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)

      nfields1=nfields1+1
       do j=1,ny
         do i=1,nx
           tmp(i,j)=sangle_wall_sfc(i,j)
         end do
       end do
       name='SANGLE'
       long_name='Solid angle for IR flux from walls'
       units='srad'
       call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)

    end if
  end if
end if

if(.not.SFC_FLX_FXD) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      tmp(i,j)=shf_xy(i,j)*rhow(k)*cp*coef
    end do
   end do
  name='SHF'
  long_name='Sensible Heat Flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      k = k_terra(i,j)
      tmp(i,j)=lhf_xy(i,j)*rhow(k)*lcond*coef
    end do
   end do
  name='LHF'
  long_name='Latent Heat Flux'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
  nfields1=nfields1+1
  do j=1,ny
   do i=1,nx
     k = k_terra(i,j)
     tmp(i,j)=taux_xy(i,j)*rhow(k)*coef
   end do
  end do
  name='TAUX'
  long_name='Surface Stress in x'
  units='N/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
  nfields1=nfields1+1
  do j=1,ny
   do i=1,nx
     k = k_terra(i,j)
     tmp(i,j)=tauy_xy(i,j)*rhow(k)*coef
   end do
  end do
  name='TAUY'
  long_name='Surface Stress in y'
  units='N/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
  if(dodns) then
     nfields1=nfields1+1
     do j=1,ny
      do i=1,nx
        tmp(i,j)=shf_top_xy(i,j)*rhow(nz)*cp*coef
      end do
     end do
    name='SHF_TOP'
    long_name='Sensible Heat Flux at top'
    units='W/m2'
    call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)
     nfields1=nfields1+1
     do j=1,ny
      do i=1,nx
        tmp(i,j)=lhf_top_xy(i,j)*rhow(nz)*lcond*coef
      end do
     end do
    name='LHF_TOP'
    long_name='Latent Heat Flux at top'
    units='W/m2'
    call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)
    nfields1=nfields1+1
    do j=1,ny
     do i=1,nx
       tmp(i,j)=taux_top_xy(i,j)*rhow(nz)*coef
     end do
    end do
   name='TAUX_TOP'
   long_name='Surface Stress in x at top'
   units='N/m2'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                                save2Dbin,dompi,rank,nsubdomains)
    nfields1=nfields1+1
    do j=1,ny
     do i=1,nx
       tmp(i,j)=tauy_top_xy(i,j)*rhow(nz)*coef
     end do
    end do
   name='TAUY_TOP'
   long_name='Surface Stress in y at top'
   units='N/m2'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                                save2Dbin,dompi,rank,nsubdomains)

  end if

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=t2m_xy(i,j)*coef
    end do
   end do
  name='T2m'
  long_name='2-m temperature'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

!   nfields1=nfields1+1
!   do j=1,ny
!    do i=1,nx
!      k = k_terra(i,j)
!      tmp(i,j)=heat_index(t2m_xy(i,j)*coef,q2m_xy(i,j)*coef,pres(k))
!    end do
!   end do
!  name='HI2m'
!  long_name='2-m heat index'
!  units='F'
!  call compress2D(tmp,nx,ny,name,long_name,units, &
!                               save2Dbin,dompi,rank,nsubdomains)
!
!   nfields1=nfields1+1
!   do j=1,ny
!    do i=1,nx
!      k = k_terra(i,j)
!      tmp(i,j)=q2m_xy(i,j)*coef/qsatw(t2m_xy(i,j)*coef,pres(k))*100.
!    end do
!   end do
!  name='RH2m'
!  long_name='2-m relative himidity'
!  units='%'
!  call compress2D(tmp,nx,ny,name,long_name,units, &
!                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=q2m_xy(i,j)*coef*1.e3
    end do
   end do
  name='Q2m'
  long_name='2-m humidity'
  units='g/kg'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)


   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=u10m_xy(i,j)*coef
    end do
   end do
  name='U10m'
  long_name='10-m zonal wind'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=v10m_xy(i,j)*coef
    end do
   end do
  name='V10m'
  long_name='10-m meridional wind'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=u10ma_xy(i,j)*coefa
    end do
   end do
  name='U10ma'
  long_name='Mean 10-m zonal wind'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=v10ma_xy(i,j)*coefa
    end do
   end do
  name='V10ma'
  long_name='Mean 10-m meridional wind'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)


   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=sqrt(gust10m_xy(i,j))
    end do
   end do
  name='GUST10m'
  long_name='10-m gust'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=ra_xy(i,j)*coef
    end do
   end do
  name='Ra'
  long_name='surface flux resistance'
  units='s/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)


   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=landtype(i,j)
    end do
   end do
  name='LAND'
  long_name='Land Type'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(dolongwave) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwns_xy(i,j)*coef
     end do
   end do
  name='LWNS'
  long_name='Net LW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwnsc_xy(i,j)*coef
     end do
   end do
  name='LWNSC'
  long_name='Net clear-sky LW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwds_xy(i,j)*coef
     end do
   end do
  name='LWDS'
  long_name='Downward LW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwdsc_xy(i,j)*coef
     end do
   end do
  name='LWDSC'
  long_name='Downward clear-sky LW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwnt_xy(i,j)*coef
     end do
   end do
  name='LWNT'
  long_name='Net LW at TOA'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwntc_xy(i,j)*coef
     end do
   end do
  name='LWNTC'
  long_name='Clear-Sky Net LW at TOA'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

  if(save2Dradac) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwnsac_xy(i,j)
     end do
   end do
   name='LWNSAC'
   long_name='Net LW at the surface (accum)'
   units='J/m2'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwdsac_xy(i,j)
     end do
   end do
   name='LWDSAC'
   long_name='Downward LW at the surface (accum)'
   units='J/m2'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwntac_xy(i,j)
     end do
   end do
   name='LWNTAC'
   long_name='Net LW at TOA (accum)'
   units='J/m2'
   call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)

  end if

  if(.not.save2Davg.and.save2Drada) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwnsa_xy(i,j)*coefa
     end do
   end do
  name='LWNSA'
  long_name='Net LW at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwdsa_xy(i,j)*coefa
     end do
   end do
  name='LWDSA'
  long_name='Downward LW at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwnta_xy(i,j)*coefa
     end do
   end do
  name='LWNTA'
  long_name='Net LW at TOA (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwnsca_xy(i,j)*coefa
     end do
   end do
  name='LWNSCA'
  long_name='Net clear-sky LW at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwdsca_xy(i,j)*coefa
     end do
   end do
  name='LWDSCA'
  long_name='Downward clear-sky LW at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=lwntca_xy(i,j)*coefa
     end do
   end do
  name='LWNTCA'
  long_name='Net Clear-sky LW at TOA (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
  end if
end if

if(doshortwave) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=solin_xy(i,j)*coef
     end do
   end do
  name='SOLIN'
  long_name='Solar TOA insolation'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swns_xy(i,j)*coef
     end do
   end do
  name='SWNS'
  long_name='Net SW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swds_xy(i,j)*coef
     end do
   end do
  name='SWDS'
  long_name='Downward SW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swnsc_xy(i,j)*coef
     end do
   end do
  name='SWNSC'
  long_name='Net Clear-sky SW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swnt_xy(i,j)*coef
     end do
   end do
  name='SWNT'
  long_name='Net SW at TOA'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swntc_xy(i,j)*coef
     end do
   end do
  name='SWNTC'
  long_name='Net Clear-Sky SW at TOA'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swdsc_xy(i,j)*coef
     end do
   end do
  name='SWDSC'
  long_name='Downward clear-sky SW at the surface'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
  if(save2Dradac) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swnsac_xy(i,j)
     end do
   end do
  name='SWNSAC'
  long_name='Net SW at the surface (accum)'
  units='J/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swdsac_xy(i,j)
     end do
   end do
  name='SWDSAC'
  long_name='Downward SW at the surface (accum)'
  units='J/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swntac_xy(i,j)
     end do
   end do
  name='SWNTAC'
  long_name='Net SW at TOA (accum)'
  units='J/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
  end if

  if(.not.save2Davg.and.save2Drada) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swnsa_xy(i,j)*coefa
     end do
   end do
  name='SWNSA'
  long_name='Net SW at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swdsa_xy(i,j)*coefa
     end do
   end do
  name='SWDSA'
  long_name='Downward SW at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swnta_xy(i,j)*coefa
     end do
   end do
  name='SWNTA'
  long_name='Net SW at TOA (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swnsca_xy(i,j)*coefa
     end do
   end do
  name='SWNSCA'
  long_name='Net SW clear-sky at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swdsca_xy(i,j)*coefa
     end do
   end do
  name='SWDSCA'
  long_name='Downward clear-sky SW at the surface (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swntca_xy(i,j)*coefa
     end do
   end do
  name='SWNTCA'
  long_name='Net clear-sky SW at TOA (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=solina_xy(i,j)*coefa
     end do
   end do
  name='SOLINA'
  long_name='SOLAR TOA Insolation (aver)'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
  end if
end if

if(docloud) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=cw_xy(i,j)*coef
     end do
   end do
  name='CWP'
  long_name='Cloud Water Path'
  units='kg/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=iw_xy(i,j)*coef
     end do
   end do
  name='IWP'
  long_name='Ice Path'
  units='kg/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=cld_xy(i,j)*coef
     end do
   end do
  name='CLD'
  long_name='Cloud Cover'
  units=' '
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=eis_xy(i,j)*coef
     end do
   end do
  name='EIS'
  long_name='Estimated Inversion Strenth'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

  ! compute CAPE and CIN for deep domains only. -MK 
  if(.not.ncycle.gt.ncycle_max.and.pres(nzm).lt.100.) then
    call cape_cin(cape,cin)
  else
    cape = 0.
    cin = 0.
  end if 

  nfields1=nfields1+1
  do j=1,ny
    do i=1,nx
       tmp(i,j)=cape(i,j)
    end do
  end do
  name='CAPE'
  long_name='Convectively Avalable Potential Energy'
  units='J/kg'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do j=1,ny
    do i=1,nx
       tmp(i,j)=cin(i,j)
    end do
  end do
  name='CIN'
  long_name='Convective Inhibition'
  units='J/kg'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(doprecip) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=rw_xy(i,j)*coef
     end do
   end do
  name='RWP'
  long_name='Rain Water Path'
  units='kg/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=sw_xy(i,j)*coef
     end do
   end do
  name='SWP'
  long_name='Snow Water Path'
  units='kg/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=gw_xy(i,j)*coef
     end do
   end do
  name='GWP'
  long_name='Graupel Water Path'
  units='kg/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=ZdBZ_xy(i,j)
     end do
   end do
  name='ZdBZ'
  long_name='Composite Radar Reflectivity (inst)'
  units='dBZ'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)


end if

if(docloud) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=pws_xy(i,j)*coef
     end do
   end do
  name='PWS'
  long_name='Saturated Precipitable Water'
  units='kg/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=fse_xy(i,j)*coef
     end do
   end do
  name='FSE'
  long_name='Frozen Mosist Static Energy'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=pw_xy(i,j)*coef
     end do
   end do
  name='PW'
  long_name='Precipitable Water'
  units='mm/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=rh200_xy(i,j)*coef
     end do
   end do
  name='RH200'
  long_name='Relative Humidity 200mb'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=rh500_xy(i,j)*coef
     end do
   end do
  name='RH500'
  long_name='Relative Humidity 500mb'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=rh700_xy(i,j)*coef
     end do
   end do
  name='RH700'
  long_name='Relative Humidity 700mb'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=rh850_xy(i,j)*coef
     end do
   end do
  name='RH850'
  long_name='Relative Humidity 850mb'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

if(donudge3D.or.doregion) then

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=pwobs_xy(i,j)*coef
     end do
   end do
  name='PWOBS'
  long_name='Observed Precipitable Water'
  units='mm/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=uobs200_xy(i,j)*coef
     end do
   end do
  name='UOBS200'
  long_name='Observed U at 200 mb'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=vobs200_xy(i,j)*coef
     end do
   end do
  name='VOBS200'
  long_name='Observed V at 200 mb'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=uobs850_xy(i,j)*coef
     end do
   end do
  name='UOBS850'
  long_name='Observed U at 850 mb'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=vobs850_xy(i,j)*coef
     end do
   end do
  name='VOBS850'
  long_name='Observed V at 850 mb'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=tau_nudge3Du(i,j)
     end do
   end do
  name='TAU_3DU'
  long_name='Nudging coef'
  units='1/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains) 

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=tau_nudge3Dv(i,j)
     end do
   end do
  name='TAU_3DV'
  long_name='Nudging coef'
  units='1/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=tau_nudge3D(i,j)
     end do
   end do
  name='TAU_3D'
  long_name='Nudging coef'
  units='1/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=ta_xy(i,j)*coef
     end do
   end do
  name='TA'
  long_name='Mass-weighted Column Temperature'
  units='mm'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=usfc_xy(i,j)*coef + ug
     end do
   end do
  name='USFC'
  long_name='U at the surface'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=tsfc_xy(i,j)*coef 
     end do
   end do
  name='TSFC'
  long_name='Temperature near the surface'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=qsfc_xy(i,j)*coef *1.e3
     end do
   end do
  name='QSFC'
  long_name='Vapor near the surface'
  units='g/kg'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=u200_xy(i,j)*coef + ug
     end do
   end do
  name='U200'
  long_name='U at 200 mb'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
      tmp(i,j)=vsfc_xy(i,j)*coef + vg
     end do
   end do
  name='VSFC'
  long_name='V at the surface'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=sqrt(gustsfc_xy(i,j))
    end do
   end do
  name='GUSTSFC'
  long_name='Gust at the surface'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=v200_xy(i,j)*coef + vg
     end do
   end do
  name='V200'
  long_name='V at 200 mb'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)


  nfields1=nfields1+1
  tmp = 0.
  do j=1,ny
    do i=1,nx
      tmp(i,j)=vor200_xy(i,j)*coef
    end do
  end do
  name='VOR200'
  long_name=' Vorticity at 200 mb'
  units='1/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)


  nfields1=nfields1+1
  tmp = 0.
  do j=1,ny
    do i=1,nx
      tmp(i,j)=vor850_xy(i,j)*coef
    end do
  end do
  name='VOR850'
  long_name=' Vorticity at 850 mb'
  units='1/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                                 save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=w500_xy(i,j)*coef
     end do
   end do
  name='W500'
  long_name='W at 500 mb'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   factor = 1./ggr
   do j=1,ny
     do i=1,nx
       tmp(i,j)=phi500_xy(i,j)*coef*factor + z(k500)
     end do
   end do
  name='PHI500'
  long_name='Geopotential 500 mb'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=omega200_xy(i,j)*coef
     end do
   end do
  name='OM200'
  long_name='Pressure velocity at 200 mb'
  units='Pa/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=omega500_xy(i,j)*coef
     end do
   end do
  name='OM500'
  long_name='Pressure velocity at 500 mb'
  units='Pa/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=omega700_xy(i,j)*coef
     end do
   end do
  name='OM700'
  long_name='Pressure velocity at 700 mb'
  units='Pa/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=omega850_xy(i,j)*coef
     end do
   end do
  name='OM850'
  long_name='Pressure velocity at 850 mb'
  units='Pa/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)


if(doslabocean) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=qocean_xy(i,j)*coef
    end do
   end do
  name='QOCN'
  long_name='Deep Ocean Cooling'
  units='W/m2'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if((ocean_type.ne.0.or.doslabocean.or.doequilocean.or.dodynamicocean).and..not.dossthomo.or. &
    (ISLAND.or.LAND).and.SLM.or.readsst) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      tmp(i,j)=sst_xy(i,j)*coef+t00
    end do
   end do
  name='SKT'
  long_name='Skin Surface Temperature'
  units='K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

if(ISLAND.or.doseaice) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
      if(doseaice) then
        tmp(i,j)=landmask(i,j)+landicemask(i,j)+3*icemask(i,j)
      else
        tmp(i,j)=landmask(i,j)+landicemask(i,j)+3*seaicemask(i,j)
      end if
    end do
   end do
  name='LANDMASK'
  long_name='Landmask:0-ocean,1-land,2-landice,3-seaice'
  units=' '
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

end if
if(doseaice.and.doseaiceevol) then
   nfields1=nfields1+1
   do j=1,ny
    do i=1,nx
        tmp(i,j)=seaice_h(i,j)
    end do
   end do
  name='SEAICEH'
  long_name='Sea-ice thinkness'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
end if

  if(dotracers) then
   do itr = 1,ntracers
    nfields1=nfields1+1
    do j=1,ny
      do i=1,nx
        tmp(i,j)=tr_xy(i,j,itr)*coef
      end do
    end do
    name=trim(tracername(itr))
    long_name=trim(tracername(itr))
    units=trim(tracerunits(itr))
    call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
   end do
   if(dotrsfcflux) then
     do itr = 1,ntracers
      nfields1=nfields1+1
      do j=1,ny
        do i=1,nx
          tmp(i,j)=tr_ac_xy(i,j,itr)
        end do
      end do
      name=trim(tracername(itr))//'_ACS'
      long_name='Accumulated on sfc '//trim(tracername(itr))
      units=trim(tracerunits(itr))//' kg/m3'
      call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
     end do
   end if
   do itr = 1,ntracers
    nfields1=nfields1+1
    do j=1,ny
      do i=1,nx
        tmp(i,j)=tr_acs_xy(i,j,itr)
      end do
    end do
    name=trim(tracername(itr))//'_AC'
    long_name='Accumulated near sfc '//trim(tracername(itr))
    units=trim(tracerunits(itr))//' kg/m3*s'
    call compress2D(tmp,nx,ny,name,long_name,units, &
                               .true.,dompi,rank,nsubdomains)
   end do

  end if

!=====================================================
! UW ADDITIONS
 
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=psfc_xy(i,j)*coef
     end do
   end do
  name='PSFC'
  long_name='Surface Pressure'
  units='hPa'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=swvp_xy(i,j)*coef
     end do
   end do
  name='SWVP'
  long_name='Saturated Water Vapor Path'
  units='mm'
  call compress2D(tmp,nx,ny,name,long_name,units, &
       save2Dbin,dompi,rank,nsubdomains)

   ! 850 mbar zonal velocity
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=u850_xy(i,j)*coef + ug
     end do
   end do
  name='U850'
  long_name='850 mbar zonal velocity'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)

   ! meridional wind at 850 mbar
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=v850_xy(i,j)*coef + vg
     end do
   end do
  name='V850'
  long_name='850 mbar meridional velocity'
  units='m/s'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2Dbin,dompi,rank,nsubdomains)
if(docloud) then

  ! cloud top height
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j)=cloudtopheight(i,j)/1000.
      end do
   end do
   name='ZC'
   long_name='Cloud top height (Instantaneous)'
   units='km'
   call compress2D(tmp,nx,ny,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   ! cloud top temperature
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j)=cloudtoptemp(i,j)
      end do
   end do
   name='TB'
   long_name='Cloud top temperature (Instantaneous)'
   units='K'
   call compress2D(tmp,nx,ny,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   ! echo top height
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j)=echotopheight(i,j)/1000.
      end do
   end do
   name='ZE'
   long_name='Echo top height (Instantaneous)'
   units='km'
   call compress2D(tmp,nx,ny,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)

   if(save2Davg) then
   ! cloud cover
   nfields1=nfields1+1
   do j=1,ny
      do i=1,nx
         tmp(i,j)=cloudcover(i,j)
      end do
   end do
   name='CLDC'
   long_name='Cloud cover (Instantaneous)'
   units=''
   call compress2D(tmp,nx,ny,name,long_name,units, &
        save2Dbin,dompi,rank,nsubdomains)
   end if

end if

! END UW ADDITIONS
!=====================================================

call task_barrier()

!===================================================================


if(nfields.ne.nfields1) then
  if(masterproc) print*,'write_fields2D: error in nfields!!',' nfields=',nfields,'nfields1=',nfields1
  call task_abort()
end if

if(dopnetcdf) then
     call close_file_pnetcdf()
     dopnetcdf = .false.
else
  if(masterproc) then
     close(46)
     if(save2Dsep.and.dogzip2D) call systemf('gzip -f '//filename)
  endif
end if
if(masterproc) print*, 'Done. ', nfields1,'fields'


end
