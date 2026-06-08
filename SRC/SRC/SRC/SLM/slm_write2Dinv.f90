#include "./../fppmacros"

!=====================================================================
! Based on write_fields2D.f90 of SAM
! Write 2D invariant variables (don;t change) into the output file in OUT_INV
!	
!=====================================================================
subroutine slm_write2Dinv

use grid
use slm_vars
use params, only: dolatlon, docheck
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

if(docheck) return

dopnetcdf = save2DLnetcdf

nfields = 17+(nsoil*13)

nfields1 = 0 

if(.not.dopnetcdf) then

 if(masterproc) then

! Make sure that the new run doesn't overwrite the file from the old run

    if(.not.save2DLsep.and.notopened2DL.and.(nrestart.eq.0.or.nrestart.eq.2)) then
      filestatus='new'
    else
      filestatus='unknown'
    end if

    filetype = '.2D'

    filename='./OUT_INV/'//trim(case)//'_'//trim(caseid)//'_INV'//filetype
          open(46,file=filename,status='unknown',form='unformatted',BUFFEREDYES ACTION='WRITE')
    print*, 'Writting to file: '//trim(filename)
    write(46) 0, dolatlon, time, datechar, timeUTsec
    write(46) 'lnd'
    write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
    write(46) real(dx,4),lat_gl,latv_gl,real(y_gl,4),real(yv_gl,4)
    write(46) real(dy,4),lon_gl,lonu_gl
    write(46) 0._8

 end if! masterproc

else

! open netcdf file:

 write(timechar,'(i10)') nstep
 do k=1,11-lenstr(timechar)-1
   timechar(k:k)='0'
 end do

 filename ='./OUT_INV/'//trim(case)//'_'//trim(caseid)//'_INV.2D_lnd.nc'

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
!! Albedos
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(albedovis_v(i,j),4)
     end do
   end do
  name='ALBVV'
  long_name='Albedo Vegetation Visible'
  units=' '
 call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(albedonir_v(i,j),4)
     end do
   end do
  name='ALBVN'
  long_name='Albedo Vegetation Near-infrared'
  units=' '
 call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(albedovis_s(i,j),4)
     end do
   end do
  name='ALBSV'
  long_name='Albedo Soil Visible'
  units=' '
 call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(albedonir_s(i,j),4)
     end do
   end do
  name='ALBSN'
  long_name='Albedo SOIL Near-Infrared'
  units=' '
 call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!====================================================================
!! Vegetation Height
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(ztop(i,j),4)
     end do
   end do
  name='ZTOP'
  long_name='Vegetation Height'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
!! Vegetation roughness
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(z0_sfc(i,j),4)
     end do
   end do
  name='Z0_SFC'
  long_name='Vegetation Roughness'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!====================================================================
! Displacement height 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(disp_hgt(i,j),4)
     end do
   end do
  name='DISPHGT'
  long_name='Displacement Height'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! Parameter Khai_L
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(Khai_L(i,j),4)
     end do
   end do
  name='KHAI_L'
  long_name='Parameter Khai_L'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! iRoot Depth 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(rootL(i,j),4)
     end do
   end do
  name='ROOTL'
  long_name='Root Depth'
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! parameter root_a 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(root_a(i,j),4)
     end do
   end do
  name='ROOT_A'
  long_name='parameter root_a'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! parameter root_b   
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(root_b(i,j),4)
     end do
   end do
  name='ROOT_B'
  long_name='parameter root_b'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!====================================================================
! Minimum r_c 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(Rc_min(i,j),4)
     end do
   end do
  name='RC_MIN'
  long_name='Minimum Canopy Resistance'
  units='s/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
!!Parameter for stomatal Resitance Rgl
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(Rgl(i,j),4)
     end do
   end do
  name='RGL'
  long_name='Parameter Rgl'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! Parameter for stomatal Resistance hs_rc
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(hs_rc(i,j),4)
     end do
   end do
  name='HS_RC'
  long_name='Parameter hs_rc'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! Basal area index 
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(BAI(i,j),4)
     end do
   end do
  name='BAI'
  long_name='Basal area index'
  units='sqf/acre'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! imperviousness
!====================================================================
    nfields1=nfields1+1
    do j=1,ny
      do i=1,nx
        tmp(i,j)=real(IMPERV(i,j),4)
      end do
    end do
   name='IMPERV'
   long_name='Impreviousness'
   units=''
   call compress2D(tmp,nx,ny,name,long_name,units, &
                                save2DLbin,dompi,rank,nsubdomains)
!====================================================================
! Vegetation Mask
!====================================================================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       if(vegetated(i,j)) then
         tmp(i,j)=1.
       else
         tmp(i,j)=0.
       end if  
     end do
   end do
  name='VEGMASK'
  long_name='Vegetation Mask 1 - veg, 0 - not'
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)


do k = 1,nsoil
write(numchar,'(i0)') k
!========================
! Clay Content
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(CLAY(k,i,j),4)
     end do
   end do
  name='CLAY'//numchar
  long_name='Clay Content layer'//numchar
  units='%'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! Sand Content
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(SAND(k,i,j),4)
     end do
   end do
  name='SAND'//numchar
  long_name='Sand Content layer'//numchar
  units='%'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! soil wetness 
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(sst_cond(k,i,j),4)
     end do
   end do
  name='COND'//numchar
  long_name='Dry Soil thermal Conductivity layer'//numchar
  units='W/mK'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil Porocity
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(poro_soil(k,i,j),4)
     end do
   end do
  name='PORO'//numchar
  long_name='Soil Porocity layer'//numchar
  units='m/m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil Bconst
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(Bconst(k,i,j),4)
     end do
   end do
  name='BCONST'//numchar
  long_name='Soil Bconst layer'//numchar
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! soil moisture potential at saturation, mm
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(m_pot_sat(k,i,j),4)
     end do
   end do
  name='POTSAT'//numchar
  long_name='moisture potential at sat layer'//numchar
  units='mm'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil soil heat capacity
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(sst_capa(k,i,j),4)
     end do
   end do
  name='CAPA'//numchar
  long_name='soil heat capacity layer'//numchar
  units='J/m3/K'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! volumetric moisture content at field capacity
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(theta_FC(k,i,j),4)
     end do
   end do
  name='THETA_FC'//numchar
  long_name='volum. moisture content at field capacity layer'//numchar
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! volumetric moisture content at wilting point
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(theta_WP(k,i,j),4)
     end do                    
   end do
  name='THETA_WP'//numchar
  long_name='volum. moisture content at wilting point layer'//numchar
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! soil wetness at field capacity
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(w_s_FC(k,i,j),4)
     end do
   end do
  name='W_FC'//numchar
  long_name='soil wetness at field capacity layer'//numchar
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil wetness at wilting point
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(w_s_WP(k,i,j),4)
     end do
   end do
  name='W_WP'//numchar
  long_name='soil wetness at wilting point layer'//numchar
  units=''
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

!========================
! soil layer-center depth
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(node_z(k,i,j),4)
     end do
   end do
  name='S_Z'//numchar
  long_name='Soil mid-layer depth'//numchar
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)
!========================
! soil layer thickness
!========================
   nfields1=nfields1+1
   do j=1,ny
     do i=1,nx
       tmp(i,j)=real(s_depth_mm(k,i,j),4)*0.001
     end do
   end do
  name='S_DZ'//numchar
  long_name='Soil layer thickness'//numchar
  units='m'
  call compress2D(tmp,nx,ny,name,long_name,units, &
                               save2DLbin,dompi,rank,nsubdomains)

end do ! k

call task_barrier()

if(nfields.ne.nfields1) then
  if(masterproc) print*,'land_fields2D: error in nfields!!',' nfields=',nfields,'nfields1=',nfields1
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


end subroutine slm_write2Dinv
