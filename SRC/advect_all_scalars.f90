subroutine advect_all_scalars()

  use vars
  use microphysics
  use chemistry, only: gchem_field, ngchem_fields, aqchem_field, aqchem_gasprod_field, archem_field, gchwle, gchadv, &
       aqchwle, aqchadv, archwle, archadv
  use chemistry_params, only : do_iepox_droplet_chem, do_iepox_aero_chem
  use chem_aqueous, only: naqchem_fields
  use chem_aerosol, only: narchem_fields
  use sgs
  use tracers
  use params, only: dotracers, docup, dochem
  use stat_coars
  use terrain, only: terraw
  implicit none
  real dummy(nz)
  real dtdx, dtdy, dtdz, rhox, rhoy, rhoz , a1, a2
  real, allocatable :: f(:,:,:)
  integer i,j,k
  logical compute_variance_stats

call t_startf ('advect_scalars')

!---------------------------------------------------------
!      velocities at n+1/2 time step :

dtdx = dtn/dx
dtdy = dtn/dy
dtdz = dtn/dz
a1 = 0.5
a2 = 0.5
if(nstep.eq.1.and.ncycle.eq.1) then
 a1 = 1.
 a2 = 0.
end if

do k=1,nzm
  do j=dimy1_u,dimy2_u
   rhox = rho(k)*dtdx*adz(k)*ady(j)
   do i=dimx1_u,dimx2_u
     u1(i,j,k) = a1*u(i,j,k)*rhox+a2*u1(i,j,k)*dtn
   end do
  end do
end do
do k=1,nzm
  do j=dimy1_v,dimy2_v
   rhoy = rho(k)*dtdy*muv(j)*adz(k)
   do i=dimx1_v,dimx2_v
     v1(i,j,k) = a1*v(i,j,k)*rhoy+a2*v1(i,j,k)*dtn
   end do
  end do
end do
do k=1,nzm
  do j=dimy1_w,dimy2_w
   rhoz = rhow(k)*dtdz*ady(j)*mu(j)
   do i=dimx1_w,dimx2_w
     w1(i,j,k) = a1*w(i,j,k)*rhoz+a2*w1(i,j,k)*dtn
   end do
  end do
end do


!---------------------------------------------------------
!      advection of scalars :

if(docloud) total_water_adv = total_water_adv - total_water()

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(k,dummy)
     do k = 0,nmicro_fields
       if(k.eq.0) then
        if(collect_coars) then
          allocate(f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm))
          f(:,:,:) = t(:,:,:)
          call advect_scalar(f,dummy,dummy,dummy,dummy,dummy,.true.,.false.,.true.,.false.,250.)
          call coars_fld(misc,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm),fld_flux(:,:,:,4))
          deallocate(f)
        end if
        call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.,.false.,.true.,.false.,250.)
       else
        if(collect_coars) misc = 0.
        if(k.eq.index_water_vapor) then ! transport water-vapor variable no metter what
            call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k), &
                                q2leadv,q2legrad,qwleadv,.true.,.true.,.false.,.false.,0.)
        else if(docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
         .or. doprecip.and.flag_precip(k).eq.1 ) then
             if (flag_advect(k).eq.1) then
                call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k), &
                                        dummy,dummy,dummy,.false.,.true.,.false.,.false.,0.)
             end if
        end if 
        if(collect_coars) call coars_fld(misc,mu(1:ny),ady(1:ny), &
                                                 terraw(1:nx,1:ny,1:nzm),fld_flux(:,:,:,4+k))
       end if 
     end do
!$OMP END PARALLEL DO 

if(docloud) total_water_adv = total_water_adv + total_water()

!
!    Advection of sgs prognostics:
!

     if(dosgs.and.advect_sgs) then
       do k = 1,nsgs_fields
           call advect_scalar(sgs_field(:,:,:,k),sgsadv(:,k),sgswle(:,k), &
                                         dummy,dummy,dummy,.false.,.true.,.false.,.false.,0.)
       end do
     end if


 ! advection of tracers:

     if(dotracers) then
       do k = 1,ntracers
         call advect_scalar(tracer(:,:,:,k),tradv(:,k),trwle(:,k), &
                                         dummy,dummy,dummy,.false.,.true.,.false.,.false.,0.)
       end do
     end if
  
  ! advection of chemistry:
    compute_variance_stats = .false.

    if ( dochem ) then
      do k = 1,ngchem_fields
        call advect_scalar(gchem_field(:,:,:,k),gchadv(:,k),gchwle(:,k),dummy,dummy,dummy, &
              compute_variance_stats, .true., .false., .false., 0.)
          
      end do
    endif

    if ( dochem .and. do_iepox_droplet_chem ) then
      do k = 1,naqchem_fields
        call advect_scalar(aqchem_field(:,:,:,k),aqchadv(:,k),aqchwle(:,k),dummy,dummy,dummy, &
              compute_variance_stats, .true., .false., .false., 0.)
      end do

      do k = 1,naqchem_fields
        call advect_scalar(aqchem_gasprod_field(:,:,:,k),aqchadv(:,k),aqchwle(:,k),dummy,dummy,dummy, &
              compute_variance_stats, .true., .false., .false., 0.)
      end do
    endif

    if (dochem .and. do_iepox_aero_chem ) then 
      do k = 1,narchem_fields
        call advect_scalar(archem_field(:,:,:,k),archadv(:,k),archwle(:,k),dummy,dummy,dummy, &
              compute_variance_stats, .true., .false., .false., 0.)
      end do
    endif

call t_stopf ('advect_scalars')

call t_startf ('precip_fall')
!---------------------------------------------------------
!
!   Precipitation fallout in SAM1MOM, KH3, and DRIZZLE  microphysics
!

    precinst(:,:) = 0.

    if(collect_coars) then
      misc = 0.
      misc1 = 0.
    end if
    if(.not.docup.and.doprecip) then
       total_water_prec = total_water_prec + total_water()
       call micro_precip_fall()
       total_water_prec = total_water_prec - total_water()
       if(collect_coars) then
          call coars_fld(misc,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm), &
                                          fld_flux(:,:,:,4+nmicro_fields+1))
          call coars_fld(misc1,mu(1:ny),ady(1:ny),terraw(1:nx,1:ny,1:nzm), &
                                          fld_flux(:,:,:,4+nmicro_fields+2))
       end if
    end if

!-----------------------------------------------------------
!      sedimentation of cloud-ice for SAM1MOM and KH3 microphysics:

    if(.not.docup.and.docloud.and.doprecip.and..not.dowarmcloud) call ice_fall()
    if(.not.docup.and.docloud.and.docloudfall) call cloud_fall()


call t_stopf ('precip_fall')

end subroutine advect_all_scalars
