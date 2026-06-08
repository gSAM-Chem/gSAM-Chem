subroutine nudging()
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor, nmicro_fields
use sgs, only: sgs_field
use terrain
use stat_coars
use tracers
use check_energy
implicit none

real(8) coef, coef1, xxx, yyy, ddd1, ddd2, iddd, tmp
integer i,j,k,it,jt
real fld, fld1
real tend 
real tau(nzm), taum(nzm), taubound, taubound1, lamda, lamda1, smooth, sss
real nu, tauz(nzm),zzz, tau_max
real(8) buf1(2), buf2(2), veladj, dxy
logical uninitialized 
real tend_in(nx,ny,nzm)
integer, parameter :: kblock = 8
smooth(sss) = 0.5*(1.+cos(3.14159265*sss)) ! smooth transition function 0 < sss < 1

call t_startf ('nudging')

call sumenergy(1,'nudging')

if(donudge3D.or.doregion) then

 if(docolumn) then
  if(masterproc) print*,'ERROR: donudge3D=T or doregion=T, but docolumn is also T! Exitting...'
  call task_abort()
 end if
! ------------------------------------------------------------
! precompute damping coefficients around lateral  boundaries if required

 taubound = 10.*dt     ! nudging time scale at the boundary
 taubound1 = nudge3D_tau ! in the inner domain

 lamda = 1./taubound
 if(taubound1.lt.taubound) then
   lamda1 = 0.
 else
   lamda1 = 1./taubound1
 end if

 if(doregion.and..not.allocated(tau_nudge3D)) then
    if(nudge3D_tau.eq.0..and..not.doregion) then
     if(masterproc) print*,'ERROR: donudge3D=T, but nudge3D_tau is not defined.'
     call task_abort()
    end if
    allocate(tau_nudge3D(nx,ny),tau_nudge3Du(nx+1,ny),tau_nudge3Dv(nx,ny+YES3D))
    tau_nudge3D(:,:) = 0.
    tau_nudge3Du(:,:) = 0.
    tau_nudge3Dv(:,:) = 0.
    if(dobufferzonex) then
      call task_rank_to_index(rank,it,jt)
      ddd1 = bufferzonex/2.*nx_gl*dx*muv_gl(ny_gl/(YES3D+1))
      ddd2 = (nx_gl-1)*dx*muv_gl(ny_gl/(YES3D+1)) - ddd1
      iddd = 1./ddd1
      do i=1,nx
        xxx = (i+it-1.)*dx*muv_gl(ny_gl/(YES3D+1))
        if(xxx.le.ddd1) then
          tmp = lamda1+(lamda-lamda1)*smooth(xxx*iddd)
          tau_nudge3D(i,:) = tmp
          tau_nudge3D(i,:) = tmp
          tau_nudge3Dv(i,:) = tmp
        else if(xxx.ge.ddd2) then
          tmp = lamda1+(lamda-lamda1)*smooth(1.-(xxx-ddd2)*iddd)
          tau_nudge3D(i,:) = tmp
          tau_nudge3D(i,:) = tmp
          tau_nudge3Dv(i,:) = tmp
        end if
      end do
      ddd1 = bufferzonex/2.*(nx_gl+1)*dx*muv_gl(ny_gl/(YES3D+1))
      ddd2 = nx_gl*dx*muv_gl(ny_gl/(YES3D+1)) - ddd1
      iddd = 1./ddd1
      do i=1,nx+1
        xxx = (i+it-1.)*dx*muv_gl(ny_gl/(YES3D+1))
        if(xxx.le.ddd1) then
          tau_nudge3Du(i,:) = lamda1+(lamda-lamda1)*smooth(xxx*iddd)
        else if(xxx.ge.ddd2) then
          tau_nudge3Du(i,:) = lamda1+(lamda-lamda1)*smooth(1.-(xxx-ddd2)*iddd)
        end if
      end do

    end if
    if(RUN3D.and.dobufferzoney) then
      call task_rank_to_index(rank,it,jt)
      ddd1 = bufferzoney/2.*(y_gl(ny_gl)-y_gl(1))
      ddd2 = (1.-bufferzoney/2.)*(y_gl(ny_gl)-y_gl(1))
      iddd = 1./ddd1
      do j=1,ny
        yyy = y_gl(j+jt)-y_gl(1)
        if(yyy.le.ddd1) then
         tmp = lamda1+(lamda-lamda1)*smooth(yyy*iddd)
         do i=1,nx
           tau_nudge3D(i,j) = max(tau_nudge3D(i,j),tmp) 
         end do  
         do i=1,nx+1
           tau_nudge3Du(i,j) = max(tau_nudge3Du(i,j),tmp)
         end do
        else if(yyy.ge.ddd2) then
         tmp = lamda1+(lamda-lamda1)*smooth(1.-(yyy-ddd2)*iddd)
         do i=1,nx
           tau_nudge3D(i,j) = max(tau_nudge3D(i,j),tmp)
         end do
         do i=1,nx+1
           tau_nudge3Du(i,j) = max(tau_nudge3Du(i,j),tmp)
         end do
        end if
      end do
      ddd1 = bufferzoney/2.*(yv_gl(ny_gl+1)-yv_gl(1))
      ddd2 = (1.-bufferzoney/2.)*(yv_gl(ny_gl+1)-yv_gl(1))
      iddd = 1./ddd1
      do j=1,ny+YES3D
        yyy = yv_gl(j+jt)-yv_gl(1)
        if(yyy.le.ddd1) then
          tmp = lamda1+(lamda-lamda1)*smooth(yyy*iddd)
          do i=1,nx
           tau_nudge3Dv(i,j) = max(tau_nudge3Dv(i,j),tmp)
          end do
        else if(yyy.ge.ddd2) then
          tmp = lamda1+(lamda-lamda1)*smooth(1.-(yyy-ddd2)*iddd)
          do i=1,nx
           tau_nudge3Dv(i,j) = max(tau_nudge3Dv(i,j),tmp)
          end do
        end if
      end do
    endif
  end if

!--------------------------------------------------------
! prepare damping to obs at the domain top:

  tauz = 0.
  if(doregion) then
    tau_max = 1./taubound
    do k=1,nzm
     nu = (zi(k)-zi(1))/(zi(nzm)-zi(1))
     if(nu.gt.nub) then
        tauz(k) = tau_max*sin(0.5*pi*(nu-nub)/(1.-nub))**2
    !    zzz = 100.*((nu-nub)/(1.-nub))**2
    !    tauz(k) = tau_max*zzz/(1.+zzz)
     end if
    end do
  end if

!--------------------------------------------------------
! nidging to 3D data, for example to ERA5

  uninitialized = .false.
  if(.not.allocated(uobs)) then
      allocate(uobs(nx+1,ny,nzm,2))
      allocate(vobs(nx,ny+YES3D,nzm,2))
      allocate(wobs(nx,ny,nz,2))
      allocate(tobs(nx,ny,nzm,2))
      allocate(qobs(nx,ny,nzm,2))
      uninitialized = .true.
  end if
  if(donudge3D) then
    if(docheck.or.uninitialized.or.(dayfld3Dobs(1).eq.0.or.dayfld3Dobs(2).gt.0..and.day.gt.dayfld3Dobs(2)).and.icycle.eq.1) &
          call read_nudging(nudge3D_dir,nudge3D_file,uobs,vobs,wobs,tobs,qobs,dayfld3Dobs)
    if(dayfld3Dobs(1).gt.0.) then
      coef=(day-dayfld3Dobs(1))/(dayfld3Dobs(2)-dayfld3Dobs(1))
    else
      coef = 0.
    end if
  else
   coef = 0.
  end if
  if(docheck) return
!====================================================================================
! Nudging U

  if(donudging_uv.and.nstep.ge.nudge3Dstep_start.and.nstep.le.nudge3Dstep_end) then
    do k=1,nzm
     if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        taum(k) = 1./nudge3D_tau
     else
        taum(k) = 0. 
     end if
    end do
    if(.not.allocated(tend_in_u)) allocate(tend_in_u(nx,ny,nzm))
    if(.not.dospectralnudging.or.dospectralnudging.and. &
                (icycle.eq.1.and.mod(nstep-1,nstep_spectral).eq.0)) then
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            fld = uobs(i,j,k,1)+(uobs(i,j,k,2)-uobs(i,j,k,1))*coef
            tend_in_u(i,j,k) = (fld-u(i,j,k))*taum(k)*terrau(i,j,k)
          end do
        end do
      end do 
      if(dospectralnudging) call box_smooth_3d(tend_in_u,terrau(1:nx,1:ny,1:nzm),&
           nx,ny,nzm,nsubdomains_x,nsubdomains_y,nx_spectral,ny_spectral,comm,kblock)
     end if
     tend_in(:,:,:) = tend_in_u(:,:,:)
  else
    taum = 0.
    tend_in = 0.
  end if
  do k=1,nzm
   do j=1,ny
     do i=1,nx
       if(taum(k).gt.max(tau_nudge3Du(i,j),tauz(k))) then
        tend = tend_in(i,j,k)
       else
        fld = uobs(i,j,k,1)+(uobs(i,j,k,2)-uobs(i,j,k,1))*coef
        tend = (fld-u(i,j,k))*max(tau_nudge3Du(i,j),tauz(k))*terrau(i,j,k)
       end if
       dudt(i,j,k,na)=dudt(i,j,k,na)+tend
     end do
   end do
  end do
!====================================================================================
! Nudging V

  if(donudging_uv.and.nstep.ge.nudge3Dstep_start.and.nstep.le.nudge3Dstep_end) then
    do k=1,nzm
     if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        taum(k) = 1./nudge3D_tau
     else
        taum(k) = 0. 
     end if
    end do
    if(.not.allocated(tend_in_v)) allocate(tend_in_v(nx,ny,nzm))
    if(.not.dospectralnudging.or.dospectralnudging.and. &
                (icycle.eq.1.and.mod(nstep-1,nstep_spectral).eq.0)) then
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            fld = vobs(i,j,k,1)+(vobs(i,j,k,2)-vobs(i,j,k,1))*coef
            tend_in_v(i,j,k) = (fld-v(i,j,k))*taum(k)*terrav(i,j,k)
          end do
        end do
      end do 
      if(dospectralnudging) call box_smooth_3d(tend_in_v,terrav(1:nx,1:ny,1:nzm),&
          nx,ny,nzm,nsubdomains_x,nsubdomains_y,nx_spectral,ny_spectral,comm,kblock)
     end if
     tend_in(:,:,:) = tend_in_v(:,:,:)
  else
    taum = 0.
    tend_in = 0.
  end if
  do k=1,nzm
   do j=1,ny
     do i=1,nx
       if(taum(k).gt.max(tau_nudge3Dv(i,j),tauz(k))) then
        tend = tend_in(i,j,k)
       else
        fld = vobs(i,j,k,1)+(vobs(i,j,k,2)-vobs(i,j,k,1))*coef
        tend = (fld-v(i,j,k))*max(tau_nudge3Dv(i,j),tauz(k))*terrav(i,j,k)
       end if
       dvdt(i,j,k,na)=dvdt(i,j,k,na)+tend
     end do
   end do
  end do
!====================================================================================
! Nudging W

  if(donudging_w.and.nstep.ge.nudge3Dstep_start.and.nstep.le.nudge3Dstep_end) then
    do k=1,nzm
     if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        taum(k) = 1./nudge3D_tau
     else
        taum(k) = 0.
     end if
    end do
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          fld = wobs(i,j,k,1)+(wobs(i,j,k,2)-wobs(i,j,k,1))*coef
          tend_in(i,j,k) = (fld-w(i,j,k))*taum(k)*terraw(i,j,k)
        end do
      end do
    end do
    if(dospectralnudging) call box_smooth_3d(tend_in,terraw(1:nx,1:ny,1:nzm),&
        nx,ny,nzm,nsubdomains_x,nsubdomains_y,nx_spectral,ny_spectral,comm,kblock)
  else
    taum = 0.
    tend_in = 0.
  end if
  do k=1,nzm
    do j=1,ny
      do i=1,nx
       if(taum(k).gt.max(tau_nudge3D(i,j),tauz(k))) then
        tend = tend_in(i,j,k)
       else
         fld = wobs(i,j,k,1)+(wobs(i,j,k,2)-wobs(i,j,k,1))*coef
         tend = (fld-w(i,j,k))*max(tau_nudge3D(i,j),tauz(k))*terraw(i,j,k)
       end if
       dwdt(i,j,k,na)=dwdt(i,j,k,na)+tend
      end do
    end do
  end do

!====================================================================================
! nudge velocities normal to the lateral boundaries for doregion=.true.
! Don;t just specify them, nudge over a small time scale

  if(doregion) then
    call task_rank_to_index(rank,it,jt)
    if(it.eq.0) then
      do k=1,nzm
       do j=1,ny
          fld = uobs(1,j,k,1)+(uobs(1,j,k,2)-uobs(1,j,k,1))*coef
          u_w(j,k)=u_w(j,k) - dtn*(u_w(j,k)-fld)*tau_nudge3Du(1,j)*terrau(1,j,k)
       end do
      end do
    end if
    if(it+nx.eq.nx_gl) then
      do k=1,nzm
       do j=1,ny
          fld = uobs(nx+1,j,k,1)+(uobs(nx+1,j,k,2)-uobs(nx+1,j,k,1))*coef
          u_e(j,k)=u_e(j,k) - dtn*(u_e(j,k)-fld)*tau_nudge3Du(nx+1,j)*terrau(nx+1,j,k)
       end do
      end do
    end if
    if(jt.eq.0) then
      do k=1,nzm
       do i=1,nx
          fld = vobs(i,1,k,1)+(vobs(i,1,k,2)-vobs(i,1,k,1))*coef
          v_s(i,k)=v_s(i,k) - dtn*(v_s(i,k)-fld)*tau_nudge3Dv(i,1)*terrav(i,1,k)
       end do
      end do
    end if
    if(jt+ny.eq.ny_gl) then
      do k=1,nzm
       do i=1,nx
          fld = vobs(i,ny+YES3D,k,1)+(vobs(i,ny+YES3D,k,2)-vobs(i,ny+YES3D,k,1))*coef
          v_n(i,k)=v_n(i,k) - dtn*(v_n(i,k)-fld)*tau_nudge3Dv(i,ny+YES3D)*terrav(i,ny+YES3D,k)
       end do
      end do
    end if
   end if ! doregion
!====================================================================================
! Nudging T

    if((donudging_tq.or.donudging_t) &
          .and.nstep.ge.nudge3Dstep_start.and.nstep.le.nudge3Dstep_end) then
      do k=1,nzm
        if(z(k).ge.nudging_t_z1.and.z(k).le.nudging_t_z2) then
           taum(k) = 1./nudge3D_tau
        else
           taum(k) = 0.
        end if
      end do
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            fld = tobs(i,j,k,1)+(tobs(i,j,k,2)-tobs(i,j,k,1))*coef
            tend_in(i,j,k) = (fld-t(i,j,k))*taum(k)
          end do
        end do
      end do
      if(dospectralnudging) call box_smooth_3d(tend_in,terra(1:nx,1:ny,1:nzm),&
          nx,ny,nzm,nsubdomains_x,nsubdomains_y,nx_spectral,ny_spectral,comm,kblock)
    else
      taum = 0.
      tend_in = 0.
    end if
    do k=1,nzm
     do j=1,ny
       do i=1,nx
         if(taum(k).gt.max(tau_nudge3D(i,j),tauz(k))) then
          tend = tend_in(i,j,k)
         else
          fld = tobs(i,j,k,1)+(tobs(i,j,k,2)-tobs(i,j,k,1))*coef
          tend = (fld+gamaz(k)-t(i,j,k))*max(tau_nudge3D(i,j),tauz(k))
         end if
         t(i,j,k)=t(i,j,k)+tend*dtn
       end do
     end do
    end do
  
!====================================================================================
! Nudging Q

    if((donudging_tq.or.donudging_q) &
          .and.nstep.ge.nudge3Dstep_start.and.nstep.le.nudge3Dstep_end) then
      do k=1,nzm
        if(z(k).ge.nudging_q_z1.and.z(k).le.nudging_q_z2) then
           taum(k) = 1./nudge3D_tau
        else
           taum(k) = 0.
        end if
      end do
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            fld = qobs(i,j,k,1)+(qobs(i,j,k,2)-qobs(i,j,k,1))*coef
            tend_in(i,j,k) = (fld-micro_field(i,j,k,index_water_vapor))*taum(k)
          end do
        end do
      end do
      if(dospectralnudging) call box_smooth_3d(tend_in,terra(1:nx,1:ny,1:nzm),&
          nx,ny,nzm,nsubdomains_x,nsubdomains_y,nx_spectral,ny_spectral,comm,kblock)
    else
      taum = 0.
      tend_in = 0.
    end if
    do k=1,nzm
     do j=1,ny
       do i=1,nx
         if(taum(k).gt.max(tau_nudge3D(i,j),tauz(k))) then
          tend = tend_in(i,j,k)
         else
          fld = qobs(i,j,k,1)+(qobs(i,j,k,2)-qobs(i,j,k,1))*coef
          tend = (fld-micro_field(i,j,k,index_water_vapor))*max(tau_nudge3D(i,j),tauz(k))
         end if
         micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor)+tend*dtn
       end do
     end do
    end do

else

!==================================================================
! nidging to 1D data in snd file

  coef = 1./tauls

  if(donudging_uv) then
    unudge = 0.
    vnudge = 0.
    do k=1,nzm
      if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        unudge(k)=unudge(k) - (u0(k)-ug0(k))*coef
        vnudge(k)=vnudge(k) - (v0(k)-vg0(k))*coef
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ug0(k))*coef*terrau(i,j,k)
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))*coef*terrav(i,j,k)
          end do
        end do
      end if
    end do
  endif

  coef = 1./tautqls

  if(donudging_tq.or.donudging_t) then
    coef1 = dtn / tautqls
    tnudge = 0.
    do k=1,nzm
      if(z(k).ge.nudging_t_z1.and.z(k).le.nudging_t_z2) then
        tnudge(k)=tnudge(k) -(t0(k)-tg0(k)-gamaz(k))*coef
        do j=1,ny
          do i=1,nx
             t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*coef1*terra(i,j,k)
          end do
        end do
      end if
    end do
  endif

  if(donudging_tq.or.donudging_q) then
    qnudge = 0.
    coef1 = dtn / tautqls
    do k=1,nzm
      if(z(k).ge.nudging_q_z1.and.z(k).le.nudging_q_z2) then
        qnudge(k)=qnudge(k) -(q0(k)-qg0(k))*coef
        do j=1,ny
          do i=1,nx
             micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor) &
                                -(q0(k)-qg0(k))*coef1*terra(i,j,k)
          end do
        end do
      end if
    end do
  endif

  if(docolumn.and.(dobufferzonex.or.dobufferzoney)) then
   if(masterproc) print*,'ERROR: dobufferzonex dobufferzoney cannot be set to T when docolumn=T! Exitting...'
   call task_abort()
  end if
  tau(:) = 1.e6
  if(dobufferzonex) then
    call task_rank_to_index(rank,it,jt)
    do k=1,nzm
       tau(k) = min(tau(k),0.3*bufferzonex*nx_gl*dx/(abs(u0(k))+1.))
       coef = 1./tau(k)
       coef1 = dtn / tau(k)
       do j=1,ny
           do i=1,nx
             xxx = dx*(i+it-1)
             if(xxx.le.bufferzonex/2.*nx_gl*dx.or.xxx.ge.(1.-bufferzonex/2.)*nx_gl*dx) then
               dudt(i,j,k,na)=dudt(i,j,k,na)-(u(i,j,k)-ug0(k))*coef*terrau(i,j,k)
               dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v(i,j,k)-vg0(k))*coef*terrav(i,j,k)
               dwdt(i,j,k,na)=dwdt(i,j,k,na)-w(i,j,k)*coef*terraw(i,j,k)
               t(i,j,k)=t(i,j,k)-(t(i,j,k)-tg0(k)-gamaz(k))*coef1*terra(i,j,k)
               micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)- &
                     (micro_field(i,j,k,index_water_vapor)-qg0(k))*coef1*terra(i,j,k)
               micro_field(i,j,k,1:index_water_vapor-1)= &
                      micro_field(i,j,k,1:index_water_vapor-1)- &
                      micro_field(i,j,k,1:index_water_vapor-1)*coef1*terra(i,j,k)
               micro_field(i,j,k,index_water_vapor+1:nmicro_fields)= &
                      micro_field(i,j,k,index_water_vapor+1:nmicro_fields)- &
                      micro_field(i,j,k,index_water_vapor+1:nmicro_fields)*coef1*terra(i,j,k)
               sgs_field(i,j,k,:)=sgs_field(i,j,k,:)- sgs_field(i,j,k,:)*coef1*terra(i,j,k)
               if(dotracers) tracer(i,j,k,:) = 0.
             end if
           end do
       end do
    end do
  end if
  if(RUN3D.and.dobufferzoney) then
      call task_rank_to_index(rank,it,jt)
      tau(:) = 1.e6
      do k=1,nzm
         tau(k) = min(tau(k),0.3*bufferzoney*(y_gl(ny_gl)-y_gl(1))/(abs(v0(k))+1.))
         coef = 1./tau(k)
         coef1 = dtn / tau(k)
         do j=1,ny
           yyy = y_gl(j+jt)-y_gl(1)
           if(yyy.le.bufferzoney/2.*(yv_gl(ny_gl+1)-yv_gl(1)) &
              .or.yyy.ge.(1.-bufferzoney/2.)*(yv_gl(ny_gl+1)-yv_gl(1))) then
             do i=1,nx
                 dudt(i,j,k,na)=dudt(i,j,k,na)-(u(i,j,k)-ug0(k))*coef*terrau(i,j,k)
                 dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v(i,j,k)-vg0(k))*coef*terrav(i,j,k)
                 dwdt(i,j,k,na)=dwdt(i,j,k,na)-w(i,j,k)*coef*terraw(i,j,k)
                 t(i,j,k)=t(i,j,k)-(t(i,j,k)-tg0(k)-gamaz(k))*coef1*terra(i,j,k)
                 micro_field(i,j,k,index_water_vapor)= &
                      micro_field(i,j,k,index_water_vapor)- &
                     (micro_field(i,j,k,index_water_vapor)-qg0(k))*coef1*terra(i,j,k)
                 micro_field(i,j,k,1:index_water_vapor-1)= &
                      micro_field(i,j,k,1:index_water_vapor-1)- &
                      micro_field(i,j,k,1:index_water_vapor-1)*coef1*terra(i,j,k)
                 micro_field(i,j,k,index_water_vapor+1:nmicro_fields)= &
                      micro_field(i,j,k,index_water_vapor+1:nmicro_fields)- &
                      micro_field(i,j,k,index_water_vapor+1:nmicro_fields)*coef1*terra(i,j,k)
                 sgs_field(i,j,k,:)=sgs_field(i,j,k,:)- sgs_field(i,j,k,:)*coef1*terra(i,j,k)
                 if(dotracers) tracer(i,j,k,:) = 0.
             end do
           end if
         end do
      end do
  endif

end if

call sumenergy(-1,'nudging')

call t_stopf('nudging')

end subroutine nudging
