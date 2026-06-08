
subroutine boundaries(flag)

! exchange boundary informations between subdomains        	

use vars
use microphysics
use chemistry
use chemistry_params, only : do_iepox_aero_chem, do_iepox_droplet_chem
use chem_aqueous, only: naqchem_fields
use chem_aerosol, only: narchem_fields
use sgs
use params, only: dotracers, dosgs, dowallx, dowally, doregion, docolumn, dodebug, docheck, dochem
use tracers
implicit none
integer i,j,k,i1,i2,j1,j2,it,jt,im
real factor
integer flag
real tmp(1:ny,1:nzm)

if(docheck) return

if(docolumn) return

call t_startf ('boundaries')

if(flag.eq.1) then

 if(dompi) then
   call task_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV)
   call task_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3)
   call task_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV)
 else
   call bound_exchange(u,dimx1_u,dimx2_u,dimy1_u,dimy2_u,nzm,2,3,2+NADV,2+NADV)
   call bound_exchange(v,dimx1_v,dimx2_v,dimy1_v,dimy2_v,nzm,2+NADV,2+NADV,2,3)
   call bound_exchange(w,dimx1_w,dimx2_w,dimy1_w,dimy2_w,nz,2+NADV,2+NADV,2+NADV,2+NADV)
 end if

 if(doregion) then
  if(dompi) then
    call task_exchange_2D(u_e, dimy1_u, dimy2_u, nzm, 2+NADV,2+NADV, 2)
    call task_exchange_2D(v_n, dimx1_v, dimx2_v, nzm, 2+NADV,2+NADV, 1)
  end if
 end if

 if(dowallx) then
  if(mod(rank,nsubdomains_x).eq.0) then
    if(doregion) then
      do i=dimx1_u,0
       u(i,:,:) = u(1,:,:)
      end do
      do i=dimx1_w,0
       w(i,:,:) = w(1,:,:) 
      end do
    else
     if(dodns) then
      u(:1,:,:) = 0.
      w(:0,:,:) = 0. 
     else
      u(1,:,:) = 0.
      u(dimx1_u:0,:,:) = -u(-dimx1_u+2:2:-1,:,:)
      w(dimx1_w:0,:,:) = w(-dimx1_w+1:1:-1,:,:) 
     end if
    end if
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    if(doregion) then
      u(nx+1,:,1:nzm) = u_e(:,1:nzm)
      do i=nx+2,dimx2_u
       u(i,:,:) = u(nx+1,:,:)
      end do
      do i=nx+1,dimx2_w
       w(i,:,:) = w(nx,:,:) 
      end do
    else
     if(dodns) then
      u(nx+1:,:,:) = 0.
      w(nx+1:,:,:) = 0. 
     else 
      u(nx+1,:,:) = 0.
      u(nx+2:dimx2_u,:,:) = -u(nx:2*(nx+1)-dimx2_u:-1,:,:)
      w(nx+1:dimx2_w,:,:) =  w(nx:2*nx-dimx2_w+1:-1,:,:)
     end if
    end if
  end if
 end if
 if(RUN3D.and.dowally) then
  if(rank.lt.nsubdomains_x) then
    if(doregion) then
      do j=dimy1_v,0
       v(:,j,:) = v(:,1,:)
      end do
      do j=dimy1_u,0
       u(:,j,:) = u(:,1,:)
      end do
      do j=dimy1_w,0
       w(:,j,:) = w(:,1,:) 
      end do
    else
     if(dodns) then
      v(:,:1,:) = 0.
      u(:,:0,:) = 0.
      w(:,:0,:) = 0.
     else
      v(:,1,:) = 0.
      v(:,dimy1_v:0,:) = -v(:,-dimy1_v+2:2:-1,:)
      u(:,dimy1_u:0,:) = u(:,-dimy1_u+1:1:-1,:)
      w(:,dimy1_w:0,:) = w(:,-dimy1_w+1:1:-1,:)
     end if
    end if
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
    if(doregion) then
      v(:,ny+YES3D,:) = v_n(:,:)
      do j=ny+2,dimy2_v
       v(:,j,:) = v(:,ny+YES3D,:)
      end do
      do j=ny+YES3D,dimy2_u
       u(:,j,:) = u(:,ny,:)
      end do
      do j=ny+YES3D,dimy2_w
       w(:,j,:) = w(:,ny,:)
      end do
    else
     if(dodns) then 
      v(:,ny+YES3D:,:) = 0.
      u(:,ny+YES3D:,:) = 0.
      w(:,ny+YES3D:,:) = 0.
     else
      v(:,ny+YES3D,:) = 0
      v(:,ny+2*YES3D:dimy2_v,:) = -v(:,ny:2*(ny+YES3D)-dimy2_v:-1,:)
      u(:,ny+YES3D:dimy2_u,:) = u(:,ny:2*ny-dimy2_u+YES3D:-1,:)
      w(:,ny+YES3D:dimy2_w,:) = w(:,ny:2*ny-dimy2_w+YES3D:-1,:)
     end if
    end if
  end if
 end if
 if(dowallx) then
  if(mod(rank,nsubdomains_x).eq.0) then
    if(doregion) then
      do i=dimx1_v,0
       v(i,:,:) = v(1,:,:)
      end do
    else
     if(dodns) then
      v(:0,:,:) = 0.
     else
      v(dimx1_v:0,:,:) = v(-dimx1_v+1:1:-1,:,:)
     end if
    end if
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    if(doregion) then
      do i=nx+1,dimx2_v
       v(i,:,:) = v(nx,:,:)
      end do
    else
     if(dodns) then 
      v(nx+1:,:,:) = 0.
     else
      v(nx+1:dimx2_v,:,:) =  v(nx:2*nx-dimx2_v+1:-1,:,:)
     end if 
    end if
  end if
 end if

end if

if(flag.eq.2.or.flag.eq.3) then

 if(flag.eq.2) then
  i1 = 3+NADVS
  i2 = 3+NADVS
  j1 = 3+NADVS
  j2 = 3+NADVS
 else
  i1 = 1
  i2 = 1
  j1 = 1
  j2 = 1
 end if

 if(dompi) then
   call task_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
 else
   call bound_exchange(t,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
 end if

 if(dowallx) then
  if(mod(rank,nsubdomains_x).eq.0) then
   if(doregion) then
    do i=dimx1_s,0
     t(i,:,:) = t(1,:,:)
    end do
   else
    if(.not.dodns) t(dimx1_s:0,:,:) = t(-dimx1_s+1:1:-1,:,:)
   end if
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
   if(doregion) then
    do i=nx+1,dimx2_s
     t(i,:,:) = t(nx,:,:)
    end do
   else
    if(.not.dodns) t(nx+1:dimx2_s,:,:) =  t(nx:2*nx-dimx2_s+1:-1,:,:)
   end if
  end if
 end if
 if(RUN3D.and.dowally) then
  if(rank.lt.nsubdomains_x) then
   if(doregion) then
    do j=dimy1_s,0
     t(:,j,:) = t(:,1,:)
    end do
   else
    if(.not.dodns) t(:,dimy1_s:0,:) = t(:,-dimy1_s+1:1:-1,:)
   end if
  end if
  if(rank.gt.nsubdomains-nsubdomains_x-1) then
   if(doregion) then
    do j=ny+1,dimy2_s
     t(:,j,:) = t(:,ny,:)
    end do
   else
    if(.not.dodns) t(:,ny+1:dimy2_s,:) = t(:,ny:2*ny-dimy2_s+1:-1,:)
   end if
  end if
 end if

 if(dosgs.and.advect_sgs) then
   do im = 1,nsgs_fields
     if(dompi) then
      call task_exchange(sgs_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     else
      call bound_exchange(sgs_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     end if
     if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
       if(doregion) then
        do i=dimx1_s,0
         sgs_field(i,:,:,im) = sgs_field(1,:,:,im)
        end do
       else
         sgs_field(dimx1_s:0,:,:,im) = sgs_field(-dimx1_s+1:1:-1,:,:,im)
       end if
       end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
       if(doregion) then
         do i=nx+1,dimx2_s
          sgs_field(i,:,:,im) = sgs_field(nx,:,:,im)
         end do
       else
        sgs_field(nx+1:dimx2_s,:,:,im) =  sgs_field(nx:2*nx-dimx2_s+1:-1,:,:,im)
       end if
      end if
     end if
     if(RUN3D.and.dowally) then
      if(rank.lt.nsubdomains_x) then
       if(doregion) then
        do j=dimy1_s,0
         sgs_field(:,j,:,im) = sgs_field(:,1,:,im)
        end do
       else
        sgs_field(:,dimy1_s:0,:,im) = sgs_field(:,-dimy1_s+1:1:-1,:,im)
       end if
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
       if(doregion) then
        do j=ny+1,dimy2_s
         sgs_field(:,j,:,im) = sgs_field(:,ny,:,im)
        end do
       else
        sgs_field(:,ny+1:dimy2_s,:,im) = sgs_field(:,ny:2*ny-dimy2_s+1:-1,:,im)
       end if
      end if
     end if
   end do
 end if

 do im = 1,nmicro_fields
  if(flag_advect(im).eq.1) then
    if(   im.eq.index_water_vapor             &
     .or. docloud.and.flag_precip(im).ne.1    &
     .or. doprecip.and.flag_precip(im).eq.1 ) then
     if(dompi) then
       call task_exchange(micro_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     else
       call bound_exchange(micro_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     end if
     if(im.eq.index_water_vapor) then
       factor = 1.
     else
       factor = 0.
     end if
     if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
       if(doregion) then
        do i=dimx1_s,0
         micro_field(i,:,:,im) = micro_field(1,:,:,im)*factor
        end do
       else
         if(.not.dodns) micro_field(dimx1_s:0,:,:,im) = micro_field(-dimx1_s+1:1:-1,:,:,im)
       end if
       end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
       if(doregion) then
         do i=nx+1,dimx2_s
          micro_field(i,:,:,im) = micro_field(nx,:,:,im)*factor
         end do
       else
        if(.not.dodns) micro_field(nx+1:dimx2_s,:,:,im) =  micro_field(nx:2*nx-dimx2_s+1:-1,:,:,im)
       end if
      end if
     end if
     if(RUN3D.and.dowally) then
      if(rank.lt.nsubdomains_x) then
       if(doregion) then
        do j=dimy1_s,0
         micro_field(:,j,:,im) = micro_field(:,1,:,im)*factor
        end do
       else
        if(.not.dodns) micro_field(:,dimy1_s:0,:,im) = micro_field(:,-dimy1_s+1:1:-1,:,im)
       end if
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
       if(doregion) then
        do j=ny+1,dimy2_s
         micro_field(:,j,:,im) = micro_field(:,ny,:,im)*factor
        end do
       else
         if(.not.dodns) micro_field(:,ny+1:dimy2_s,:,im) = micro_field(:,ny:2*ny-dimy2_s+1:-1,:,im)
       end if
      end if
     end if
    end if
  end if
 end do

 if(dotracers) then
   do im=1,ntracers
     if(dompi) then
      call task_exchange(tracer(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     else
      call bound_exchange(tracer(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     end if
     if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
       if(doregion) then
        do i=dimx1_s,0
         tracer(i,:,:,im) = 0.
        end do
       else
         tracer(dimx1_s:0,:,:,im) = tracer(-dimx1_s+1:1:-1,:,:,im)
       end if
       end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
       if(doregion) then
         do i=nx+1,dimx2_s
          tracer(i,:,:,im) = 0.
         end do
       else
        tracer(nx+1:dimx2_s,:,:,im) =  tracer(nx:2*nx-dimx2_s+1:-1,:,:,im)
       end if
      end if
     end if
     if(RUN3D.and.dowally) then
      if(rank.lt.nsubdomains_x) then
       if(doregion) then
        do j=dimy1_s,0
         tracer(:,j,:,im) = 0.
        end do
       else
        tracer(:,dimy1_s:0,:,im) = tracer(:,-dimy1_s+1:1:-1,:,im)
       end if
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
       if(doregion) then
        do j=ny+1,dimy2_s
         tracer(:,j,:,im) = 0.
        end do
       else
        tracer(:,ny+1:dimy2_s,:,im) = tracer(:,ny:2*ny-dimy2_s+1:-1,:,im)
       end if
      end if
     end if
   end do
 end if

 if( dochem ) then
   do im = 1, ngchem_fields
     if(dompi) then
      call task_exchange(gchem_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     else
      call bound_exchange(gchem_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     end if
     if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
       if(doregion) then
        do i=dimx1_s,0
         gchem_field(i,:,:,im) = 0.
        end do
       else
         gchem_field(dimx1_s:0,:,:,im) = gchem_field(-dimx1_s+1:1:-1,:,:,im)
       end if
       end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
       if(doregion) then
         do i=nx+1,dimx2_s
          gchem_field(i,:,:,im) = 0.
         end do
       else
        gchem_field(nx+1:dimx2_s,:,:,im) =  gchem_field(nx:2*nx-dimx2_s+1:-1,:,:,im)
       end if
      end if
     end if
     if(RUN3D.and.dowally) then
      if(rank.lt.nsubdomains_x) then
       if(doregion) then
        do j=dimy1_s,0
         gchem_field(:,j,:,im) = 0.
        end do
       else
        gchem_field(:,dimy1_s:0,:,im) = gchem_field(:,-dimy1_s+1:1:-1,:,im)
       end if
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
       if(doregion) then
        do j=ny+1,dimy2_s
         gchem_field(:,j,:,im) = 0.
        end do
       else
        gchem_field(:,ny+1:dimy2_s,:,im) = gchem_field(:,ny:2*ny-dimy2_s+1:-1,:,im)
       end if
      end if
     end if
   end do
 end if

 if( dochem .and. do_iepox_droplet_chem ) then
   do im = 1, naqchem_fields
     if(dompi) then
      call task_exchange(aqchem_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
      call task_exchange(aqchem_gasprod_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)

     else
      call bound_exchange(aqchem_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
      call bound_exchange(aqchem_gasprod_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     end if
     if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
       if(doregion) then
        do i=dimx1_s,0
         aqchem_field(i,:,:,im) = 0.
         aqchem_gasprod_field(i,:,:,im) = 0.
        end do
       else
         aqchem_field(dimx1_s:0,:,:,im) = aqchem_field(-dimx1_s+1:1:-1,:,:,im)
         aqchem_gasprod_field(dimx1_s:0,:,:,im) = aqchem_gasprod_field(-dimx1_s+1:1:-1,:,:,im)
       end if
       end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
       if(doregion) then
         do i=nx+1,dimx2_s
          aqchem_field(i,:,:,im) = 0.
          aqchem_gasprod_field(i,:,:,im) = 0.
         end do
       else
        aqchem_field(nx+1:dimx2_s,:,:,im) =  aqchem_field(nx:2*nx-dimx2_s+1:-1,:,:,im)
        aqchem_gasprod_field(nx+1:dimx2_s,:,:,im) =  aqchem_gasprod_field(nx:2*nx-dimx2_s+1:-1,:,:,im)
       end if
      end if
     end if
     if(RUN3D.and.dowally) then
      if(rank.lt.nsubdomains_x) then
       if(doregion) then
        do j=dimy1_s,0
         aqchem_field(:,j,:,im) = 0.
         aqchem_gasprod_field(:,j,:,im) = 0.
        end do
       else
        aqchem_field(:,dimy1_s:0,:,im) = aqchem_field(:,-dimy1_s+1:1:-1,:,im)
        aqchem_gasprod_field(:,dimy1_s:0,:,im) = aqchem_gasprod_field(:,-dimy1_s+1:1:-1,:,im)
       end if
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
       if(doregion) then
        do j=ny+1,dimy2_s
         aqchem_field(:,j,:,im) = 0.
         aqchem_gasprod_field(:,j,:,im) = 0.
        end do
       else
        aqchem_field(:,ny+1:dimy2_s,:,im) = aqchem_field(:,ny:2*ny-dimy2_s+1:-1,:,im)
        aqchem_gasprod_field(:,ny+1:dimy2_s,:,im) = aqchem_gasprod_field(:,ny:2*ny-dimy2_s+1:-1,:,im)
       end if
      end if
     end if
   end do
 end if

 if( dochem .and. do_iepox_aero_chem ) then
   do im = 1, narchem_fields
     if(dompi) then
      call task_exchange(archem_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     else
      call bound_exchange(archem_field(:,:,:,im),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,i1,i2,j1,j2)
     end if
     if(dowallx) then
      if(mod(rank,nsubdomains_x).eq.0) then
       if(doregion) then
        do i=dimx1_s,0
         archem_field(i,:,:,im) = 0.
        end do
       else
         archem_field(dimx1_s:0,:,:,im) = archem_field(-dimx1_s+1:1:-1,:,:,im)
       end if
       end if
      if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
       if(doregion) then
         do i=nx+1,dimx2_s
          archem_field(i,:,:,im) = 0.
         end do
       else
        archem_field(nx+1:dimx2_s,:,:,im) =  archem_field(nx:2*nx-dimx2_s+1:-1,:,:,im)
       end if
      end if
     end if
     if(RUN3D.and.dowally) then
      if(rank.lt.nsubdomains_x) then
       if(doregion) then
        do j=dimy1_s,0
         archem_field(:,j,:,im) = 0.
        end do
       else
        archem_field(:,dimy1_s:0,:,im) = archem_field(:,-dimy1_s+1:1:-1,:,im)
       end if
      end if
      if(rank.gt.nsubdomains-nsubdomains_x-1) then
       if(doregion) then
        do j=ny+1,dimy2_s
         archem_field(:,j,:,im) = 0.
        end do
       else
        archem_field(:,ny+1:dimy2_s,:,im) = archem_field(:,ny:2*ny-dimy2_s+1:-1,:,im)
       end if
      end if
     end if
   end do
 end if

endif


if(flag.eq.4) then

 if(dosgs) then
   do i = 1,nsgs_fields_diag
    if(dompi) then
       call task_exchange(sgs_field_diag(:,:,:,i),dimx1_d,dimx2_d,dimy1_d,dimy2_d,nzm, &
                   1-dimx1_d,dimx2_d-nx,1-YES3D+1-dimy1_d,1-YES3D+dimy2_d-ny)
    else
       call bound_exchange(sgs_field_diag(:,:,:,i),dimx1_d,dimx2_d,dimy1_d,dimy2_d,nzm, &
                   1-dimx1_d,dimx2_d-nx,1-YES3D+1-dimy1_d,1-YES3D+dimy2_d-ny)
    end if 
    if(dowallx) then
     if(mod(rank,nsubdomains_x).eq.0) then
       sgs_field_diag(dimx1_d:0,:,:,i) = sgs_field_diag(-dimx1_d+1:1:-1,:,:,i)
     end if
     if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
       sgs_field_diag(nx+1:dimx2_d,:,:,i) =  sgs_field_diag(nx:2*nx-dimx2_d+1:-1,:,:,i)
     end if
    end if
    if(RUN3D.and.dowally) then
     if(rank.lt.nsubdomains_x) then
       sgs_field_diag(:,dimy1_d:0,:,i) = sgs_field_diag(:,-dimy1_d+1:1:-1,:,i)
     end if
     if(rank.gt.nsubdomains-nsubdomains_x-1) then
       sgs_field_diag(:,ny+1:dimy2_d,:,i) = sgs_field_diag(:,ny:2*ny-dimy2_d+1:-1,:,i)
     end if
    end if
   end do
 end if

end if

call t_stopf ('boundaries')

end subroutine boundaries
	
	
