subroutine diffuse_damping_mom_z
	
!        momentum SGS diffusion in the vertical

use vars
use sgs
use params, only: docolumn, dowallx, dowally, dosgs, doimplicitdiff, doregion, &
                  z0, dodamping, nub, dodebug, dodamping_w, damping_w_cu, &
                  dodamping_u, dolatlon, dodamping_poles, damping_u_cu, doflat
use terrain
use check_energy
implicit none

real rdz2,rdz,rdx2u,rdx2v,rdx2w,rdy2u,rdy2v,rdy2w,rdz25

integer i,j,k,ib,jb,kb,kc,kcu
real(8) tkx, tky, tkz, iadzw, iadz, a, b, c, d, e, wnd
real(8) alpha(nx,ny,nzm),beta(nx,ny,nzm)
real nu, tau_max, tauy(ny), umax(ny), tauz(nzm), tau_vel0, wmax, wmax1, tau
real vel0(nx,ny,nzm), tau_vel(nx,ny,nzm)
real, allocatable :: fuz(:,:,:),fvz(:,:,:)
real zzz, velmax(nx,ny), velmin(nx,ny)

!-----------------------------------------
if(.not.(dosgs.and.doimplicitdiff)) return

call t_startf('diffuse_damping_mom_z')

if(dodebug) then
  if(masterproc) print*,'before diffuse_mom_z:'
  call fminmax_print('u:',u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nz),1,nx,1,ny,nz)
end if

tau_max=1./dt
tau_vel0=1./dt

!------------------------------------------------
! for statistics (fluxes)
if (dostatis) then
 rdz=1./dz
 allocate(fuz(nx,ny,nz),fvz(nx,ny,nz))
 do k=1,nzm-1
  kc=k+1
  rdz2 = rdz*rdz*rho(k)/adz(k)
  rdz25 = rdz*rdz*rhow(kc)/adzw(kc)
   do j=1,ny
    jb=j-YES3D
    do i=1,nx
     ib=i-1
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
              /(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
     fuz(i,j,kc)=-0.5*rdz25*tkz*(u(i,j,kc)-u(i,j,k))*min(terrau(i,j,kc),terrau(i,j,k))
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
              /(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
     fvz(i,j,kc)=-0.5*rdz25*tkz*(v(i,j,kc)-v(i,j,k))*min(terrav(i,j,kc),terrav(i,j,k))
   end do
  end do
 end do
 fuz(:,:,1) = 0.
 fvz(:,:,1) = 0.
 do j=1,ny
  do i=1,nx
    fuz(i,j,k_terrau(i,j))=fluxbu(i,j) * rdz * rhow(k_terrau(i,j))
    fvz(i,j,k_terrav(i,j))=fluxbv(i,j) * rdz * rhow(k_terrav(i,j))
    tkz=rdz2*tk(i,j,nzm)
    fuz(i,j,nz)=fluxtu(i,j) * rdz * rhow(nz)
    fvz(i,j,nz)=fluxtv(i,j) * rdz * rhow(nz)
  end do
 end do
end if

!----------------------------------------------------

rdz=1./dz
rdz2 = dtn*rdz*rdz 

!----------------------------------------------------
! Diffusion/damping of u:

vel0(:,:,:) = 0.
tau_vel(:,:,:) = 0.
if(dolatlon.and..not.doflat.and.dodamping_poles) then
  do j=1,ny
    tauy(j) = tau_max*(1.-mu(j)**2)**200
    umax(j) = damping_u_cu*dx*mu(j)/dtn
  end do
  do k = 1,nzm
   if(dodamping_u.and.pres(k).lt.70.) then
    tau = tau_max
    do j=1,ny
      do i=1,nx
          if(terrau(i,j,k)*u(i,j,k).gt.umax(j)) then
            tau_vel(i,j,k)= tau
            vel0(i,j,k) = umax(j)
          else if(terrau(i,j,k)*u(i,j,k).lt.-umax(j)) then
            tau_vel(i,j,k)= tau
            vel0(i,j,k) = -umax(j)
          end if
      end do
    end do! j
   else
    do j=1,ny
      do i=1,nx
          if(terrau(i,j,k)*u(i,j,k).gt.umax(j)) then
            tau_vel(i,j,k)= tauy(j)
            vel0(i,j,k) = umax(j)
          else if(terrau(i,j,k)*u(i,j,k).lt.-umax(j)) then
            tau_vel(i,j,k)= tauy(j)
            vel0(i,j,k) = -umax(j)
          end if
      end do
    end do
   end if
  end do ! k
end if

do j=1,ny
 do i=1,nx
  k = k_terrau(i,j)
  ib=i-1
  kc =min(nzm,k+1)
  iadzw = rdz2*rhow(kc)/(adzw(kc)*adz(k)*rho(k))
  tkz = (tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
           +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
          /(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
  tkz = tkz*min(terrau(i,j,kc),terrau(i,j,k))
  c = -tkz*iadzw 
  b = 1. + dtn*tau_vel(i,j,k) - c
  d = u(i,j,k) + dtn*tau_vel(i,j,k)*vel0(i,j,k) &
               + dtn*rhow(k)/(dz*adz(k)*rho(k))*fluxbu(i,j)
  alpha(i,j,k) = -c/b
  beta(i,j,k) = d/b 
 end do
end do
do k=2,nzm-1
 kc=k+1
 kb=k-1
 iadz = rdz2*rhow(k)/(adzw(k)*adz(k)*rho(k))
 iadzw = rdz2*rhow(kc)/(adzw(kc)*adz(k)*rho(k))
 do j=1,ny
   do i=1,nx
    if(k.gt.k_terrau(i,j)) then
      ib=i-1
      tkz = (tk(i,j,kb)*terra(i,j,kb)+tk(ib,j,kb)*terra(ib,j,kb) &
              +tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k)) &
             /(terra(i,j,kb)+terra(ib,j,kb)+terra(i,j,k)+terra(ib,j,k)+1.e-10)
      tkz = tkz*min(terrau(i,j,kb),terrau(i,j,k))
      a = -tkz*iadz
      tkz = (tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
              +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
             /(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
      tkz = tkz*min(terrau(i,j,kc),terrau(i,j,k))
      c = -tkz*iadzw
      b = 1. + dtn*tau_vel(i,j,k) - a - c
      d = u(i,j,k)+ dtn*tau_vel(i,j,k)*vel0(i,j,k) 
      e = b + a*alpha(i,j,k-1)
      alpha(i,j,k) = -c/e
      beta(i,j,k) = (d-a*beta(i,j,k-1))/e
    end if
   end do
 end do
end do 
k = nzm
kb = k-1
kc = k+1
do j=1,ny
 do i=1,nx
  ib=i-1
  tkz = (tk(i,j,kb)*terra(i,j,kb)+tk(ib,j,kb)*terra(ib,j,kb) &
          +tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k)) &
         /(terra(i,j,kb)+terra(ib,j,kb)+terra(i,j,k)+terra(ib,j,k)+1.e-10)
  tkz = tkz*min(terrau(i,j,kb),terrau(i,j,k))
  iadz = rdz2*rhow(k)/(adzw(k)*adz(k)*rho(k))
  a = -tkz*iadz
  b = 1. + dtn*tau_vel(i,j,k) - a
  d = u(i,j,k)+ dtn*tau_vel(i,j,k)*vel0(i,j,k) &
              - dtn*rhow(kc)/(rho(k)*dz*adz(k))*fluxtu(i,j)
  e = b + a*alpha(i,j,k-1)
  u(i,j,k) = (d-a*beta(i,j,k-1))/e
 end do
end do
do k=nzm-1,1,-1
 kc=k+1
 do j=1,ny
   do i=1,nx
    if(k.ge.k_terrau(i,j)) then
      u(i,j,k) = alpha(i,j,k)*u(i,j,kc) + beta(i,j,k)
    end if
   end do
 end do
end do 

! -----------------------------------------
! Diffusion/damping of v:
  
vel0(:,:,:) = 0.
tau_vel(:,:,:) = 0.
if(dolatlon.and..not.doflat.and.dodamping_poles) then
  do j=1,ny
    tauy(j) = tau_max*(1.-mu(j)**2)**200
    umax(j) = damping_u_cu*dx*mu(j)/dtn
  end do
  do k = 1,nzm
   if(dodamping_u.and.pres(k).lt.70.) then
    tau = tau_max
    do j=1,ny
      do i=1,nx
          if(terrav(i,j,k)*v(i,j,k).gt.umax(j)) then
            tau_vel(i,j,k)= tau
            vel0(i,j,k) = umax(j)
          else if(terrav(i,j,k)*v(i,j,k).lt.-umax(j)) then
            tau_vel(i,j,k)= tau
            vel0(i,j,k) = -umax(j)
          end if
      end do
    end do! j
   else
    do j=1,ny
      do i=1,nx
          if(terrav(i,j,k)*v(i,j,k).gt.umax(j)) then
            tau_vel(i,j,k)= tauy(j)
            vel0(i,j,k) = umax(j)
          else if(terrav(i,j,k)*v(i,j,k).lt.-umax(j)) then
            tau_vel(i,j,k)= tauy(j)
            vel0(i,j,k) = -umax(j)
          end if
      end do
    end do
   end if
  end do ! k
end if

do j=1,ny
 do i=1,nx
   velmax(i,j) = maxval(v(i,j,:))
 end do
end do

do j=1,ny
 jb=j-YES3D
 do i=1,nx
  k = k_terrav(i,j)
  kc =min(nzm,k+1)
  iadzw = rdz2*rhow(kc)/(adzw(kc)*adz(k)*rho(k))
  tkz=(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
             +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
            /(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
  tkz = tkz*min(terrav(i,j,kc),terrav(i,j,k))
  c = -tkz*iadzw
  b = 1. + dtn*tau_vel(i,j,k) - c
  d = v(i,j,k) + dtn*tau_vel(i,j,k)*vel0(i,j,k) &
               + dtn*rhow(k)/(rho(k)*dz*adz(k))*fluxbv(i,j)
  alpha(i,j,k) = -c/b
  beta(i,j,k) = d/b
 end do
end do
do k=2,nzm-1
 kc=k+1
 kb=k-1
 iadz = rdz2*rhow(k)/(adzw(k)*adz(k)*rho(k))
 iadzw = rdz2*rhow(kc)/(adzw(kc)*adz(k)*rho(k))
 do j=1,ny
   jb=j-YES3D
   do i=1,nx
    if(k.gt.k_terrav(i,j)) then
      tkz=(tk(i,j,kb)*terra(i,j,kb)+tk(i,jb,kb)*terra(i,jb,kb) &
              +tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k)) &
             /(terra(i,j,kb)+terra(i,jb,kb)+terra(i,j,k)+terra(i,jb,k)+1.e-10)
      tkz = tkz*min(terrav(i,j,kb),terrav(i,j,k))
      a = -tkz*iadz
      tkz=(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
              +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
             /(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
      tkz = tkz*min(terrav(i,j,kc),terrav(i,j,k))
      c = -tkz*iadzw
      b = 1. + dtn*tau_vel(i,j,k) - a - c
      d = v(i,j,k) + dtn*tau_vel(i,j,k)*vel0(i,j,k) 
      e = b + a*alpha(i,j,k-1)
      alpha(i,j,k) = -c/e
      beta(i,j,k) = (d-a*beta(i,j,k-1))/e
    end if
   end do
 end do
end do 
k = nzm
kb = k-1
kc = k+1
do j=1,ny
 jb=j-YES3D
 do i=1,nx
  tkz=(tk(i,j,kb)*terra(i,j,kb)+tk(i,jb,kb)*terra(i,jb,kb) &
              +tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k)) &
             /(terra(i,j,kb)+terra(i,jb,kb)+terra(i,j,k)+terra(i,jb,k)+1.e-10)
  tkz = tkz*min(terrav(i,j,kb),terrav(i,j,k))
  iadz = rdz2*rhow(k)/(adzw(k)*adz(k)*rho(k))
  a = -tkz*iadz
  b = 1. + dtn*tau_vel(i,j,k) - a
  d = v(i,j,k) + dtn*tau_vel(i,j,k)*vel0(i,j,k) &
               - dtn*rhow(kc)/(rho(k)*dz*adz(k))*fluxtv(i,j)
  e = b + a*alpha(i,j,k-1)
  v(i,j,k)  = (d-a*beta(i,j,k-1))/e
 end do
end do
do k=nzm-1,1,-1
 kc=k+1
 do j=1,ny
   do i=1,nx
    if(k.ge.k_terrav(i,j)) then
      v(i,j,k) = alpha(i,j,k)*v(i,j,kc) + beta(i,j,k)
    end if
   end do
 end do
end do 

! -----------------------------------------
! Diffusion/damping of w:

! compute inverse damping scale for damping of w
tauz(:) = 0.
if(dodamping.and..not.doregion) then  ! damping already done by nudging for doregion=T
 i=0
 do k=1,nzm
   nu = (zi(k)-zi(1))/(zi(nzm)-zi(1))
   if(nu.gt.nub) then
   !   tauz(k) = tau_max*sin(0.5*pi*(nu-nub)/(1.-nub))**2
      zzz = 100.*((nu-nub)/(1.-nub))**2
      tauz(k) = tau_max*zzz/(1.+zzz)
   end if
 end do
end if

vel0(:,:,:) = 0.
tau_vel(:,:,:) = 0.
if(dodamping_w) then
  do k = 1,nzm
   if(tauz(k).eq.0.) then
     wmax1 = damping_w_cu*dz*adzw(k)/dtn
     do j=1,ny
       wmax = wmax1 !*dx*mu(j)/(dy*ady(j))
       do i=1,nx
          if(terraw(i,j,k)*w(i,j,k).gt.wmax) then
            tau_vel(i,j,k)= tau_vel0
            vel0(i,j,k) = wmax
          else if(terraw(i,j,k)*w(i,j,k).lt.-wmax) then
            tau_vel(i,j,k)= tau_vel0
            vel0(i,j,k) = -wmax
          end if
       end do! i
     end do! j
   end if
  end do ! k
end if

do j=1,ny
 do i=1,nx
  k = k_terra(i,j)
  alpha(i,j,k) = 0.
  beta(i,j,k) = 0.
 end do
end do
do k=2,nzm
 kc=k+1
 kb=k-1
 iadz = rdz2*rho(kb)/(adzw(k)*adz(kb)*rhow(k))
 iadzw = rdz2*rho(k)/(adzw(k)*adz(k)*rhow(k))
 do j=1,ny
   do i=1,nx
    if(k.gt.k_terra(i,j)) then
      tkz=tk(i,j,kb)
      tkz = tkz*min(terraw(i,j,kb),terraw(i,j,k))
      a = -tkz*iadz
      tkz=tk(i,j,k)
      tkz = tkz*min(terraw(i,j,kc),terraw(i,j,k))
      c = -tkz*iadzw
      b = 1. + dtn*(tauz(k)+tau_vel(i,j,k)) - a - c
      d = w(i,j,k) + dtn*tau_vel(i,j,k)*vel0(i,j,k) 
      e = b + a*alpha(i,j,k-1)
      alpha(i,j,k) = -c/e
      beta(i,j,k) = (d-a*beta(i,j,k-1))/e
    end if
   end do
 end do
end do 
do k=nzm,2,-1
 kc=k+1
 do j=1,ny
   do i=1,nx
     if(k.ge.k_terra(i,j)) then
      w(i,j,k) = alpha(i,j,k)*w(i,j,kc) + beta(i,j,k)
    end if
   end do
 end do
end do

!---------------------------------------------------------------------------
! for statistics (fluxes)
if (dostatis) then
 rdz=1./dz
 do k=1,nzm-1
  kc=k+1
  rdz2 = rdz*rdz*rho(k)/adz(k)
  rdz25 = rdz*rdz*rhow(kc)/adzw(kc)
   do j=1,ny
    jb=j-YES3D
    do i=1,nx
     ib=i-1
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
              /(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
     fuz(i,j,kc)=fuz(i,j,kc)-0.5*rdz25*tkz*(u(i,j,kc)-u(i,j,k))*min(terrau(i,j,kc),terrau(i,j,k))
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
              /(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
     fvz(i,j,kc)=fvz(i,j,kc)-0.5*rdz25*tkz*(v(i,j,kc)-v(i,j,k))*min(terrav(i,j,kc),terrav(i,j,k))
   end do
  end do
 end do
 do k=1,nzm-1
   kc=k+1
   do j=1,ny
    jb=j-YES3D
    do i=1,nx
     uwsb(kc)=uwsb(kc)+fuz(i,j,kc)*wgtu(j,kc)*terrau(i,j,kc)
     vwsb(kc)=vwsb(kc)+fvz(i,j,kc)*wgtv(j,kc)*terrav(i,j,kc)
   end do
  end do
 end do
 do j=1,ny
  do i=1,nx
    uwsb(k_terrau(i,j)) = uwsb(k_terrau(i,j)) + fuz(i,j,k_terrau(i,j)) &
                             *wgtu(j,k_terrau(i,j))*terrau(i,j,k_terrau(i,j))
    vwsb(k_terrav(i,j)) = vwsb(k_terrav(i,j)) + fvz(i,j,k_terrav(i,j)) &
                             *wgtv(j,k_terrav(i,j))*terrav(i,j,k_terrav(i,j))
  end do
 end do
 deallocate(fuz,fvz)
end if

if(dodebug) then
  if(masterproc) print*,'after diffuse_mom_z:'
  call fminmax_print('u:',u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nz),1,nx,1,ny,nz)
end if
call t_stopf('diffuse_damping_mom_z')
end subroutine diffuse_damping_mom_z
