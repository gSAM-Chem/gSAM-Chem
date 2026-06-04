subroutine diffuse_mom_z
	
!        momentum SGS diffusion in the vertical

use vars
use sgs
use params, only: docolumn, dowallx, dowally, dosgs, doimplicitdiff, &
                  z0, dodamping, nub, dodebug
use terrain
use check_energy
implicit none

real rdz2,rdz,rdx2u,rdx2v,rdx2w,rdy2u,rdy2v,rdy2w,rdz25

integer i,j,k,ic,ib,jb,jc,kb,kc,kcu
real(8) tkx, tky, tkz, iadzw, iadz, a, b, c, d, e, wnd
real(8) alpha(nx,ny,nzm),beta(nx,ny,nzm)
real nu, tau_max, tauz(nzm)
real, allocatable ::  fuz(:,:,:),fvz(:,:,:)

!-----------------------------------------
if(.not.(dosgs.and.doimplicitdiff)) return

call t_startf('diffuse_mom_z')

if(dodebug) then
  if(masterproc) print*,'before diffuse_mom_z:'
  call fminmax_print('u:',u(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('v:',v(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('w:',w(1:nx,1:ny,1:nz),1,nx,1,ny,nz)
end if

! compute inverse damping scale for damping of w
tauz(:) = 0.
if(dodamping) then
 tau_max=1./10.
 do k=1,nzm
   nu = (zi(k)-zi(1))/(zi(nzm)-zi(1))
   if(nu.gt.nub) then
      tauz(k) = tau_max*sin(0.5*pi*(nu-nub)/(1.-nub))**2
   end if
 end do
end if

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
    jb=j-1
    do i=1,nx
     ib=i-1
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
              /(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
     tkz=max(tkminu(j,k),min(tkmaxu(j,k),tkz))
     fuz(i,j,kc)=-0.5*rdz25*tkz*(u(i,j,kc)-u(i,j,k))*min(terrau(i,j,kc),terrau(i,j,k))
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
              /(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
     tkz=max(tkminv(j,k),min(tkmaxv(j,k),tkz))
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
rdz2 = 0.5*dtn*rdz*rdz 

!----------------------------------------------------
! Diffusion of u:

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
  b = 1. - c
  d = u(i,j,k) +c*(u(i,j,k)-u(i,j,kc)) + dtn*rhow(k)/(dz*adz(k)*rho(k))*fluxbu(i,j)
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
      b = 1. - a - c
      d = u(i,j,k)+a*(u(i,j,k)-u(i,j,kb))+c*(u(i,j,k)-u(i,j,kc))
      e = b + a*alpha(i,j,k-1)
      alpha(i,j,k) = -c/e
      beta(i,j,k) = (d-a*beta(i,j,k-1))/e
    end if
   end do
 end do
end do 
k = nzm
kb = k-1
do j=1,ny
 do i=1,nx
  ib=i-1
  tkz = (tk(i,j,kb)*terra(i,j,kb)+tk(ib,j,kb)*terra(ib,j,kb) &
          +tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k)) &
         /(terra(i,j,kb)+terra(ib,j,kb)+terra(i,j,k)+terra(ib,j,k)+1.e-10)
  tkz = tkz*min(terrau(i,j,kb),terrau(i,j,k))
  iadz = rdz2*rhow(k)/(adzw(k)*adz(k)*rho(k))
  a = -tkz*iadz
  b = 1. - a
  d = u(i,j,k)+a*(u(i,j,k)-u(i,j,kb)) - dtn*rhow(kc)/(rho(k)*dz*adz(k))*fluxtu(i,j)
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
! Diffusion of v:
  
do j=1,ny
 jb=j-1
 do i=1,nx
  k = k_terrav(i,j)
  kc =min(nzm,k+1)
  iadzw = rdz2*rhow(kc)/(adzw(kc)*adz(k)*rho(k))
  tkz=(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
             +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
            /(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
  tkz = tkz*min(terrav(i,j,kc),terrav(i,j,k))
  c = -tkz*iadzw
  b = 1. - c
  d = v(i,j,k)+c*(v(i,j,k)-v(i,j,kc)) + dtn*rhow(k)/(rho(k)*dz*adz(k))*fluxbv(i,j)
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
   jb=j-1
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
      b = 1. - a - c
      d = v(i,j,k)+a*(v(i,j,k)-v(i,j,kb))+c*(v(i,j,k)-v(i,j,kc))
      e = b + a*alpha(i,j,k-1)
      alpha(i,j,k) = -c/e
      beta(i,j,k) = (d-a*beta(i,j,k-1))/e
    end if
   end do
 end do
end do 
k = nzm
kb = k-1
do j=1,ny
 jb=j-1
 do i=1,nx
  tkz=(tk(i,j,kb)*terra(i,j,kb)+tk(i,jb,kb)*terra(i,jb,kb) &
              +tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k)) &
             /(terra(i,j,kb)+terra(i,jb,kb)+terra(i,j,k)+terra(i,jb,k)+1.e-10)
  tkz = tkz*min(terrav(i,j,kb),terrav(i,j,k))
  iadz = rdz2*rhow(k)/(adzw(k)*adz(k)*rho(k))
  a = -tkz*iadz
  b = 1. - a
  d = v(i,j,k) + a*(v(i,j,k)-v(i,j,kb)) - dtn*rhow(kc)/(rho(k)*dz*adz(k))*fluxtv(i,j)
  e = b + a*alpha(i,j,k-1)
  v(i,j,k)  = (d-a*beta(i,j,k-1))/e
 end do
end do
do k=nzm-1,1,-1
 kc=k+1
 do j=1,ny
   do i=1,nx
    if(k_terrav(i,j).ge.nzm) cycle
    if(k.ge.k_terrav(i,j)) then
      v(i,j,k) = alpha(i,j,k)*v(i,j,kc) + beta(i,j,k)
    end if
   end do
 end do
end do 

! -----------------------------------------
! Diffusion of w:

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
      b = 1. + dtn*tauz(k) - a - c
      d = w(i,j,k) +a*(w(i,j,k)-w(i,j,kb))+c*(w(i,j,k)-w(i,j,kc))
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
    jb=j-1
    do i=1,nx
     ib=i-1
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(ib,j,k)*terra(ib,j,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(ib,j,kc)*terra(ib,j,kc)) &
              /(terra(i,j,k)+terra(ib,j,k)+terra(i,j,kc)+terra(ib,j,kc)+1.e-10)
     tkz=max(tkminu(j,k),min(tkmaxu(j,k),tkz))
     fuz(i,j,kc)=fuz(i,j,kc)-0.5*rdz25*tkz*(u(i,j,kc)-u(i,j,k))*min(terrau(i,j,kc),terrau(i,j,k))
     tkz=(tk(i,j,k)*terra(i,j,k)+tk(i,jb,k)*terra(i,jb,k) &
               +tk(i,j,kc)*terra(i,j,kc)+tk(i,jb,kc)*terra(i,jb,kc)) &
              /(terra(i,j,k)+terra(i,jb,k)+terra(i,j,kc)+terra(i,jb,kc)+1.e-10)
     tkz=max(tkminv(j,k),min(tkmaxv(j,k),tkz))
     fvz(i,j,kc)=fvz(i,j,kc)-0.5*rdz25*tkz*(v(i,j,kc)-v(i,j,k))*min(terrav(i,j,kc),terrav(i,j,k))
   end do
  end do
 end do
 do k=1,nzm-1
   kc=k+1
   do j=1,ny
    jb=j-1
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
call t_stopf('diffuse_mom_z')
end subroutine diffuse_mom_z
