#include "fppmacros"

subroutine setdata()
	
use vars
use params
use terrain, only: terrau, terrav, terraw
!use micro_params
use microphysics, only: micro_init, micro_proc, micro_set
use chemistry, only: chem_init
use sgs, only: sgs_init, sgs_proc, dodns, RHO_DNS, PRES_DNS
use simple_ocean, only: set_sst

implicit none
	
integer ndmax,n,i,j,k,iz,it,jt,m
real presr(nz)	
parameter (ndmax = 1000)
real zz(ndmax),tt(ndmax),qq(ndmax),uu(ndmax),vv(ndmax) 
real zz1(ndmax),tt1(ndmax),qq1(ndmax),uu1(ndmax),vv1(ndmax) 
real tmp,rrr1,rrr2, pres1, pz(ndmax),ta(ndmax)
real pz1(ndmax)
real ratio_t1,ratio_t2,ratio_p1,ratio_p2
real tpert0(ndmax), qpert0(ndmax)
real latit,long
logical:: zgrid, tabssound=.false.
integer status
integer nx1,ny1,nz1
real  coef
real(4), allocatable :: presin(:), zin(:) 
real, external :: qsatw,qsati
real(8), allocatable :: xx(:), yy(:), xx_u(:), yy_v(:)
real(8), external :: utTimeSeconds
integer itime(6)

call t_startf ('setdata')

!-------------------------------------------------------------
!	read subensemble perturbation file first:

if(doensemble) then
	
open(76,file=trim(rundatadir)//'/tqpert',status='old',form='formatted')  
read(76,*)
  do j=0,nensemble
    read(76,*) i,n
    do i=1,n
      read(76,end=766,fmt=*) pz(i),tpert0(i),qpert0(i)
      tpert0(i)=tpert0(i)*(1000./pz(i))**(rgas/cp)
    end do
  end do
  close(76)
  if(masterproc) then
    print*,'Subensemble run. nensemble=',nensemble
    print*,'tpert:',(tpert0(i),i=1,n)
    print*,'qpert:',(qpert0(i),i=1,n)
  end if
  goto 767
766  print*,'Error: nensemble is too large.'  
  call task_abort()
767  continue
else
  do i=1,ndmax
    tpert0(i)=0.
    qpert0(i)=0.
  end do
end if

w = 0.
qcl = 0.
qci = 0.
qpl = 0.
qpi = 0.
qpg = 0.
p = 0.
dudt = 0.
dvdt = 0.
dwdt = 0.	   
	
!**************************************************************
! compute pressures and densities for dodns = .true.
if(dodns) then
 rho(:) = RHO_DNS
 rhow(:) = RHO_DNS
 pres0 = PRES_DNS
 presi(1) = pres0
 do k=2,nz
  presi(k) = presi(k-1) - ggr*RHO_DNS*adz(k-1)*dz*0.01
 end do
 pres(1) = presi(1) - ggr*RHO_DNS*0.5*dz*0.01
 do k=2,nzm
  pres(k) = pres(k-1) - ggr*RHO_DNS*adzw(k)*dz*0.01
 end do
end if

!**************************************************************
if(.not.readinit) then

!	Read Initial Sounding

 if(doscamiopdata) then

   !bloss: doensemble not implemented in conjunction with doscamiopdata yet
   if(doensemble) then
      if(masterproc)print *,'doensemble does not work with doscamiopdata yet'
      call task_abort()
   end if

   !bloss: doradforcing not implemented in conjunction with doscamiopdata yet
   if(doradforcing) then
      if(masterproc)print *,'doradforcing does not work with doscamiopdata yet'
      call task_abort()
   end if
   if(masterproc) print*,'Reading soundings from SCAM input file.'

   !bloss: read sounding/forcing data from SCAM input file.
   call readiopdata(status)

   isInitialized_scamiopdata = .true.

   !bloss: Interpolate sounding data to initial time.
   !       It has already been interpolated onto the model's z grid
   !         within readiopdata.
   do i = 1,nsnd-1
      if((day.ge.daysnd(i)).and.(day.lt.daysnd(i+1))) then
         coef = (day-daysnd(i)) / (daysnd(i+1)-daysnd(i))
         pres0 = (1-coef)*pres0ls(i) + coef*pres0ls(i+1) !surface pressure [mb]
         do k = 1,nzsnd
            pz(k) = (1-coef)*psnd(k,i) + coef*psnd(k,i+1) !absolute temp [K]
            tt(k) = (1-coef)*tsnd(k,i) + coef*tsnd(k,i+1) !absolute temp [K]
            qq(k) = (1-coef)*qsnd(k,i) + coef*qsnd(k,i+1) !tot water [g/kg]
            uu(k) = (1-coef)*usnd(k,i) + coef*usnd(k,i+1) !u wind [m/s]
            vv(k) = (1-coef)*vsnd(k,i) + coef*vsnd(k,i+1) !v wind [m/s]
            ta(k)=  tt(k)*(pz(k)/1000.)**(rgas/cp)
         end do
         exit !break out of do i=1,nsnd-1
      elseif(i.eq.nsnd-1) then
         if(masterproc) print*,'Error: day is beyond the sounding time range'
         call task_abort()
      end if
   end do

   zgrid = .false. ! SCAM input based on pressure.

   n = nzsnd

 else if(dosimplesnd) then

  open(77,file='./CASES/'//trim(case)//'/snd',status='old',form='formatted')
  if(masterproc) print*,'dosimplesnd=',dosimplesnd
  if(masterproc) print*,'Reading simple soundings from '//'./CASES/'//trim(case)//'/snd'
  read(77,*)
  read(77,*) n, pres0
  if(n.ne.nzm) then
   if(masterproc) print*,'for dosimplesnd=T, number of levels in snd should be the same as nzm!'
   call task_abort()
  end if
  presi(1) = pres0
  do k=1,nzm
   read(77,*) zz(k),pres(k),tt(k),qq(k),uu(k),vv(k)
   if(abs(zz(k)-z(k)).gt.1.) then
     if(masterproc) print*,'setdata:dosimplesnd=T: wrong snd file! heights and z() should match!'
     call task_abort()
   end if
   if(docap_snd_cu) call snd_cap(uu(k),vv(k))
   presi(k+1) = 2.*pres(k)-presi(k)
  end do	      
  close (77)
  if(tt(1).lt.0.) tabssound = .true.
  zgrid = .true.

 else

 open(77,file='./CASES/'//trim(case)//'/snd',status='old',form='formatted')
 if(masterproc) print*,'Reading soundings from '//'./CASES/'//trim(case)//'/snd'
 read(77,*)

 do while(.true.)

   read(77,err=55,end=55,fmt=*) rrr1, n, pres0
   do i=1,n
       read(77,*) zz(i),pz(i),tt(i),qq(i),uu(i),vv(i)
       if(docap_snd_cu) call snd_cap(uu(i),vv(i))
   end do	      
   read(77,err=55,end=55,fmt=*) rrr2, n, pres1
   do i=1,n
       read(77,*) zz1(i),pz1(i),tt1(i),qq1(i),uu1(i),vv1(i)
       if(docap_snd_cu) call snd_cap(uu1(i),vv1(i))
   end do	      
 
   if(day.ge.rrr1.and.day.le.rrr2) then
     if(tt(1).lt.0.) tabssound = .true.
     if(zz(2).gt.zz(1)) then
       zgrid = .true.
       do i=1,n
         zz(i)=zz(i)+(zz1(i)-zz(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         tt(i)=tt(i)+(tt1(i)-tt(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         qq(i)=qq(i)+(qq1(i)-qq(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         uu(i)=uu(i)+(uu1(i)-uu(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         vv(i)=vv(i)+(vv1(i)-vv(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         tt(i)=tt(i)+tpert0(i)
         qq(i)=qq(i)+qpert0(i) 
       end do
     else if(pz(2).lt.pz(1)) then
       zgrid = .false.
       do i=1,n
         pz(i)=pz(i)+(pz1(i)-pz(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         tt(i)=tt(i)+(tt1(i)-tt(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         qq(i)=qq(i)+(qq1(i)-qq(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         uu(i)=uu(i)+(uu1(i)-uu(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         vv(i)=vv(i)+(vv1(i)-vv(i))/(rrr2-rrr1+1.e-5)*(day-rrr1)
         tt(i)=tt(i)+tpert0(i)
         qq(i)=qq(i)+qpert0(i) 
         if(tabssound) then
          ta(i) = -tt(i)
         else
          ta(i)=tt(i)*(pz(i)/1000.)**(rgas/cp)
         end if
       end do
     else  
       if(masterproc) print*,'vertical grid is undefined...'
     end if
     pres0=pres0+(pres1-pres0)/(rrr2-rrr1+1.e-5)*(day-rrr1)      
     goto 56
   endif
   do i=1,n+1
     backspace(77)
!     read(77,*)      ! a bug in gfortran compiler
!     backspace(77) ! these two lines were addedf because of 
   end do	      
 
 end do
 
 55 continue
 if(masterproc) then
   print*,'Error: day is beyond the sounding time range'
   print*,day,rrr1,rrr2
 end if
 call task_abort()
 56 continue	
 
 close (77)
 
 end if ! if(doscamiopdata)
 
 if(masterproc) then
   print *	
   print *,'surface pressure from snd file: ',pres0
 endif  
 
 ! compute heights from pressure:
 
 if(.not.zgrid) then
   zz(1) = rgas/ggr*ta(1)*log(pres0/pz(1))
   do i=2,n
    zz(i)=zz(i-1)+0.5*rgas/ggr*(ta(i)+ta(i-1))*log(pz(i-1)/pz(i))
   end do
 end if  	
 !-----------------------------------------------------------
 !       Interpolate sounding into vertical grid:
 
 if(dodns) pres0 = PRES_DNS
 if(.not.(dodns.or.dosimplesnd)) then
  presr(1)=(pres0/1000.)**(rgas/cp)
  presi(1)=pres0
 end if
 do k= 1,nzm
  do iz = 2,n
   if(z(k).le.zz(iz)) then
     t0(k)=tt(iz-1)+(tt(iz)-tt(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1))	
     q0(k)=qq(iz-1)+(qq(iz)-qq(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1))
     u0(k)=uu(iz-1)+(uu(iz)-uu(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1))  
     v0(k)=vv(iz-1)+(vv(iz)-vv(iz-1))/(zz(iz)-zz(iz-1))*(z(k)-zz(iz-1)) 
     goto 12
   endif
  end do
 
 !  Utilize 1976 standard atmosphere for points above sounding:
 
 
  if(k.gt.1.) then
   call atmosphere(z(k-1)/1000.,ratio_p1,rrr1,ratio_t1)
   call atmosphere(z(k)/1000.,ratio_p2,rrr1,ratio_t2)
 
   tabs0(k)=ratio_t2/ratio_t1*tabs0(k-1)
   if(.not.(dodns.or.dosimplesnd)) then
    presi(k+1)=presi(k)*exp(-ggr/rgas/tabs0(k)*(zi(k+1)-zi(k)))
    pres(k) = 0.5*(presi(k)+presi(k+1))
   end if
   prespot(k)=(1000./pres(k))**(rgas/cp)
   q0(k) = q0(k-1)*exp(-(z(k)-z(k-1))/3000.) ! always decrease q0 with height
   u0(k)=u0(k-1)
   v0(k)=v0(k-1)
  else
   if(masterproc) print*,'in setdata, starndard atmosphere is called when no valid sounding is given'
   call task_abort()
  end if
  goto 13
 12 continue
  if(tabssound) then
    tabs0(k) = -t0(k)
    if(.not.(dodns.or.dosimplesnd)) then
     presi(k+1)=presi(k)*exp(-ggr/rgas/tabs0(k)*(zi(k+1)-zi(k)))
     pres(k) = exp(log(presi(k))+log(presi(k+1)/presi(k))* &
                              (z(k)-zi(k))/(zi(k+1)-zi(k)))
    end if
    prespot(k)=(1000./pres(k))**(rgas/cp)
  else
    tv0(k)=t0(k)*(1.+epsv*q0(k)*1.e-3)
!    tv0(k)=t0(k)
    if(.not.(dodns.or.dosimplesnd)) then
     presr(k+1)=presr(k)-ggr/cp/tv0(k)*(zi(k+1)-zi(k))
     presi(k+1)=1000.*presr(k+1)**(cp/rgas)
     pres(k) = exp(log(presi(k))+log(presi(k+1)/presi(k))* &
                              (z(k)-zi(k))/(zi(k+1)-zi(k)))
    end if
    prespot(k)=(1000./pres(k))**(rgas/cp)
    tabs0(k)=t0(k)/prespot(k)
  end if
  if(q0(k).lt.0.) then ! convert from relative humdity to mixing ratio
   q0(k) = -q0(k)/100.*qsatw(tabs0(k),pres(k))
  else
   q0(k)=q0(k)*1.e-3
  end if
 13 continue
  ug0(k)=u0(k)
  vg0(k)=v0(k)
  t0(k) = tabs0(k)+gamaz(k) 
!  t0(k) = tabs0(k)+ggr*z(k)/(cp*(1.-cpvf*qv0(k))+cpv*cpvf*qv0(k)) 
  qv0(k) = q0(k)
  qn0(k) = 0.
  qp0(k) = 0.
 end do
 
  do k=1,nzm
     u(:,:,k)= u0(k)
     v(:,:,k)= v0(k)
     t(:,:,k)= t0(k)
     tabs(:,:,k) = tabs0(k)
     qcl(:,:,k)=0.
     qci(:,:,k)=0.
     qpl(:,:,k)=0.
     qpi(:,:,k)=0.
     qpg(:,:,k)=0.
  end do
  w(:,:,:)= 0.

else
 
  open(11,file=initfile,form='unformatted',BUFFEREDYES ACTION='READ')
  read(11) nx1,ny1,nz1
  if(masterproc) print*,'reading initfile=',trim(initfile)
  if(masterproc) print*,'Is vertical velocity assume pressure velocity? w3D_pressure=',w3D_pressure
  if(masterproc) print*,'data size:',nx1,ny1,nz1
  allocate(zin(nz1),presin(nz1))
  read(11) zin(1:nz1)
  read(11) presin(1:nz1)
  if(masterproc) print*,'zin=',zin
  if(masterproc) print*,'presin=',presin
  if(.not.dodns) then  
   do k=1,nzm
    do m=1,nz1
     if(z(k).le.zin(1).or.(z(k).ge.zin(m).and.z(k).lt.zin(m+1))) then
       pres(k) = exp(log(presin(m))+(log(presin(m+1))-log(presin(m)))/(zin(m+1)-zin(m))*(z(k)-zin(m)))
       exit
     end if
    end do
   end do
   pres0 = exp(log(presin(1))+(log(presin(2))-log(presin(1)))/(zin(2)-zin(1))*(0.-zin(1)))
  end if ! .not.dodns
  deallocate(zin,presin)
  allocate(xx(nx_gl),yy(ny_gl),xx_u(nx_gl+1),yy_v(ny_gl+1))
  if(read_meters) then
    xx(:) = x_gl(:)
    yy(:) = y_gl(:)
    xx_u(:) = xu_gl(:)
    yy_v(:) = yv_gl(:)
  else
    xx(:) = lon_gl(:)
    yy(:) = lat_gl(:)
    xx_u(:) = lonu_gl(:)
    yy_v(:) = latv_gl(:)
  end if
  if(masterproc) print*,'pres0=',pres0
  if(masterproc) print*,'pres=',pres
  if(masterproc) print*,'reading u:'
  call read_field3D (11,u(1:nx+1,1:ny,1:nzm),nx1,ny1,nz1,xx_u,yy,z,nx+1,ny,nzm)
  call fminmax_print('u:',u(1:nx+1,1:ny,1:nzm),1,nx+1,1,ny,nzm)
  if(masterproc) print*,'reading v:'
  call read_field3D (11,v(1:nx,1:ny+YES3D,1:nzm),nx1,ny1,nz1,xx,yy_v,z,nx,ny+YES3D,nzm)
  call fminmax_print('v:',v(1:nx,1:ny+YES3D,1:nzm),1,nx,1,ny+YES3D,nzm)
  if(masterproc) print*,'reading omega:'
  call read_field3D (11,w(1:nx,1:ny,2:nzm),nx1,ny1,nz1,xx,yy,zi(2:nzm),nx,ny,nzm-1)
  w(:,:,1) = 0.
  w(:,:,nz) = 0.
  call fminmax_print('omega:',w(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  if(masterproc) print*,'reading tabs:'
  call read_field3D (11,tabs(1:nx,1:ny,1:nzm),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call fminmax_print('tabs:',tabs(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  if(masterproc) print*,'reading q:'
  call read_field3D (11,qv(1:nx,1:ny,1:nzm),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call fminmax_print('qv:',qv(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  if(masterproc) print*,'reading qcl:'
  call read_field3D (11,qcl(1:nx,1:ny,1:nzm),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call fminmax_print('qcl:',qcl(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  if(masterproc) print*,'reading qci:'
  call read_field3D (11,qci(1:nx,1:ny,1:nzm),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call fminmax_print('qci:',qci(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  if(masterproc) print*,'reading qpl:'
  call read_field3D (11,qpl(1:nx,1:ny,1:nzm),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call fminmax_print('qpl:',qpl(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  if(masterproc) print*,'reading qpi:'
  call read_field3D (11,qpi(1:nx,1:ny,1:nzm),nx1,ny1,nz1,xx,yy,z,nx,ny,nzm)
  call fminmax_print('qpi:',qpi(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  deallocate(xx, yy, xx_u, yy_v)
  call micro_set(tabs,t(1:nx,1:ny,1:nzm),qv,qcl,qci,qpl,qpi,qpg) ! make consistent fields
  if(masterproc) print*,'after micro_set:'
  call fminmax_print('t:',t(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('qv:',qv(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('qcl:',qcl(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('qci:',qci(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('qpl:',qpl(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('qpi:',qpi(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call fminmax_print('qpg:',qpg(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  call diagnose()
  call fminmax_print('tabs:',tabs(1:nx,1:ny,1:nzm),1,nx,1,ny,nzm)
  if(.not.dodns) then
   presr(1)=(pres0/1000.)**(rgas/cp)
   presi(1)=pres0
   do k=1,nzm
     prespot(k)=(1000./pres(k))**(rgas/cp)
     t0(k)=tabs0(k)*prespot(k)
     presr(k+1)=presr(k)-ggr/cp/t0(k)*(zi(k+1)-zi(k))
     presi(k+1)=1000.*presr(k+1)**(cp/rgas)
     pres(k) = exp(log(presi(k))+log(presi(k+1)/presi(k))* &
                               (z(k)-zi(k))/(zi(k+1)-zi(k)))
   end do
  end if
  do k=1,nzm
    prespot(k)=(1000./pres(k))**(rgas/cp)
    ug0(k)=0.
    vg0(k)=0.
  end do

end if 

prespoti(:) = (1000./presi(:))**(rgas/cp)

! recompute pressure levels (for consistancy):

!	call pressz()
        
!-------------------------------------------------------------
	
! add domain translation velocities:
u(1:nx,1:ny,:) = u(1:nx,1:ny,:) - ug
v(1:nx,1:ny,:) = v(1:nx,1:ny,:) - vg

do k=1,nzm
  qn0(k) = 0.
  qp0(k) = 0.
  u0(k) = u0(k) - ug
  v0(k) = v0(k) - vg
  ug0(k) = ug0(k) - ug
  vg0(k) = vg0(k) - vg
  if(.not.dodns) rho(k) = (presi(k)-presi(k+1))/(zi(k+1)-zi(k))/ggr*100.
end do

do k=1,nzm
 pp(:,:,k) = pres(k)
end do

if(.not.dodns) then
 if(dolatlon) then
  ! requred by GMG pressure solver
  do k=2,nzm
    rhow(k) = (rho(k-1)*adz(k)+rho(k)*adz(k-1))/(adz(k)+adz(k-1))
  end do
 else
  do k=2,nzm
    rhow(k) =  (pres(k-1)-pres(k))/(z(k)-z(k-1))/ggr*100.
  end do
 end if
 rhow(1) = 2*rho(1)-rhow(2)
 rhow(nz)= 2*rho(nzm)-rhow(nzm)
end if

if(readinit.and.w3D_pressure) then
  ! convert from Ps/s to m/s 
  do k=nzm,2,-1
    w(:,:,k) = -(adz(k)*w(:,:,k)+adz(k-1)*w(:,:,k-1)) &
               /(adz(k)+adz(k-1))/(rhow(k)*ggr)
  end do
end if
w(:,:,1) = 0.
w(:,:,nz) = 0.
p(:,:,:,:)=0.
eis(:,:) = 0.


if(docloud.or.dosmoke) call micro_init()  !initialize microphysics

if (dochem) call chem_init()

if(doregion) then
      allocate(uobs(nx+1,ny,nzm,2))
      allocate(vobs(nx,ny+YES3D,nzm,2))
      allocate(wobs(nx,ny,nz,2))
      allocate(tobs(nx,ny,nzm,2))
      allocate(qobs(nx,ny,nzm,2))
      if(masterproc) print*,'repairing divergent wind field for doregion=.true.:'
      u(1:nx+1,1:ny,1:nzm) = u(1:nx+1,1:ny,1:nzm) * terrau(1:nx+1,1:ny,1:nzm)
      v(1:nx,1:ny+YES3D,1:nzm) = v(1:nx,1:ny+YES3D,1:nzm) * terrav(1:nx,1:ny+YES3D,1:nzm)
      w(1:nx,1:ny,1:nzm) = w(1:nx,1:ny,1:nzm) * terraw(1:nx,1:ny,1:nzm)
      if(masterproc) print*,'divergence before repair:'
      call check_div(u(1:nx+1,1:ny,1:nzm),v(1:nx,1:ny+YES3D,1:nzm),w(1:nx,1:ny,1:nz))
      call nondivergent(u(1:nx+1,1:ny,1:nzm),v(1:nx,1:ny+YES3D,1:nzm),w(1:nx,1:ny,1:nz))
      uobs(1:nx+1,1:ny,1:nzm,1) = u(1:nx+1,1:ny,1:nzm)
      vobs(1:nx,1:ny+YES3D,1:nzm,1) = v(1:nx,1:ny+YES3D,1:nzm)
      wobs(1:nx,1:ny,1:nz,1) = w(1:nx,1:ny,1:nz)
      tobs(1:nx,1:ny,1:nzm,1) = tabs(1:nx,1:ny,1:nzm)
      qobs(1:nx,1:ny,1:nzm,1) = qv(1:nx,1:ny,1:nzm)
      uobs(1:nx+1,1:ny,1:nzm,2) = u(1:nx+1,1:ny,1:nzm)
      vobs(1:nx,1:ny+YES3D,1:nzm,2) = v(1:nx,1:ny+YES3D,1:nzm)
      wobs(1:nx,1:ny,1:nz,2) = w(1:nx,1:ny,1:nz)
      tobs(1:nx,1:ny,1:nzm,2) = tabs(1:nx,1:ny,1:nzm)
      qobs(1:nx,1:ny,1:nzm,2) = qv(1:nx,1:ny,1:nzm)
      call task_rank_to_index(rank,it,jt)
      if(it.eq.0) u_w(1:ny,1:nzm) = u(1,1:ny,1:nzm)
      if(it+nx.eq.nx_gl) u_e(1:ny,1:nzm) = u(nx+1,1:ny,1:nzm)
      if(jt.eq.0) v_s(1:nx,1:nzm) = v(1:nx,1,1:nzm)
      if(jt+ny.eq.ny_gl) v_n(1:nx,1:nzm) = v(1:nx,ny+YES3D,1:nzm)
      call boundaries(1)
      if(masterproc) print*,'divergence after repair:'
      call check_div(uobs,vobs,wobs)
end if

call date_from_dayofyear(day,year0,datechar)
if(dodatefilename) then
 date_pr = datechar(1:4)//"-"//datechar(5:6)//"-"//datechar(7:8)//"-"//datechar(9:10)&
          //"-"//datechar(11:12)//"-"//datechar(13:14)
 read(datechar,'(i4,i2,i2,i2,i2,i2)') itime(:)
 if(itime(1).ge.1900) timeUTsec = &
                  utTimeSeconds(itime(1),itime(2),itime(3),itime(4),itime(5),itime(6))
else
  date_pr = ''
end if

!print*,'>>>',(pres(k),',',k=1,nzm)

if(dofixdynamics) call setfixdynamics()

call setperturb()

call diagnose()
bet(:) = ggr/tabs0(:)

if(dosgs) call sgs_init()

!initialize surface:

if(.not.dosfcforcing.and.(OCEAN.or.ISLAND).and..not.readsst) call set_sst()

! Initialize BUILDINGS
call buildings_init() 


call print_profiles

if(masterproc) then
 print*,'k200=',k200
 print*,'k500=',k500
 print*,'k700=',k700
 print*,'k850=',k850
 print*,'k950=',k950
end if

! exchange boundaries:
call boundaries(1)
call boundaries(4)

if(nstep.eq.0) call write_fields3D()
!if(nstep.eq.0) call write_fields2D()
!if(nstep.eq.0) call write_fields2DZ()

if(doterrain.and.nstep.eq.0) then
  call task_barrier()
  call write_terrain()
  call task_barrier()
  if(masterproc) print*,'saved...'
end if

call t_stopf ('setdata')


end


subroutine print_profiles
use vars
use params
implicit none
real coef
integer k
real, external :: qsatw,qsati


if(masterproc) then
 print *,'Initial Sounding:'
 print *, ' k      z    rho     rhoi    s      h     h*      t      u      v   adz  bet'
 do k=nzm,1,-1
   write(6,'(i4,1x,g11.4,2f7.4,f7.2,f7.2,4f7.2,6g11.4)') k,z(k),rho(k),rhow(k),tabs0(k)+ggr*z(k)/cp, &
          t0(k)+lcond/cp*qv0(k), t0(k)+lcond/cp*qsatw(tabs0(k),pres(k)), &
                t0(k),u0(k)+ug,v0(k)+vg, adz(k), bet(k)
 end do
 print *, ' k      z    rho     rhoi    s      h     h*      qt      u      v   adz  bet'

 print *
 print *,'  k      z      dz     pres   presi   Tabs     tp     qt      Qn     REL   RELi'
 coef=1.
 if(dosmoke) coef = 0.
 do k = nzm,1,-1
  write(6,'(i4,1x,6f8.2,2f8.4,2f8.2)')   k,z(k),zi(k+1)-zi(k),pres(k),presi(k),tabs0(k), &
     tabs0(k)*prespot(k),q0(k)*1.e3, qn0(k)*1.e3, &
     coef*100.*qv0(k)/qsatw(tabs0(k),pres(k)), coef*100.*qv0(k)/qsati(tabs0(k),pres(k))
 end do
 print *,'  k      z      dz     pres   presi   Tabs     tp     qt      Qn     REL   RELi'
endif


end

subroutine snd_cap(u,v)
 use params, only: cap_snd_cu,ug,vg
 use grid, only: dx, dy, dt
 implicit none
 real, intent(inout) ::  u,v
 real cu, coef
 integer k
 cu = sqrt(((u-ug)/dx)**2+((v-vg)/dy)**2)*dt
 if(cu.gt.cap_snd_cu) then
    coef = cap_snd_cu/cu
    u = ug+(u-ug)*coef
    v = vg+(v-vg)*coef
 end if
end


!subroutine wind_cap(u,v)
! use params, only: cap_wind
! use grid, only: nx,ny,nzm
! implicit none
! real, intent(inout) ::  u(:,:,:),v(:,:,:)
! real vel, coef
! integer i,j,k
! do k=1,nzm
!  do j=1,ny
!   do i=1,nx
!     vel = sqrt((u(i,j,k)+u(i+1,j,k))**2+(v(i,j,k)**2)
!     if(vel.gt.cap_wind) then
!       coef = cap_wind/vel
!       u(i,j,k) = u(i,j,k)*coef
!       v(i,j,k) = v(i,j,k)*coef
!     end if
!   end do
!  end do
! end do
!end


