
subroutine shear_prod3D(def2)
	
use vars
use terrain, only: k_terra
implicit none
	
real def2(nx,ny,nzm)
	
real idx,rdx,rdx1,rdx2
real idy,rdy,rdy1,rdy2
real rdz,rdz1,rdz2
integer i,j,k,ib,ic,jb,jc,kb,kc

idx=1./dx 
idy=1./dy 

do k=2,nzm-1  

 kb=k-1
 kc=k+1
 rdz = 1./(dz*adz(k))
 rdz1 = 1./(dz*adzw(kc))
 rdz2 = 1./(dz*adzw(k))

 do j=1,ny
   jb=j-YES3D
   jc=j+YES3D
   rdx=idx*imu(j) ! take into account grid anisotropy
   rdx1=idx*imuv(jc) ! take into account grid anisotropy
   rdx2=idx*imuv(j) ! take into account grid anisotropy
   rdy=1./(dy*ady(j))
   rdy1=1./(dy*adyv(jc))
   rdy2=1./(dy*adyv(j))
   do i=1,nx
     ib=i-1
     ic=i+1
      def2(i,j,k)=2.* ( &
          ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
          ( (v(i,jc,k)-v(i,j,k))*rdy)**2+ &
          ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
        + 0.25 * ( &
          ( (u(ic,jc,k)-u(ic,j ,k))*rdy1+(v(ic,jc,k)-v(i ,jc,k))*rdx1 )**2 +  &
          ( (u(i ,jc,k)-u(i ,j ,k))*rdy1+(v(i ,jc,k)-v(ib,jc,k))*rdx1 )**2 +  &
          ( (u(ic,j ,k)-u(ic,jb,k))*rdy2+(v(ic,j ,k)-v(i ,j ,k))*rdx2 )**2 +  &
          ( (u(i ,j ,k)-u(i ,jb,k))*rdy2+(v(i ,j ,k)-v(ib,j ,k))*rdx2 )**2 )   
      def2(i,j,k)=def2(i,j,k) &
        + 0.25 * ( &
          ( (u(ic,j,kc)-u(ic,j, k))*rdz1+(w(ic,j,kc)-w(i ,j,kc))*rdx )**2 + &
          ( (u(i ,j,kc)-u(i ,j, k))*rdz1+(w(i ,j,kc)-w(ib,j,kc))*rdx )**2 + &
          ( (u(ic,j,k )-u(ic,j,kb))*rdz2+(w(ic,j,k )-w(i ,j,k ))*rdx )**2 + &
          ( (u(i ,j,k )-u(i ,j,kb))*rdz2+(w(i ,j,k )-w(ib,j,k ))*rdx )**2 )
      def2(i,j,k)=def2(i,j,k) &
        + 0.25 * ( &
          ( (v(i,jc,kc)-v(i,jc, k))*rdz1+(w(i,jc,kc)-w(i,j ,kc))*rdy1 )**2 + &
          ( (v(i,j ,kc)-v(i,j , k))*rdz1+(w(i,j ,kc)-w(i,jb,kc))*rdy1 )**2 + &
          ( (v(i,jc,k )-v(i,jc,kb))*rdz2+(w(i,jc,k )-w(i,j ,k ))*rdy2 )**2 + &
          ( (v(i,j ,k )-v(i,j ,kb))*rdz2+(w(i,j ,k )-w(i,jb,k ))*rdy2 )**2 )

    end do
 end do
end do ! k


do j=1,ny
  jb=j-YES3D
  jc=j+YES3D
  do i=1,nx
     k=k_terra(i,j)
     kc=k+1
     rdz = 1./(dz*adz(k))
     rdz1 = 1./(dz*adzw(kc))
     rdx=idx*imu(j) ! take into account grid anisotropy
     rdx1=idx*imuv(jc) ! take into account grid anisotropy
     rdx2=idx*imuv(j) ! take into account grid anisotropy
     rdy=1./(dy*ady(j))
     rdy1=1./(dy*adyv(jc))
     rdy2=1./(dy*adyv(j))
     ib=i-1
     ic=i+1
      def2(i,j,k)=2.* ( &
          ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
          ( (v(i,jc,k)-v(i,j,k))*rdy)**2+ &
          ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
        + 0.25 * ( &
          ( (u(ic,jc,k)-u(ic,j ,k))*rdy1+(v(ic,jc,k)-v(i ,jc,k))*rdx1 )**2 +  &
          ( (u(i ,jc,k)-u(i ,j ,k))*rdy1+(v(i ,jc,k)-v(ib,jc,k))*rdx1 )**2 +  &
          ( (u(ic,j ,k)-u(ic,jb,k))*rdy2+(v(ic,j ,k)-v(i ,j ,k))*rdx2 )**2 +  &
          ( (u(i ,j ,k)-u(i ,jb,k))*rdy2+(v(i ,j ,k)-v(ib,j ,k))*rdx2 )**2 )   &
         + 0.5 * ( &
          ( (u(ic,j,kc)-u(ic,j, k))*rdz1+(w(ic,j,kc)-w(i ,j,kc))*rdx )**2 + &
          ( (u(i ,j,kc)-u(i ,j, k))*rdz1+(w(i ,j,kc)-w(ib,j,kc))*rdx )**2 ) &
         + 0.5 * ( &
          ( (v(i,jc,kc)-v(i,jc, k))*rdz1+(w(i,jc,kc)-w(i,j ,kc))*rdy1 )**2 + &
          ( (v(i,j ,kc)-v(i,j , k))*rdz1+(w(i,j ,kc)-w(i,jb,kc))*rdy1 )**2 ) 
   end do 
end do
	 
	
k=nzm
kc=k+1
kb=k-1

rdz = 1./(dz*adz(k))
rdz2 = 1./(dz*adzw(k))

do j=1,ny
  jb=j-1*YES3D
  jc=j+1*YES3D
  rdx=idx*imu(j) ! take into account grid anisotropy
  rdx1=idx*imuv(jc) ! take into account grid anisotropy
  rdx2=idx*imuv(j) ! take into account grid anisotropy
  rdy=1./(dy*ady(j))
  rdy1=1./(dy*adyv(jc))
  rdy2=1./(dy*adyv(j))
  do i=1,nx
      ib=i-1
      ic=i+1
      def2(i,j,k)=2.* ( &
           ( (u(ic,j,k)-u(i,j,k))*rdx)**2+ &
           ( (v(i,jc,k)-v(i,j,k))*rdy)**2+ &
           ( (w(i,j,kc)-w(i,j,k))*rdz)**2 ) &
       + 0.25 * ( &
           ( (u(ic,jc,k)-u(ic,j ,k))*rdy1+(v(ic,jc,k)-v(i ,jc,k))*rdx1 )**2 +  &
           ( (u(i ,jc,k)-u(i ,j ,k))*rdy1+(v(i ,jc,k)-v(ib,jc,k))*rdx1 )**2 +  &
           ( (u(ic,j ,k)-u(ic,jb,k))*rdy2+(v(ic,j ,k)-v(i ,j ,k))*rdx2 )**2 +  &
           ( (u(i ,j ,k)-u(i ,jb,k))*rdy2+(v(i ,j ,k)-v(ib,j ,k))*rdx2 )**2 )   &
 	+ 0.5 * ( &
           ( (u(ic,j,k )-u(ic,j,kb))*rdz2+(w(ic,j,k )-w(i ,j,k ))*rdx )**2 + &
           ( (u(i ,j,k )-u(i ,j,kb))*rdz2+(w(i ,j,k )-w(ib,j,k ))*rdx )**2 ) &
       + 0.5 * ( &
           ( (v(i,jc,k )-v(i,jc,kb))*rdz2+(w(i,jc,k )-w(i,j ,k ))*rdy2 )**2 + &
           ( (v(i,j ,k )-v(i,j ,kb))*rdz2+(w(i,j ,k )-w(i,jb,k ))*rdy2 )**2 ) 
  end do 
end do
	
end

