! compute acceleration due to Coriolis force

subroutine coriolis

use vars
use params, only:  dowallx,dowally,docoriolis,dolatlon,docoriolisz,dometric,docolumn
use check_energy

implicit none
	
integer i,j,k,jc,jb,it,jt
real iady, iadz
real u_av, v_av
	
call t_startf("coriolis")

if(RUN3D) then

  call sumenergy(1,'coriolis')

  if(docoriolis) then 

   if(dolatlon) then
     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
      do j=1,ny
       iady = 1./ady(j)
       jb = j
       jc = j+1
       if(j+jt.eq.ny_gl) jc=ny
       if(j+jt.eq.1) jb=2
       do i=1,nx
         dudt(i,j,k,na)=dudt(i,j,k,na) + 0.25*(fcory(j)+u(i,j,k)*tanr(j))*iady* &
            (adyv(jb)*(v(i,jb,k)+v(i-1,jb,k))+adyv(jc)*(v(i,jc,k)+v(i-1,jc,k)))
         dvdt(i,j,k,na)=dvdt(i,j,k,na) - 0.25*imuv(j)* &
           ((fcory(j)+u(i,j,k)*tanr(j))*mu(j)*(u(i,j,k)) + &
            (fcory(j)+u(i+1,j,k)*tanr(j))*mu(j)*(u(i+1,j,k)) + &
            (fcory(j-1)+u(i,j-1,k)*tanr(j-1))*mu(j-1)*(u(i,j-1,k)) + &
            (fcory(j-1)+u(i+1,j-1,k)*tanr(j-1))*mu(j-1)*(u(i+1,j-1,k)))
        end do ! i
      end do ! j
     end do ! k
   else
     do k=1,nzm
      do j=1,ny
       do i=1,nx
         v_av=0.25*(v(i,j,k)+v(i,j+1,k)+v(i-1,j,k)+v(i-1,j+1,k))
         dudt(i,j,k,na)=dudt(i,j,k,na)+fcory(j)*(v_av-vg0(k))
         u_av=0.25*(u(i,j,k)+u(i+1,j,k)+u(i,j-1,k)+u(i+1,j-1,k))
         dvdt(i,j,k,na)=dvdt(i,j,k,na)-0.5*(fcory(j)+fcory(j-1))*(u_av-ug0(k))
       end do ! i
      end do ! j
     end do ! k
   end if ! dolatlon 

   if(docoriolisz) then
     do k=1,nzm
      iadz = 1./adz(k)
      do j=1,ny
       do i=1,nx
         dudt(i,j,k,na)=dudt(i,j,k,na) - 0.25*fcorzy(j)*iadz* &
           (adzw(k)*(w(i,j,k)+w(i-1,j,k))+adzw(k+1)*(w(i,j,k+1)+w(i-1,j,k+1)))
       end do ! i
      end do ! j
     end do ! k
     do k=2,nzm
      do j=1,ny
       do i=1,nx
         dwdt(i,j,k,na)=dwdt(i,j,k,na) + 0.25*fcorzy(j)* &
           (u(i,j,k)+u(i+1,j,k)+u(i,j,k-1)+u(i+1,j,k-1))
       end do ! i
      end do ! j
     end do ! k
   end if

  else

   if(dolatlon.and.dometric) then
     call task_rank_to_index(rank,it,jt)
     do k=1,nzm
      do j=1,ny
       iady = 1./ady(j)
       jb = j
       jc = j+1
       if(j+jt.eq.ny_gl) jc=ny
       if(j+jt.eq.1) jb=2
       do i=1,nx
         dudt(i,j,k,na)=dudt(i,j,k,na) + 0.25*u(i,j,k)*tanr(j)*iady* &
         (adyv(jb)*(v(i,jb,k)+v(i-1,jb,k))+adyv(jc)*(v(i,jc,k)+v(i-1,jc,k))-4.*vg0(k))
         dvdt(i,j,k,na)=dvdt(i,j,k,na) - 0.25*imuv(j)* &
           (u(i,j,k)*tanr(j)*mu(j)*(u(i,j,k)-ug0(k)) + &
            u(i+1,j,k)*tanr(j)*mu(j)*(u(i+1,j,k)-ug0(k)) + &
            u(i,j-1,k)*tanr(j-1)*mu(j-1)*(u(i,j-1,k)-ug0(k)) + &
            u(i+1,j-1,k)*tanr(j-1)*mu(j-1)*(u(i+1,j-1,k)-ug0(k)))
       end do ! i
      end do ! j
     end do ! k
     end if ! dolatlon

  end if ! docoriolis

  call sumenergy(-1,'coriolis')

else  ! 2D run or single-column run

  do k=1,nzm
   do j=1,ny
    do i=1,nx
     dudt(i,j,k,na)=dudt(i,j,k,na)+fcory(j)*(v(i,j,k)-vg0(k))
     dvdt(i,j,k,na)=dvdt(i,j,k,na)-fcory(j)*(u(i,j,k)-ug0(k))
    end do ! i
   end do ! i
  end do ! k

endif ! 3D or 2D

!-------------------------------------------------------------------------------
	
!bloss: accumulate coriolis accelerations for statistics
if(dostatis) then

   utendcor(:) = 0.
   vtendcor(:) = 0.

   if(RUN3D) then

      do k=1,nzm
         do j=1,ny
            do i=1,nx
               utendcor(k)=utendcor(k)+fcory(j)*(v(i,j,k)-vg0(k))
               vtendcor(k)=vtendcor(k)-fcory(j)*mu(j)*imuv(j)*(u(i,j,k)-ug0(k))
            end do ! i
         end do ! j
      end do ! k

   else

      do k=1,nzm
         do j=1,ny
            do i=1,nx
               utendcor(k)=utendcor(k)+fcory(j)*(v(i,j,k)-vg0(k))
               vtendcor(k)=vtendcor(k)-fcory(j)*(u(i,j,k)-ug0(k))
            end do ! i
         end do ! i
      end do ! k

   endif !if(RUN3D)

end if !if(dostatis)

call t_stopf("coriolis")

end subroutine coriolis

