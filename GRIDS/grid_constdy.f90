implicit none
integer, parameter :: ny = 180
real(8), parameter :: latmax = 90.
real(8), parameter :: dy = 2.*latmax/real(ny)
real(8), parameter :: dym = 20000.*latmax/90./real(ny)
real(8) y(ny)
integer j
y(1)=-latmax+0.5_8*dy
print*,y(1),1
do j=2,ny
 y(j) = y(j-1)+dy
 print*,y(j),j,dym
end do
print*,y
end
