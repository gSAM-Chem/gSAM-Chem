! wind based on "A Microscale Model for Air Pollutant Dispersion Simulation in Urban Areas: Presentation of the Model and Performance over a Single Building" Znahg et al (2016), Adv. in Atmos. Sci.
implicit none
real, parameter :: zmax = 150. ! m
real, parameter :: uref = 6. ! m/s
real, parameter :: href = 100. ! m
real, parameter :: dz = 1.
real z, u
integer n

z = 0.
n = 0
do while(z.le.zmax)
 n = n+1
 u = uref*(z/href)**0.21
 write(*,'(6f15.5)') z,-999.,300.,0.,u,0.
 z = z + dz 
end do
print*,'number of levels:',n

end
