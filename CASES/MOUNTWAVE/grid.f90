real zz(100)
dz0=20.
z0=0.
dz1=50.
z1=1000.
dz2=100.
z2=2000.
dz3=500.
z3=4000.
dz4=500.
z4=18000.
dz5=1500.
z5=25000.
dz6=1500.
zmax=37000.

k=1
zz(1)=0.
print*,zz(1),k
do while(z.lt.zmax)
 if(z.le.z1) then
  dz = dz0+(z-z0)/(z1-z0)*(dz1-dz0)
 else if(z.le.z2) then 
  dz = dz1+(z-z1)/(z2-z1)*(dz2-dz1)
 else if(z.le.z3) then 
  dz = dz2+(z-z2)/(z3-z2)*(dz3-dz2)
 else if(z.le.z4) then
  dz = dz3+(z-z3)/(z4-z3)*(dz4-dz3)
 else if(z.le.z5) then
  dz = dz4+(z-z4)/(z5-z4)*(dz5-dz4)
 else 
  dz = dz6
 end if
 z = z + dz
 k = k+1
 print*,z,k,dz
 zz(k) = z
end do
print*,(zz(i),', ',i=1,k)

end
