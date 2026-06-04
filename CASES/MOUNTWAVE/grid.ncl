z=0.
dz0 = 20.

n=1
dz = dz0
print(z+"  "+n)


z1 = 1000.
c1 = 1.
do while(z.le.z1)
 z = z + dz 
 n = n + 1
 print(z+"  "+n+"  "+dz)
 dz = dz*c1
end do

z1 = 3000.
c1 = 1.2
do while(z.le.z1)
 z = z + dz 
 n = n + 1
 print(z+"  "+n+"  "+dz)
 dz = dz*c1
end do

z1 = 18000.
c1 = 1.
do while(z.le.z1)
 z = z + dz 
 n = n + 1
 print(z+"  "+n+"  "+dz)
 dz = dz*c1
end do

z1 = 37000.
c1 = 1.05
do while(z.le.z1)
 z = z + dz 
 n = n + 1
 print(z+"  "+n+"  "+dz)
 dz = dz*c1
end do
