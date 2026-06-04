dz1 = 40.
dz2 = 500.
dz3 = 1000.
sgrid1 = 1.1 
sgrid2 = 1.1 
z2=18000.
z=0.
zmax = 35000.
k=1
dz = dz1
print(z+"  "+k)
do while(z.le.zmax)
 z = z + dz
 k = k+1
 print(z+"  "+k+"  "+dz)
 if(z.lt.z2) then
  if(dz*sgrid1.lt.dz2) then 
   dz = dz*sgrid1
  else
   dz = dz2
  end if
 else
  if(dz*sgrid2.lt.dz3) then 
   dz = dz*sgrid1
  else
   dz = dz3
  end if
 end if
end do
