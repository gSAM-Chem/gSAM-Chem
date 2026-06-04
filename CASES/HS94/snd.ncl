; make initial sounding:

pi = acos(-1.)
phi = 30.*pi/180.
dp = 10. ; pressure step in mb
p0 = 1000.
kmax = p0/dp

print("  z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]")
print(0.+"   "+kmax+"   "+p0)
do k=kmax,1,1
 p = k*dp
 t = (315.-60.*sin(phi)^2-10.*log(p/p0)*cos(phi)^2)*(p/p0)^0.286
 if(t.lt.200.) then
   t=200.
 end if 
 t=t*(p0/p)^0.286
 print("-999.  "+p+"   "+t+"   "+"0.   0.   0.")
end do


