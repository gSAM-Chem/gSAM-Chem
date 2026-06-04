; make isothermal atmosphere with constant wind
t0=300.
u0=20.
dp=10.

p=1000.
do while(p.gt.0.)
 print("-999.   "+p+"   "+(-t0)+"   "+"0. "+u0+"   0.")
 p = p-dp
end do

 
