implicit none
real pi, pres, nu, nuv, t, phi,tm, wt
integer i,k, nzm,ny
nzm = 300
ny=361
pi = acos(-1.)
print*,'  z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]'
print*,0.,nzm,1000.
do k=nzm,1,-1
 tm = 0.
 wt = 0.
 do i=1,ny
      phi = (i-1.)/(ny-1.)*90.*pi/180.
      pres = k*1000./nzm
      nu = pres/1000.
      nuv = pi/2.*(nu-0.252)
      if(nu.ge.0.2) then
         t = 288.*nu**(287.*0.005/9.80616)
      else
         t = 288.*nu**(287.*0.005/9.80616)+4.8e5*(0.2-nu)**5
      end if
      t = t + 0.75*nu*pi*35./287.*sin(nuv)*cos(nuv)**0.5 &
             *((-2.*sin(phi)**6*(cos(phi)**2+1./3.)+10./63.)*2.*35.*cos(nuv)**1.5 &
               +(8./5.*cos(phi)**3*(sin(phi)**2+2./3.)-pi/4.)*6371229.*7.29212e-5)
      t=t*(1000./pres)**(287./1004.)
      tm = tm + t*cos(phi)
      wt = wt + cos(phi)
 end do
 tm = tm/wt
 print*,'-999.    ',pres,tm,'   0.    0.    0.'
end do

end

