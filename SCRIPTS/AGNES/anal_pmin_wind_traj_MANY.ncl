ens_max = 1

filename = (/ \
  "../../OUT_2D/AGNES_1024x1024x72_3km_10s_JUN1812_NONUDGE"\
/)

fileout = "ens_pmin_vmax.bin"

system("rm -f "+fileout)
fbinrecwrite(fileout,-1,ens_max)

do ens = 0,ens_max-1

ncfiles=systemfunc("ls -1 "+filename(ens)+"_1972*.2D_atm.nc")

nfiles=dimsizes(ncfiles)
pmin = new((/nfiles/),float)
vel10max = pmin
landc = pmin
day = pmin
xm = pmin
ym = pmin

f = addfiles(ncfiles,"r")
lon = f[0]->lon
lat = f[0]->lat
nlon = dimsizes(lon)
nlat = dimsizes(lat)

; define initial area to search for center:
lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)
print("ens: "+ens)
print("lon:"+lonmin+" "+lonmax)
print("lat:"+latmin+" "+latmax)

; background (time averaged) pressure to make the field more uniform
; over land and terrain:

do i=0,nfiles-1
 p = f[i]->PSFC(0,{latmin:latmax},{lonmin:lonmax})
 if(i.eq.0) then
  p0 = p
  p0 = 0.
 end if
 p0 = p0 + p
end do
p0 = p0/tofloat(nfiles) 

do i=0,nfiles-1

 time = f[i]->time
 day(i) = tofloat(f[i]->day)

 p = f[i]->PSFC(0,{latmin:latmax},{lonmin:lonmax}) 
 p1 = p0({latmin:latmax},{lonmin:lonmax})
 p = p - p1
 land = f[i]->LANDMASK(0,{latmin:latmax},{lonmin:lonmax}) 
 print("land:"+min(land)+"  "+max(land)) 

 ; find coordinates of the minimum pressure:

 dims = dimsizes(p)
 p1D = ndtooned(p)
 inds  = ind_resolve(minind(p1D),dims) 
 delete(p1D)
 ym(i) = tofloat(p&lat(inds(0,0)))
 xm(i) = tofloat(p&lon(inds(0,1)))
 pmin(i) = p(inds(0,0),inds(0,1))+p1(inds(0,0),inds(0,1))
 landc(i) = land(inds(0,0),inds(0,1))
 delete(p)
 delete(p1)
 delete(land)

 ; define box to look for maximum wind:
 latmin = ym(i)-0.5
 latmax = ym(i)+0.5
 lonmin = xm(i)-0.5
 lonmax = xm(i)+0.5

 u = f[i]->U10m(0,{latmin:latmax},{lonmin:lonmax}) 
 v = f[i]->V10m(0,{latmin:latmax},{lonmin:lonmax}) 
 vel10max(i) = sqrt(max(u^2+v^2))
 j = i+1

 print(j+"  "+day(i)+"  "+pmin(i)+"  "+vel10max(i)+"  "+ym(i)+"  "+(xm(i)-360.)+"  "+ landc(i))
 delete(u)
 delete(v)

end do


fbinrecwrite(fileout,-1,j)
fbinrecwrite(fileout,-1,day)
fbinrecwrite(fileout,-1,pmin)
fbinrecwrite(fileout,-1,vel10max)
fbinrecwrite(fileout,-1,xm)
fbinrecwrite(fileout,-1,ym)
fbinrecwrite(fileout,-1,landc)

delete(ncfiles)
delete(p0)
delete(pmin)
delete(vel10max)
delete(landc)
delete(day)
delete(xm)
delete(ym)

end do ; ens

