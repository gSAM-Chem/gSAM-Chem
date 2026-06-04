filename = (/ \
  "../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_NUDGEUV12_2018" \
/)

g = addfile("/glade/work/marat/ERA5_Hourly/NOR/ERA5_SFC_JAN3-5_2018_NOR_hourly.nc","r")


nstart = 0
nend = 30

nstart_obs = 0
nend_obs = nend+nstart_obs

ens_max = dimsizes(filename)

fileout = "ens_pmin_vmax.bin"

system("rm -f "+fileout)
fbinrecwrite(fileout,-1,ens_max+1)

do ens = 0,ens_max-1

ncfiles=systemfunc("ls -1 "+filename(ens)+"*.2D_atm.nc")
print("ncfiles = "+ncfiles)

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
lonmin0 = -80.+360.
lonmax0 = -75.+360.
latmin0 = 26.
latmax0 = 32.
lonmin = lonmin0
lonmax = lonmax0
latmin = latmin0
latmax = latmax0
print("ens: "+ens)
print("lon:"+lonmin+" "+lonmax)
print("lat:"+latmin+" "+latmax)

; background (time averaged) pressure to make the field more uniform
; over land and terrain:


do i= nstart,nend
 p = f[i]->PSFC(0,:,:)
 if(i.eq.nstart) then
  p0 = p
  p0 = 0.
 end if
 p0 = p0 + p
end do
p0 = p0/tofloat(nend-nstart+1) 
delete(p)

j=0
do i=nstart,nend

 time = f[i]->time
 day(j) = tofloat(f[i]->day)

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
 ym(j) = tofloat(p&lat(inds(0,0)))
 xm(j) = tofloat(p&lon(inds(0,1)))
 pmin(j) = p(inds(0,0),inds(0,1))+p1(inds(0,0),inds(0,1))
 landc(j) = land(inds(0,0),inds(0,1))
 delete(p)
 delete(p1)
 delete(land)

 ; define box to look for maximum wind:
 latmin = ym(j)-1.5
 latmax = ym(j)+1.5
 lonmin = xm(j)-1.5
 lonmax = xm(j)+1.5

 u = f[i]->U10m(0,{latmin:latmax},{lonmin:lonmax}) 
 v = f[i]->V10m(0,{latmin:latmax},{lonmin:lonmax}) 
 vel10max(j) = sqrt(max(u^2+v^2))

 print((j+1)+"  "+day(j)+"  "+pmin(j)+"  "+vel10max(j)+"  "+ym(j)+"  "+(xm(j)-360.)+"  "+ landc(j))
 delete(u)
 delete(v)
 j = j+1

end do


fbinrecwrite(fileout,-1,j)
fbinrecwrite(fileout,-1,day)
fbinrecwrite(fileout,-1,pmin)
fbinrecwrite(fileout,-1,vel10max)
fbinrecwrite(fileout,-1,xm)
fbinrecwrite(fileout,-1,ym)

delete(ncfiles)
delete(p0)
;delete(pmin)
;delete(vel10max)
delete(landc)
;delete(day)
;delete(xm)
;delete(ym)

end do ; ens

;===============================================================
; Observed track

lonmin = lonmin0-360.
lonmax = lonmax0-360.
latmin = latmin0
latmax = latmax0

do i= nstart_obs,nend_obs
 p = g->sp(i,:,:)
 if(i.eq.nstart_obs) then
  p0 = p
  p0 = 0.
 end if
 p0 = p0 + p
end do
p0 = p0/tofloat(nend-nstart+1)
delete(p)

j=0
do i=nstart_obs,nend_obs

 pp = g->sp(i,::-1,:)
 p = pp({latmin:latmax},{lonmin:lonmax})
 p1 = p0({latmin:latmax},{lonmin:lonmax})
 p = p - p1

 ; find coordinates of the minimum pressure:

 dims = dimsizes(p)
 p1D = ndtooned(p)
 inds  = ind_resolve(minind(p1D),dims)
 delete(p1D)
 ym(j) = tofloat(p&latitude(inds(0,0)))
 xm(j) = tofloat(p&longitude(inds(0,1)))
 pmin(j) = 0.01*(p(inds(0,0),inds(0,1))+p1(inds(0,0),inds(0,1)))
 delete(p)
 delete(p1)

 ; define box to look for maximum wind:
 latmin = ym(j)-1.5
 latmax = ym(j)+1.5
 lonmin = xm(j)-1.5
 lonmax = xm(j)+1.5

 u = g->u10(0,{latmin:latmax},{lonmin:lonmax})
 v = g->v10(0,{latmin:latmax},{lonmin:lonmax})
 vel10max(j) = sqrt(max(u^2+v^2))

 print((j+1)+"  "+day(j)+"  "+pmin(j)+"  "+vel10max(j)+"  "+ym(j)+"  "+xm(j))
 delete(u)
 delete(v)
 j = j+1

end do

fbinrecwrite(fileout,-1,j)
fbinrecwrite(fileout,-1,day)
fbinrecwrite(fileout,-1,pmin)
fbinrecwrite(fileout,-1,vel10max)
fbinrecwrite(fileout,-1,xm)
fbinrecwrite(fileout,-1,ym)

