load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

zsoil_era5 = (/0.,0.035,0.175,0.64,1.945,3./) ; 1st and last added for interpolation
zsoil_thick = (/0.01,0.02,0.04,0.08,0.16,0.32,0.37,0.50,1.00/)
zsoil = zsoil_thick
clay = zsoil
sand = zsoil
clay = 13.
sand = 17.
tsoil_era5 = zsoil_era5
wsoil_era5 = tsoil_era5
tsoil = zsoil
wsoil = tsoil
nsoil = dimsizes(zsoil_thick)
z = 0.
do k=0,nsoil-1
  z = z+zsoil_thick(k)
  zsoil(k) = z - 0.5*zsoil_thick(k)
  print(k+" "+zsoil(k))
end do
time_soil = 22.

f = addfile("lev_era5.nc","r")
fs = addfile("sfc_era5.nc","r")
times = calendar_decode2(f->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(dates)
year = floattointeger(dates/1.e6)
month = floattointeger((dates-year*1000000)/10000.)
day = floattointeger((dates-year*1000000-month*10000)/100.)
hour = dates-year*1000000-month*10000-day*100
days = day_of_year(year,month,day)+hour/24.
print(days)

nzi = dimsizes(f->level)
presi = doubletofloat(f->level(nzi-1:0:1))
print("  presi: "+presi)

z3 = dim_avg(dim_avg(short2flt(f->z(:,nzi-1:0:1,:,:))))
z3 = z3/9.80665

print("reading ...")
u=dim_avg(dim_avg(short2flt(f->u(:,nzi-1:0:1,:,:))))
v=dim_avg(dim_avg(short2flt(f->v(:,nzi-1:0:1,:,:))))
w=dim_avg(dim_avg(short2flt(f->w(:,nzi-1:0:1,:,:))))
t=dim_avg(dim_avg(short2flt(f->t(:,nzi-1:0:1,:,:))))
q=dim_avg(dim_avg(short2flt(f->q(:,nzi-1:0:1,:,:))))
q=q+dim_avg(dim_avg(short2flt(f->ciwc(:,nzi-1:0:1,:,:))))
q=q+dim_avg(dim_avg(short2flt(f->crwc(:,nzi-1:0:1,:,:))))
q=q+dim_avg(dim_avg(short2flt(f->cswc(:,nzi-1:0:1,:,:))))

write_table("snd", "w", [/"  z[m] p[mb] tp[K] q[g/kg] u[m/s] v[m/s]"/], "%s")
write_table("lsf", "w", [/"  z[m] p[mb] tpls[K/s] qls[kg/kg/s] uls vls wls"/], "%s")
ntime = dimsizes(days)
sp=dim_avg(dim_avg(short2flt(fs->sp(:,:,:))))/100.
do i=0,ntime-1

  zzz=z3(i,:)
  pres0 = presi(0)+(presi(0)-presi(1))/(zzz(1)-zzz(0))*zzz(0)
;  pres0 = sp(i)
  rho = presi*100./(297.*t(i,:))
  w(i,:) = -w(i,:)/(9.91*rho)
  rho = 0.
  rho1 = rho

  write_table("snd", "a", [/days(i),nzi,pres0/], "%f  %i  %f")
  write_table("snd", "a", [/zzz(:),presi,(t(i,:)*(1000./presi)^(287./1004.)),q(i,:)*1.e3,u(i,:),v(i,:)/],"%f  %f  %f  %f  %f  %f")
  write_table("lsf", "a", [/days(i),nzi,pres0/], "%f  %i  %f")
  write_table("lsf", "a", [/zzz(:),presi,rho,rho1,u(i,:),v(i,:),w(i,:)/],"%f  %f  %f  %f  %f  %f  %f")
end do

do i=0,ntime-1
if(days(i).eq.time_soil) then
 tsoil_era5(1) = dim_avg(dim_avg(short2flt(fs->stl1(i,:,:))))
 tsoil_era5(2) = dim_avg(dim_avg(short2flt(fs->stl2(i,:,:))))
 tsoil_era5(3) = dim_avg(dim_avg(short2flt(fs->stl3(i,:,:))))
 tsoil_era5(4) = dim_avg(dim_avg(short2flt(fs->stl4(i,:,:))))
 tsoil_era5(0) = tsoil_era5(1)
 tsoil_era5(5) = tsoil_era5(4)
 tsoil = linint1(zsoil_era5,tsoil_era5,False,zsoil,0) 
 print(""+tsoil)
 wsoil_era5(1) = dim_avg(dim_avg(short2flt(fs->swvl1(i,:,:))))
 wsoil_era5(2) = dim_avg(dim_avg(short2flt(fs->swvl2(i,:,:))))
 wsoil_era5(3) = dim_avg(dim_avg(short2flt(fs->swvl3(i,:,:))))
 wsoil_era5(4) = dim_avg(dim_avg(short2flt(fs->swvl4(i,:,:))))
 wsoil_era5(0) = wsoil_era5(1)
 wsoil_era5(5) = wsoil_era5(4)
 wsoil = linint1(zsoil_era5,wsoil_era5,False,zsoil,0) 
 print(""+wsoil)
end if
end do

write_table("soil", "w", [/"  thickness[m],  soilt[K],   wetness,     SAND,        CLAY,   relax function"/], "%s")
relax = zsoil
relax = 0.
write_table("soil", "a", [/zsoil_thick,tsoil,wsoil,sand,clay,relax/],"%f  %f  %f  %f  %f  %f  %f")
