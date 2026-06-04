; make soil file from hourly era5 data at NCAR


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
date = 202001
month = 1
daysinmo = 31
range = date+"0100_"+date+""+daysinmo+"23"
day = 18
hour = 1

; chose coordinate of region for avaraging:
lat1 = 53.
lat2 = 55.
lon1 = 283.
lon2 = 286.


filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_139_stl1.ll025sc."+range+".nc"
f_st1 = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_170_stl2.ll025sc."+range+".nc"
f_st2 = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_183_stl3.ll025sc."+range+".nc"
f_st3 = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_236_stl4.ll025sc."+range+".nc"
f_st4 = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_039_swvl1.ll025sc."+range+".nc"
f_sw1 = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_040_swvl2.ll025sc."+range+".nc"
f_sw2 = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_041_swvl3.ll025sc."+range+".nc"
f_sw3 = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_042_swvl4.ll025sc."+range+".nc"
f_sw4 = addfile(filepath,"r")


zsoil_era5 = (/0.,0.035,0.175,0.64,1.945,3./) ; 1st and last added for interpolation
zsoil_thick = (/0.01,0.02,0.04,0.08,0.16,0.32,0.37,0.50,1.00/)
zsoil = zsoil_thick
clay = zsoil
sand = zsoil
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

m = (day-1)*24+hour-1

tsoil_era5(1) = dim_avg(dim_avg(f_st1->STL1(m,{lat1:lat2},{lon1:lon2})))
tsoil_era5(2) = dim_avg(dim_avg(f_st2->STL2(m,{lat1:lat2},{lon1:lon2})))
tsoil_era5(3) = dim_avg(dim_avg(f_st3->STL3(m,{lat1:lat2},{lon1:lon2})))
tsoil_era5(4) = dim_avg(dim_avg(f_st4->STL4(m,{lat1:lat2},{lon1:lon2})))
tsoil_era5(0) = tsoil_era5(1)
tsoil_era5(5) = tsoil_era5(4)
print("era5:"+tsoil_era5)
tsoil = linint1(zsoil_era5,tsoil_era5,False,zsoil,0)
print(""+tsoil)
wsoil_era5(1) = dim_avg(dim_avg(f_sw1->SWVL1(m,{lat1:lat2},{lon1:lon2})))
wsoil_era5(2) = dim_avg(dim_avg(f_sw2->SWVL2(m,{lat1:lat2},{lon1:lon2})))
wsoil_era5(3) = dim_avg(dim_avg(f_sw3->SWVL3(m,{lat1:lat2},{lon1:lon2})))
wsoil_era5(4) = dim_avg(dim_avg(f_sw4->SWVL4(m,{lat1:lat2},{lon1:lon2})))
wsoil_era5(0) = wsoil_era5(1)
wsoil_era5(5) = wsoil_era5(4)
print("era5:"+wsoil_era5)
wsoil = linint1(zsoil_era5,wsoil_era5,False,zsoil,0)
print(""+wsoil)


f = addfile("./DATA/clay_1x1.nc","r")
fld=f->Band1
fld = where(fld.lt.0..or.fld.gt.100.,0.,fld)
fld=lonFlip(fld)
fld@_FillValue = 0.
clay = dim_avg(dim_avg(fld({lat1:lat2},{lon1:lon2})))
print("clay: "+clay)
delete(fld)

f = addfile("./DATA/sand_1x1.nc","r")
fld=f->Band1
fld = where(fld.lt.0..or.fld.gt.100.,0.,fld)
fld=lonFlip(fld)
fld@_FillValue = 0.
sand = dim_avg(dim_avg(fld({lat1:lat2},{lon1:lon2})))
print("sand: "+sand)
delete(fld)

f = addfile("./DATA/LAI_0.25x0.25.nc","r")
name = "Band"+month
fld=f->$name$
fld = where(fld.eq.-88.,1.,fld)
fld = where(fld.eq.-77.,0.,fld)
fld=lonFlip(fld)
fld@_FillValue = 0.
lai = dim_avg(dim_avg(fld({lat1:lat2},{lon1:lon2})))
print("lai: "+lai)

lat = 0.5*(lat1+lat2)
lon = 0.5*(lon1+lon2)
print("latitude0 = "+lat)
print("longitude0 = "+lon)

relax = clay
relax = 0.

system("rm -f soil")
header = (/" thickness[m],  soilt[K], volumetric m3/m3,    SAND%,      CLAY%,   relax function"/)
write_table("soil", "a", [/header/], "%s")
write_table("soil", "a", [/zsoil,tsoil,wsoil,sand,clay,relax/],"%f  %f  %f  %f  %f  %f")

