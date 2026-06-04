load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"


date_start = 2024021618
date_end = 2024021718
lonmin = 265.
lonmax = 305.
latmin = 19.
latmax = 51.
identifier = "EUSA"


outdir = "./NUDGE/"

netcdf =  False

dataset = "GFS"

date_skip = 6  ; number of hours to skip
;======================================================
;======================================================
;======================================================
;======================================================

days = new((/2000/),float)
ndays = 0
date = date_start
do while (date.le.date_end)

  ndays = ndays+1
  print(""+date)
  year = floattointeger(date/1.e6)
  month = floattointeger((date-year*1000000)/10000.)
  day = floattointeger((date-year*1000000-month*10000)/100.)
  hour = date-year*1000000-month*10000-day*100
  days(ndays-1) = day_of_year(year,month,day)+hour/24.
  hour = hour+date_skip
  if(hour.gt.23) then
   hour = hour-24
   day = day+1
  end if
  if(day.gt.days_in_month(year,month)) then
   month = month+1
   day = 1
   hour = 0
  end if
  date = year*1000000+month*10000+day*100+hour

end do
days(ndays-1) = days(ndays-1)+0.01
print("ndays = "+ndays)
print(""+days(0:ndays-1))

strings = new((/5000/),string)

nstrings = 0
strings(nstrings) = ""+ndays
do k=1,ndays
nstrings = nstrings+1
strings(nstrings) = ""+days(k-1)
end do

date = date_start


do while (date.le.date_end)

fobs = addfile("/glade/work/marat/GFS/gfs.0p25."+date+".f000.grib2","r")

print(""+date)
year = floattointeger(date/1.e6)
month = floattointeger((date-year*1000000)/10000.)
day = floattointeger((date-year*1000000-month*10000)/100.)
hour = date-year*1000000-month*10000-day*100
print(""+(date-year*1000000-month*10000)/100.)
print(year+" "+month+" "+day+" "+hour)
print("day_of_year="+(day_of_year(year,month,day)+hour/24.))


nzi = dimsizes(fobs->lv_ISBL0)
grid = identifier
print("data grid:"+grid)

filename = "/nudge_"+dataset+"_"+date+"_"+identifier+".bin"
fname = outdir+filename
system("rm -f "+fname)
nstrings = nstrings+1
strings(nstrings) = filename
if(netcdf) then
  filenamenc = "./NUDGE/"+"nudge_"+dataset+"_"+date+"_"+identifier+".nc"
  system("rm -fr "+filenamenc)
  g = addfile(filenamenc,"c")
end if


; determine pressure levels at z:

nlati = dimsizes(fobs->lat_0)
presi = fobs->lv_ISBL0(nzi-1:0:1)*0.01
z3 = fobs->HGT_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 z3=lonFlip(z3)
end if
lat = doubletofloat(z3&lat_0)
lon = doubletofloat(z3&lon_0)
lat!0 = "lat"
lon!0 = "lon"
lat&lat = lat
lon&lon = lon
z33 = z3(:,{latmin:latmax},{lonmin:lonmax})
wgt = cos(z33&lat_0*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z33,wgt,0,1))
print("z3m: "+z3m+"  presi: "+presi)
nlat = dimsizes(z33&lat_0)
nlon = dimsizes(z33&lon_0)
inds = (/nlon,nlat,nzi/)
fbinrecwrite(fname,-1,inds)


print("reading U...")
fld3=fobs->UGRD_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("u: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon({lonmin:lonmax}))
fbinrecwrite(fname,-1,lat({latmin:latmax}))
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,{latmin:latmax},{lonmin:lonmax}))
end do
print("done"+" U")
if(netcdf) then
  fld3@long_name = "U"
  fld3@units = "m/s"
  g->U = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

print("reading V...")
fld3=fobs->VGRD_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("v: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon({lonmin:lonmax}))
fbinrecwrite(fname,-1,lat({latmin:latmax}))
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,{latmin:latmax},{lonmin:lonmax}))
end do
print("done"+" V")
if(netcdf) then
  fld3@long_name = "V"
  fld3@units = "m/s"
  g->V = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

print("reading W...")
fld3=fobs->VVEL_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("w: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon({lonmin:lonmax}))
fbinrecwrite(fname,-1,lat({latmin:latmax}))
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,{latmin:latmax},{lonmin:lonmax}))
end do
print("done"+" W")
if(netcdf) then
  fld3@long_name = "W"
  fld3@units = "Pa/s"
  g->W = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

print("reading T...")
fld3=fobs->TMP_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("t: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon({lonmin:lonmax}))
fbinrecwrite(fname,-1,lat({latmin:latmax}))
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,{latmin:latmax},{lonmin:lonmax}))
end do
print("done"+" T")
if(netcdf) then
  fld3@long_name = "T"
  fld3@units = "K"
  g->T = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

nzi6 = dimsizes(fobs->lv_ISBL6)

print("reading Q...")
fld3=fobs->SPFH_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
fld3(0:nzi6-1,:,:) = fld3(0:nzi6-1,:,:) + fobs->CLWMR_P0_L100_GLL0(nzi6-1:0:1,nlati-1:0:1,:)
fld3(0:nzi6-1,:,:) = fld3(0:nzi6-1,:,:) + fobs->ICMR_P0_L100_GLL0(nzi6-1:0:1,nlati-1:0:1,:)
print("q: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon({lonmin:lonmax}))
fbinrecwrite(fname,-1,lat({latmin:latmax}))
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,{latmin:latmax},{lonmin:lonmax}))
end do
print("done"+" Q")
if(netcdf) then
  fld3@long_name = "Q"
  fld3@units = "kg/kg"
  g->Q = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

hour = hour+date_skip
if(hour.gt.23) then
 hour = hour-24
 day = day+1
end if
if(day.gt.days_in_month(year,month)) then
 month = month+1
 day = 1
 hour = 0
end if

date = year*1000000+month*10000+day*100+hour

end do

fname = outdir+"/nudge_"+dataset+"_"+date_start+"-"+date_end+"_"+identifier+"_"+"FILELIST.ascii"
asciiwrite(fname,strings(0:nstrings))
