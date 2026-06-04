load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"


date_start = 2024021618
date_end = 2024021718
lonmin = 265.
lonmax = 305.
latmin = 19.
latmax = 51.
identifier = "EUSA"


outdir = "./BIN_D/"
outdirnc = "./NC_D/"

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
  days(ndays-1) = (year-1)*365.+day_of_year(year,month,day)+hour/24.
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
; make last time a bit longer (just in case)
days(ndays-1) = days(ndays-1)+0.04
print("ndays = "+ndays)
print(""+todouble(days(0:ndays-1)))



strings = new((/5000/),string)

nstrings = 0
strings(nstrings) = ""+ndays
do k=1,ndays
nstrings = nstrings+1
strings(nstrings) = ""+days(k-1)
end do

date = date_start

fname = "./"+outdir+"/ozone_"+dataset+"_"+date_start+"-"+date_end+"_"+identifier+".bin"
system("rm -f "+fname)
if(netcdf) then
  filenamenc = "./"+outdirnc+"/ozone_"+dataset+"_"+date_start+"-"+date_end+"_"+identifier+".nc"
  system("rm -fr "+filenamenc)
  g = addfile(filenamenc,"c")
end if


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
if(date.eq.date_start) then
  fbinrecwrite(fname,-1,inds)
  fbinrecwrite(fname,-1,ndays)
  fbinrecwrite(fname,-1,days(0:ndays-1))
end if


print("reading O3...")
fld3=fobs->O3MR_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("O3: "+min(fld3)+"  "+max(fld3))
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
print("done"+" O3")
if(netcdf) then
  fld3@long_name = "O3"
  fld3@units = "m/s"
  g->O3 = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if

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

