load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; CA
date_start = 2021102112
date_end = 2021102700
date_skip = 1

identifier = "CA"

outdir = "./NUDGE/"

netcdf =  True

; ERA5 data file
filedata = "./ERA5/ERA5_PRES_20-27OCT_2021_CA.nc"

;==================================================
;==================================================
;==================================================

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

print(""+date)
year = floattointeger(date/1.e6)
month = floattointeger((date-year*1000000)/10000.)
day = floattointeger((date-year*1000000-month*10000)/100.)
hour = date-year*1000000-month*10000-day*100
print(""+(date-year*1000000-month*10000)/100.)
print(year+" "+month+" "+day+" "+hour)
print("day_of_year="+(day_of_year(year,month,day)+hour/24.))

dates_obs = (year*1000000+month*10000+day*100)+"_"+(year*1000000+month*10000+day*100+23)
month_obs = ""+(year*100+month)
print(""+dates_obs)

dataset = "era5"
f_u = addfile(filedata,"r")
f_v = addfile(filedata,"r")
f_w = addfile(filedata,"r")
f_z = addfile(filedata,"r")
f_t = addfile(filedata,"r")
f_q = addfile(filedata,"r")
f_clwc = addfile(filedata,"r")
f_ciwc = addfile(filedata,"r")

times = calendar_decode2(f_u->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(""+dates)
ntime = dimsizes(f_u->time)
itime = -1
do i=0,ntime-1
 if(dates(i).eq.date) then
  if(itime.eq.-1) then
     itime = i
  end if
 end if
end do
if(itime.eq.-1) then
  print("could not find date "+date)
  exit
end if
time = day_of_year(times(itime,0),times(itime,1),times(itime,2))+times(itime,3)/24.
print("  time="+time)

nlat = dimsizes(f_u->latitude)
nlon = dimsizes(f_u->longitude)
nzi = dimsizes(f_u->level)
grid = identifier
print("data grid:"+grid)

filename = "/nudge_era5_"+date+"_"+grid+".bin"
fname = outdir+filename
system("rm -f "+fname)
nstrings = nstrings+1
strings(nstrings) = filename
if(netcdf) then
  filenamenc = "./NUDGE/"+"nudge_era5_"+date+"_"+grid+".nc"
  system("rm -fr "+filenamenc)
  g = addfile(filenamenc,"c")
end if


inds = (/nlon,nlat,nzi/)
fbinrecwrite(fname,-1,inds)

mmm = itime

; determine pressure levels at z:

nlati = dimsizes(f_u->latitude)
presi = doubletofloat(f_u->level(nzi-1:0:1))
z3 = short2flt(f_z->z(mmm,nzi-1:0:1,nlati-1:0:1,:))
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
wgt = cos(z3&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z3,wgt,0,1))/9.80665
print("z3m: "+z3m+"  presi: "+presi)

print("reading U...")
fld3=short2flt(f_u->u(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("u: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.0.) then
  lon = lon+360.
end if
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" U")
if(netcdf) then
  fld3@long_name = "U"
  fld3@units = "m/s"
  g->U = fld3
end if
delete(fld3)

print("reading V...")
fld3=short2flt(f_v->v(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("v: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.0.) then
  lon = lon+360.
end if
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" V")
if(netcdf) then
  fld3@long_name = "V"
  fld3@units = "m/s"
  g->V = fld3
end if
delete(fld3)

print("reading W...")
fld3=short2flt(f_w->w(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("w: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.0.) then
  lon = lon+360.
end if
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" W")
if(netcdf) then
  fld3@long_name = "W"
  fld3@units = "Pa/s"
  g->W = fld3
end if
delete(fld3)

print("reading T...")
fld3=short2flt(f_t->t(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("t: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.0.) then
  lon = lon+360.
end if
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" T")
if(netcdf) then
  fld3@long_name = "T"
  fld3@units = "K"
  g->T = fld3
end if
delete(fld3)

print("reading Q...")
fld3=short2flt(f_q->q(mmm,nzi-1:0:1,nlati-1:0:1,:))
fld3 = fld3 + short2flt(f_clwc->clwc(mmm,nzi-1:0:1,nlati-1:0:1,:))
fld3 = fld3 + short2flt(f_ciwc->ciwc(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("q: "+min(fld3)+"  "+max(fld3))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.0.) then
  lon = lon+360.
end if
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" Q")
if(netcdf) then
  fld3@long_name = "Q"
  fld3@units = "kg/kg"
  g->Q = fld3
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

fname = outdir+"/nudge_era5_"+date_start+"-"+date_end+"_"+grid+"_FILELIST.ascii"
asciiwrite(fname,strings(0:nstrings))
