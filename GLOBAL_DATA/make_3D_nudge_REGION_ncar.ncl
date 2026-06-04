load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

date_skip = 1

;date_start = 1972061700
;date_end = 1972062400
;lonmin = 264.
;lonmax = 299.
;latmin = 21.
;latmax = 49.
;identifier = "AGNES"

;date_start = 2022061705
;date_end = 2022061807
;lonmin = 260.
;lonmax = 272.
;latmin = 25.
;latmax = 35.
;identifier = "HOUSTON"

;date_start = 2024090200
;date_end = 2024090401
;lonmin = 287.
;lonmax = 360.
;latmin = -8.
;latmax = 28.
;identifier = "ITCZ"

;date_start = 2025070300
;date_end = 2025070500
;date_start = 2019041000
;date_end = 2019041600
;lonmin = 252.
;lonmax = 274.
;latmin = 26.
;latmax = 44.
;identifier = "USA_SM"
;lonmin = 241.
;lonmax = 285.
;latmin = 19.
;latmax = 55.
;identifier = "USA_MID"
;lonmin = 219.
;lonmax = 307.
;latmin = 19.
;latmax = 55.
;identifier = "USA"

;date_start = 2020070100
;date_end = 2020070500
;lonmin = 216.
;lonmax = 314.
;latmin = 11.
;latmax = 61.
;identifier = "USA"

;date_start = 2019090400
;date_end = 2019091100
;lonmin = 122.
;lonmax = 171.
;latmin = 9.
;latmax = 51.
;identifier = "FAXAI"

;date_start = 2003091500
;date_end = 2003092000
;lonmin = 264.
;lonmax = 299.
;latmin = 21.
;latmax = 49.
;identifier = "ISABEL"

;date_start = 2019070100
;date_end = 2019071100
;lonmin = 219.
;lonmax = 302.
;latmin = 20.
;latmax = 54.
;identifier = "CONUS"

;date_start = 2022091600
;date_end = 2022092200
;lonmin = 284.
;lonmax = 305.
;latmin = 7.
;latmax = 27.
;identifier = "FIONA"

;date_start = 2017083000
;date_end = 2017091200
;lonmin = 250.
;lonmax = 340.
;latmin = -5.
;latmax = 40.
;lonmin = 270.
;lonmax = 315.
;latmin = 0.
;latmax = 40.
;identifier = "IRMA_LG"

;lonmin = 310.
;lonmax = 340.
;latmin = -5.
;latmax = 25.
;identifier = "ITCZ"

;date_start = 2024092100
;date_end = 2024092900
;lonmin = 250.
;lonmax = 294.
;latmin = 10.
;latmax = 48.
;identifier = "HELENE"

;date_start = 2024100400
;date_end = 2024101200
;lonmin = 250.
;lonmax = 300.
;latmin = 10.
;latmax = 50.
;identifier = "MILTON"

;date_start = 2024021600
;date_end = 2024021800
;lonmin = 265.
;lonmax = 305.
;latmin = 19.
;latmax = 51.
;identifier = "EUSA"

;date_start = 2025092700
;date_end = 2025092900
;lonmin = 289.
;lonmax = 311.
;latmin = 13.
;latmax = 33.
;identifier = "HUMBERTO"

date_start = 2025090309
date_end = 2025090410
;date_start = 2018070909
;date_end = 2018071010
lonmin = 283.5
lonmax = 288.5
latmin = 39.
latmax = 43.
identifier = "LI"


outdir = "./NUDGE/"

netcdf =  False

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
dirpath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/"+month_obs+"/"
f_u = addfile(dirpath+"e5.oper.an.pl.128_131_u.ll025uv."+dates_obs+".nc","r")
f_v = addfile(dirpath+"e5.oper.an.pl.128_132_v.ll025uv."+dates_obs+".nc","r")
f_w = addfile(dirpath+"e5.oper.an.pl.128_135_w.ll025sc."+dates_obs+".nc","r")
f_z = addfile(dirpath+"e5.oper.an.pl.128_129_z.ll025sc."+dates_obs+".nc","r")
f_t = addfile(dirpath+"e5.oper.an.pl.128_130_t.ll025sc."+dates_obs+".nc","r")
f_q = addfile(dirpath+"e5.oper.an.pl.128_133_q.ll025sc."+dates_obs+".nc","r")
f_clwc = addfile(dirpath+"e5.oper.an.pl.128_246_clwc.ll025sc."+dates_obs+".nc","r")
f_ciwc = addfile(dirpath+"e5.oper.an.pl.128_247_ciwc.ll025sc."+dates_obs+".nc","r")

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

nzi = dimsizes(f_u->level)
grid = identifier
print("data grid:"+grid)

filename = "/nudge_era5_"+date+"_"+identifier+".bin"
fname = outdir+filename
system("rm -f "+fname)
nstrings = nstrings+1
strings(nstrings) = filename
if(netcdf) then
  filenamenc = "./NUDGE/"+"nudge_era5_"+date+"_"+identifier+".nc"
  system("rm -fr "+filenamenc)
  g = addfile(filenamenc,"c")
end if


mmm = itime

; determine pressure levels at z:

nlati = dimsizes(f_u->latitude)
presi = doubletofloat(f_u->level(nzi-1:0:1))
z3 = f_z->Z(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 z3=lonFlip(z3)
end if
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
lat!0 = "lat"
lon!0 = "lon"
lat&lat = lat
lon&lon = lon
z33 = z3(:,{latmin:latmax},{lonmin:lonmax})
wgt = cos(z33&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z33,wgt,0,1))/9.80665
print("z3m: "+z3m+"  presi: "+presi)
nlat = dimsizes(z33&latitude)
nlon = dimsizes(z33&longitude)
inds = (/nlon,nlat,nzi/)
fbinrecwrite(fname,-1,inds)


print("reading U...")
fld3=f_u->U(mmm,nzi-1:0:1,nlati-1:0:1,:)
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
fld3=f_v->V(mmm,nzi-1:0:1,nlati-1:0:1,:)
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
fld3=f_w->W(mmm,nzi-1:0:1,nlati-1:0:1,:)
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
fld3=f_t->T(mmm,nzi-1:0:1,nlati-1:0:1,:)
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

print("reading Q...")
fld3=f_q->Q(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
fld3 = fld3 + f_clwc->CLWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
fld3 = fld3 + f_ciwc->CIWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
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

fname = outdir+"/nudge_era5_"+date_start+"-"+date_end+"_"+identifier+"_"+"FILELIST.ascii"
asciiwrite(fname,strings(0:nstrings))
