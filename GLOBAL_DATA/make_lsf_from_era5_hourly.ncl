; make snd file from hourly era5 data at NCAR
; over a subregion with given coordinates for a date range
; MK


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
date_start = 2018070900
date_end = 2018071100
date_skip = 1

lat1 = 38.
lat2 = 40.
lon1 = 360.-72.
lon2 = 360.-70.

days = new((/1000/),float)
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
date = year*1000000+month*10000+day*100+hour

end do
;print("ndays = "+ndays)
;print(""+days(0:ndays-1))

strings = new((/1000/),string)

nstrings = 0
strings(nstrings) = ""+ndays
do k=1,ndays
nstrings = nstrings+1
strings(nstrings) = ""+days(k-1)
end do

date = date_start

write_table("lsf", "w", [/"  z[m] p[mb] tpls[K/s] qls[kg/kg/s] uls vls wls"/], "%s")

; surface elevation

fterr = addfile("/glade/campaign/collections/rda/data/d633000/e5.oper.invariant/197901/e5.oper.invariant.128_129_z.ll025sc.1979010100_1979010100.nc","r")

zsfc3 = fterr->Z(0,{lat1:lat2},{lon1:lon2})
zsfc3 = zsfc3/9.81
print("zsfc="+min(zsfc3)+"  "+max(zsfc3))
zsfc = dim_avg_Wrap(dim_avg_Wrap(zsfc3))
print("zsfc mean="+zsfc)

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
;print(""+dates_obs)

dataset = "era5"
dirpath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/"+month_obs+"/"

f_u = addfile(dirpath+"e5.oper.an.pl.128_131_u.ll025uv."+dates_obs+".nc","r")
f_v = addfile(dirpath+"e5.oper.an.pl.128_132_v.ll025uv."+dates_obs+".nc","r")
f_w = addfile(dirpath+"e5.oper.an.pl.128_135_w.ll025sc."+dates_obs+".nc","r")
f_z = addfile(dirpath+"e5.oper.an.pl.128_129_z.ll025sc."+dates_obs+".nc","r")
f_t = addfile(dirpath+"e5.oper.an.pl.128_130_t.ll025sc."+dates_obs+".nc","r")

times = calendar_decode2(f_u->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
;print(""+dates)
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

mmm = itime

; determine pressure levels at z:

nlati = dimsizes(f_u->latitude)
print("latitudes"+f_u->latitude({lat1:lat2}))
presi = doubletofloat(f_u->level(nzi-1:0:1))
z3 = dim_avg_Wrap(dim_avg_Wrap(f_z->Z(mmm,nzi-1:0:1,{lat1:lat2},{lon1:lon2})))
z3 = z3/9.81
print("z3: "+z3+"  presi: "+presi)

print("reading ...")
u=dim_avg_Wrap(dim_avg_Wrap(f_u->U(mmm,nzi-1:0:1,{lat1:lat2},{lon1:lon2})))
v=dim_avg_Wrap(dim_avg_Wrap(f_v->V(mmm,nzi-1:0:1,{lat1:lat2},{lon1:lon2})))
w=dim_avg_Wrap(dim_avg_Wrap(f_w->W(mmm,nzi-1:0:1,{lat1:lat2},{lon1:lon2})))
t=dim_avg_Wrap(dim_avg_Wrap(f_t->T(mmm,nzi-1:0:1,{lat1:lat2},{lon1:lon2})))
print(dimsizes(w))
rho = presi*100./(287.*t)
w = -w/(rho*9.81)
zzz=z3
if(zsfc.gt.zzz(0)) then
 do k=0,nzi-2
  if(zsfc.gt.zzz(k)) then
    pres0 = presi(k)-(presi(k)-presi(k+1))/(zzz(k+1)-zzz(k))*(zsfc-zzz(k))
    break
  end if
 end do
else
 pres0 = presi(0)-(presi(0)-presi(1))/(zzz(1)-zzz(0))*(zsfc-zzz(0))
end if
print("pres0="+pres0)
tmp1=(/time/)
tmp2=(/nzi/)
tmp3=(/pres0/)
zz=zzz-zsfc
zero = presi
zero = 0.
write_table("lsf", "a", [/time,nzi,pres0/], "%f  %i  %f")
write_table("lsf", "a", [/zz,presi,zero,zero,u,v,w/],"%f  %f  %f  %f  %f  %f  %f")

print(zzz+" "+presi+" "+zero+" "+zero+" "+u+" "+v+" "+w)

hour = hour+date_skip
if(hour.gt.23) then
 hour = hour-24
 day = day+1
end if
date = year*1000000+month*10000+day*100+hour

end do
