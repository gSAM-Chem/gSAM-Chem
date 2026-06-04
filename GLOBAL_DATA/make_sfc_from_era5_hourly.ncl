; make soil file from hourly era5 data at NCAR


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
date = 201807
month = 7
daysinmo = 31
range = date+"0100_"+date+""+daysinmo+"23"
day_start = 9
hour_start = 0
day_end = 11
hour_end = 0

; chose coordinate of region for avaraging:
lat1 = 38.
lat2 = 40.
lon1 = 360.-72.
lon2 = 360.-70.

filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_235_skt.ll025sc."+range+".nc"
f_skt = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_232_ie.ll025sc."+range+".nc"
f_lhf = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_231_ishf.ll025sc."+range+".nc"
f_shf = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_229_iews.ll025sc."+range+".nc"
f_taux = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_230_inss.ll025sc."+range+".nc"
f_tauy = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_167_2t.ll025sc."+range+".nc"
f_t2m = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_168_2d.ll025sc."+range+".nc"
f_td2m = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_165_10u.ll025sc."+range+".nc"
f_u10m = addfile(filepath,"r")
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc/"+date+"/e5.oper.an.sfc.128_166_10v.ll025sc."+range+".nc"
f_v10m = addfile(filepath,"r")

m1=(day_start-1)*24+hour_start
m2=(day_end-1)*24+hour_end

time = f_skt->time(m1:m2)
dates = tointeger(cd_calendar(time, 0))
doy = day_of_year(dates(:,0),dates(:,1),dates(:,2))+dates(:,3)/24.

skt = dim_avg(dim_avg(f_skt->SKT(m1:m2,{lat1:lat2},{lon1:lon2})))
lhf = -dim_avg(dim_avg(f_lhf->IE(m1:m2,{lat1:lat2},{lon1:lon2})))
shf = -dim_avg(dim_avg(f_shf->ISHF(m1:m2,{lat1:lat2},{lon1:lon2})))
lhf = lhf*2.5e6*100000./(287.*skt(0))
taux = dim_avg(dim_avg(f_taux->IEWS(m1:m2,{lat1:lat2},{lon1:lon2})))
tauy = dim_avg(dim_avg(f_tauy->INSS(m1:m2,{lat1:lat2},{lon1:lon2})))
t2m = dim_avg(dim_avg(f_t2m->VAR_2T(m1:m2,{lat1:lat2},{lon1:lon2})))
td2m = dim_avg(dim_avg(f_td2m->VAR_2D(m1:m2,{lat1:lat2},{lon1:lon2})))
u10m = dim_avg(dim_avg(f_u10m->VAR_10U(m1:m2,{lat1:lat2},{lon1:lon2})))
v10m = dim_avg(dim_avg(f_v10m->VAR_10V(m1:m2,{lat1:lat2},{lon1:lon2})))
tau = sqrt(taux^2+tauy^2)
vel10m = sqrt(u10m^2+v10m^2)


print(doy+"  "+skt+"  "+shf+"  "+lhf+" "+tau+" "+t2m+" "+td2m+" "+vel10m)

write_table("sfc", "w", [/" day sst(K) H(W/m2) LE(W/m2) TAU(m2/s2) T2m(K) Tdew2m(K) Wind10m"/], "%s")
write_table("sfc", "a", [/doy,skt,shf,lhf,tau,t2m,td2m,vel10m/],"%f  %f  %f  %f  %f  %f  %f  %f")

