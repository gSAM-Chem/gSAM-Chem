load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

;date = 2017083000
;date = 2017090400
;date = 2017090500
;lonmin = 250.
;lonmax = 340.
;latmin = -5.
;latmax = 40.
;lonmin = 270.
;lonmax = 315.
;latmin = 0.
;latmax = 40.
;identifier = "IRMA"

;date = 2024100400
;lonmin = 250.
;lonmax = 294.
;latmin = 10.
;latmax = 48.
;identifier = "MILTON"

;date = 2019070100
;lonmin = 219.
;lonmax = 302.
;latmin = 20.
;latmax = 54.
;identifier = "CONUS"

;date = 1972061818
;lonmin = 264.
;lonmax = 299.
;latmin = 21.
;latmax = 49.
;identifier = "AGNES"

;date = 2022061705
;lonmin = 260.
;lonmax = 272.
;latmin = 25.
;latmax = 35.
;identifier = "HOUSTON"

;date = 2019090500
;lonmin = 122.
;lonmax = 171.
;latmin = 9.
;latmax = 51.
;identifier = "FAXAI"

;date = 2003091600
;lonmin = 264.
;lonmax = 299.
;latmin = 21.
;latmax = 49.
;identifier = "ISABEL"

;date = 2022091612
;lonmin = 284.
;lonmax = 305.
;latmin = 7.
;latmax = 27.
;identifier = "FIONA"

;date = 2024092400
;lonmin = 250.
;lonmax = 294.
;latmin = 10.
;latmax = 48.
;identifier = "HELENE"

; COMBLE
;date = 2020032800
;lonmin = -20.
;lonmax = 50.
;latmin = 60.
;latmax = 85.
;identifier = "COMBLE"

;date = 2020070100
;lonmin = 216.
;lonmax = 314.
;latmin = 11.
;latmax = 61.
;identifier = "USA"

;date = 2024090200
;lonmin = 287.
;lonmax = 360.
;latmin = -8.
;latmax = 28.
;identifier = "ITCZ"

;date = 2025070300
;date = 2019041000
;lonmin = 252.
;lonmax = 274.
;latmin = 26.
;latmax = 44.
;identifier = "USA_SM"
;;lonmin = 241.
;lonmax = 285.
;latmin = 19.
;latmax = 55.
;identifier = "USA_MID"
;lonmin = 219.
;lonmax = 307.
;latmin = 19.
;latmax = 55.
;identifier = "USA"

;date = 2024021618
;lonmin = 265.
;lonmax = 305.
;latmin = 19.
;latmax = 51.
;identifier = "EUSA"

;date = 2025092700
;lonmin = 289.
;lonmax = 311.
;latmin = 13.
;latmax = 33.
;identifier = "HUMBERTO"

date = 2025090309
;date = 2018070909
lonmin = 283.5
lonmax = 288.5
latmin = 39.
latmax = 43.
identifier = "LI"



outdir = "BIN_D"

netcdf_out = False 
outdir_nc = "NC_D"

;======================================================
dataset = "era5"


year = floattointeger(date/1.e6)
month = floattointeger((date-year*1000000)/10000.)
day = floattointeger((date-year*1000000-month*10000)/100.)
dates_obs = (year*1000000+month*10000+day*100)+"_"+ \
           (year*1000000+month*10000+day*100+23)
year_month = year*100+month

dirpath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.pl/"+year_month+"/"
f_u = addfile(dirpath+"e5.oper.an.pl.128_131_u.ll025uv."+dates_obs+".nc","r")
f_v = addfile(dirpath+"e5.oper.an.pl.128_132_v.ll025uv."+dates_obs+".nc","r")
f_w = addfile(dirpath+"e5.oper.an.pl.128_135_w.ll025sc."+dates_obs+".nc","r")
f_z = addfile(dirpath+"e5.oper.an.pl.128_129_z.ll025sc."+dates_obs+".nc","r")
f_t = addfile(dirpath+"e5.oper.an.pl.128_130_t.ll025sc."+dates_obs+".nc","r")
f_q = addfile(dirpath+"e5.oper.an.pl.128_133_q.ll025sc."+dates_obs+".nc","r")
f_clwc = addfile(dirpath+"e5.oper.an.pl.128_246_clwc.ll025sc."+dates_obs+".nc","r")
f_ciwc = addfile(dirpath+"e5.oper.an.pl.128_247_ciwc.ll025sc."+dates_obs+".nc","r")
f_crwc = addfile(dirpath+"e5.oper.an.pl.128_075_crwc.ll025sc."+dates_obs+".nc","r")
f_cswc = addfile(dirpath+"e5.oper.an.pl.128_076_cswc.ll025sc."+dates_obs+".nc","r")

times = calendar_decode2(f_u->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(dates)
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
time = times(itime,0)*365+day_of_year(times(itime,0),times(itime,1),times(itime,2))+times(itime,3)/24.
print("  time="+time)

nzi = dimsizes(f_u->level)
grid = identifier

fname = "./"+outdir+"/init_"+dataset+"_"+date+"_"+grid+".bin"
system("rm -f "+fname)

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
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
fbinrecwrite(fname,-1,lon({lonmin:lonmax}))
fbinrecwrite(fname,-1,lat({latmin:latmax}))
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,{latmin:latmax},{lonmin:lonmax}))
end do
print("done"+" U")
if(netcdf_out) then
 fnamenc = outdir_nc+"/init_"+dataset+"_"+date+"_"+grid+".nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->u = fld3(:,{latmin:latmax},{lonmin:lonmax})
 gg->z = z3m
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
if(netcdf_out) then
 gg->v = fld3(:,{latmin:latmax},{lonmin:lonmax})
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
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" OMEGA")
if(netcdf_out) then
 gg->omega = fld3(:,{latmin:latmax},{lonmin:lonmax})
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
if(netcdf_out) then
 gg->t = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)


print("reading Q...")
fld3=f_q->Q(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
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
print("done"+" QV")
if(netcdf_out) then
 gg->qv = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

print("reading QC...")
fld3=f_clwc->CLWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("qc: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QC")
if(netcdf_out) then
 gg->qc = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

print("reading QCI..")
fld3=f_ciwc->CIWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("qi: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QI")
if(netcdf_out) then
 gg->qi = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

print("reading QR...")
fld3=f_crwc->CRWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("qr: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QR")
if(netcdf_out) then
 gg->qr = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

print("reading QS..")
fld3=f_cswc->CSWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("qs: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QS")
if(netcdf_out) then
 gg->qs = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if


