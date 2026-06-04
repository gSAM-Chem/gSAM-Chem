; This script creates 3D initial dataset
; Works only on NCAR's computers 
; Marat Khairoutdinov

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; IAN
;date = 2022092300

; IRMA
;date = 2017090500

; AMAZONIA
;date = 2014060100

; CA
;date = 2021102112

; SOCRATES
;date = 2018021512

;IND 
;date = 2020051500

;CLIMO
;date = 2023070100

;
date = 2020012000

identifier = "GLOBAL"

netcdf_out = False 

outdir = "BIN_D"
outdir_nc = "NC_D"

;=======================================================================
;=======================================================================
;=======================================================================

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

nlat = dimsizes(f_u->latitude)
nlon = dimsizes(f_u->longitude)
nzi = dimsizes(f_u->level)
grid = identifier
print("data grid:"+grid)

fname = "./"+outdir+"/init_"+dataset+"_"+date+"_"+grid+".bin"
system("rm -f "+fname)
inds = (/nlon,nlat,nzi/)
fbinrecwrite(fname,-1,inds)

mmm = itime
print(dates(mmm)+"  time="+todouble(time))

; determine pressure levels at z:

nlati = dimsizes(f_u->latitude)
presi = doubletofloat(f_u->level(nzi-1:0:1))
z3 = f_z->Z(mmm,nzi-1:0:1,nlati-1:0:1,:)
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
wgt = cos(z3&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z3,wgt,0,1))/9.80665
print("z3m: "+z3m+"  presi: "+presi)

print("reading U...")
fld3=f_u->U(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("u: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" U")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_U.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->u = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading V...")
fld3=f_v->V(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("v: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
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
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_V.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->v = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading W...")
fld3=f_w->W(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("w: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" OMEGA")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_W.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->omega = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading T...")
fld3=f_t->T(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("t: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
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
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_T.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->t = fld3
 gg->z = z3m
end if
delete(fld3)


print("reading Q...")
fld3=f_q->Q(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("q: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QV")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_Q.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qv = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QC...")
fld3=f_clwc->CLWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("qc: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QC")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_QC.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qc = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QCI..")
fld3=f_ciwc->CIWC(mmm,nzi-1:0:1,nlati-1:0:1,:)
print("qi: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QI")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_QI.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qi = fld3
 gg->z = z3m
end if

fld3 = 0.
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QR")

fld3 = 0.
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" QS")
delete(fld3)

