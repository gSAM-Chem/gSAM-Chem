load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

; IAN
;date = 2022092300

; IRMA
date = 2017090500

; AMAZONIA
;date = 2014060100

; CA
;date = 2021102112

; SOCRATES
;date = 2018021512

;identifier = "CA"
identifier = "GLOBAL"

netcdf_out = True 

; ERA5 data file
;filedata = "./ERA5/ERA5_PRES_20-27OCT_2021_CA.nc"
filedata = "./ERA5/ERA5_PRES_Sep5_2017.nc"

outdir_nc = "NC_D"
outdir = "BIN_D"

;============================================
;============================================
;============================================

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
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.0.) then
  lon = lon+360.
end if
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
fld3=short2flt(f_v->v(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("v: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
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
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_V.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->v = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading W...")
fld3=short2flt(f_w->w(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("w: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
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
fld3=short2flt(f_t->t(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("t: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.180.) then
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
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+dates(mmm)+"_"+grid+"_T.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->t = fld3
 gg->z = z3m
end if
delete(fld3)


print("reading Q...")
fld3=short2flt(f_q->q(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("q: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
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
fld3=short2flt(f_clwc->clwc(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("qc: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
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
fld3=short2flt(f_ciwc->ciwc(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("qi: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
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

