; create initial datasets from ERA5 ensembles downloaded from Copernicus site.
; Note that the pressure levels in the data is reverse than analysis, that it from 1000
; down to small. Always check the order of pressure levels, as it has changed for
; some reason.


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

dataset = "era5"

; AGNES
date = 2003091600
filedata = "/glade/work/marat/ERA5_ENS/ISABEL/era5_isabel_ens_091600.nc"
identifier = "ISABEL"
; AGNES
;date = 1972061812
;filedata = "/glade/work/marat/ERA5_ENS/AGNES/era5_agnes_ens_061812.nc"
;identifier = "AGNES"
; IRMA
;date = 2017083000
;filedata = "/glade/derecho/scratch/marat/ERA5_ENS/IRMA/era5_irma_ens_0830.nc"
;identifier = "era5_IRMA"
; HELENE
;date = 2024092400
;filedata = "/glade/derecho/scratch/marat/ERA5_ENS/HELENE/0924/era5_helene_ens_0924.nc"
;identifier = "era5_HELENE"
; MILTON
;date = 2024100500
;filedata = "/glade/derecho/scratch/marat/ERA5_ENS/MILTON/1005/data_stream-enda_stepType-instant.nc"
;identifier = "era5_MILTON"

add360 = True

outdir = "BIN_D"

netcdf_out = True
outdir_nc = "NC_D"

f = addfile(filedata,"r")

nlat = dimsizes(f->latitude)
nlon = dimsizes(f->longitude)
nzi = dimsizes(f->pressure_level)

do nens = 0,9

fname = "./"+outdir+"/init_"+dataset+"_"+date+"_"+identifier+"_0"+nens+".bin"
system("rm -f "+fname)
inds = (/nlon,nlat,nzi/)
fbinrecwrite(fname,-1,inds)

mmm = 0

; determine pressure levels at z:

nlati = dimsizes(f->latitude)
presi = doubletofloat(f->pressure_level(:))
z3 = f->z(nens,mmm,:,nlati-1:0:1,:)
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
if(add360) then
 lon = lon + 360.
end if
print("lat:"+min(lat)+"  "+max(lat))
print("lon:"+min(lon)+"  "+max(lon))
wgt = cos(z3&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z3,wgt,0,1))/9.80665
print("z3m: "+z3m+"  presi: "+presi)

print("reading U...")
fld3=f->u(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_U.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->u = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading V...")
fld3=f->v(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_V.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->v = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading W...")
fld3=f->w(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_W.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->omega = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading T...")
fld3=f->t(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_T.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->t = fld3
 gg->z = z3m
end if
delete(fld3)


print("reading Q...")
fld3=f->q(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_Q.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qv = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QC...")
fld3=f->clwc(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_QC.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qc = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QCI..")
fld3=f->ciwc(nens,mmm,:,nlati-1:0:1,:)
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
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_QI.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qi = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QR...")
fld3=f->crwc(nens,mmm,:,nlati-1:0:1,:)
print("qr: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QR")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_QR.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qr = fld3
 gg->z = z3m
end if
delete(fld3)

print("reading QS...")
fld3=f->cswc(nens,mmm,:,nlati-1:0:1,:)
print("qs: "+min(fld3)+"  "+max(fld3))
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
print("done"+" QS")
if(netcdf_out) then
 fnamenc = outdir_nc+"/"+dataset+"_"+date+"_"+"_0"+nens+"_QS.nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->qs = fld3
 gg->z = z3m
end if
delete(fld3)


end do ; nens
