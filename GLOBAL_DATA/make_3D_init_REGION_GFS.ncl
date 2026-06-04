load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"


date = 2024021618
lonmin = 265.
lonmax = 305.
latmin = 19.
latmax = 51.
identifier = "EUSA"


fobs = addfile("/glade/work/marat/GFS/gfs.0p25."+date+".f000.grib2","r")

outdir = "BIN_D"

netcdf_out = True 
outdir_nc = "NC_D"

;======================================================
dataset = "GFS"


print("date="+date)

nzi = dimsizes(fobs->lv_ISBL0)
grid = identifier

fname = "./"+outdir+"/init_"+dataset+"_"+date+"_"+grid+".bin"
system("rm -f "+fname)

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
if(netcdf_out) then
 gg->v = fld3(:,{latmin:latmax},{lonmin:lonmax})
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
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done"+" OMEGA")
if(netcdf_out) then
 gg->omega = fld3(:,{latmin:latmax},{lonmin:lonmax})
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
if(netcdf_out) then
 gg->t = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)


print("reading Q...")
fld3=fobs->SPFH_P0_L100_GLL0(nzi-1:0:1,nlati-1:0:1,:)
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

nzi6 = dimsizes(fobs->lv_ISBL6)

print("reading QC...")
fld3(0:nzi6-1,:,:)=fobs->CLWMR_P0_L100_GLL0(nzi6-1:0:1,nlati-1:0:1,:)
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

print("reading QCI..")
fld3 = 0.
fld3(0:nzi6-1,:,:)=fobs->ICMR_P0_L100_GLL0(nzi6-1:0:1,nlati-1:0:1,:)
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

print("reading QR...")
fld3 = 0.
fld3(0:nzi6-1,:,:)=fobs->RWMR_P0_L100_GLL0(nzi6-1:0:1,nlati-1:0:1,:)
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

print("reading QS..")
fld3 = 0.
fld3(0:nzi6-1,:,:)=fobs->SNMR_P0_L100_GLL0(nzi6-1:0:1,nlati-1:0:1,:)
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


