load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

;grid="768x768_EUSA"
grid="384x384_EUSA"

date=2024021618  ;  

netcdf_out = True
outdir_nc = "NC_D"

indir = "NC_D"
outdir = "BIN_D"

fobs = addfile("/glade/work/marat/GFS/gfs.0p25."+date+".f000.grib2","r")

;========================================================================
;========================================================================
;========================================================================

dataset="GFS"

fland=addfile("./"+indir+"/landmask_"+grid+".nc","r")


lat = fland->lat
lon = fland->lon
nlon = dimsizes(lon)
nlat = dimsizes(lat)
print("nlat="+nlat)
print("nlon="+nlon)
land = fland->LANDMASK(:,:)
print("landmask>>>"+min(land)+"  "+max(land))

nnn=dimsizes(fobs->lat_0)-1

nsoil = 4
zsoil = (/0.05,0.25,0.7,1.5/) ; midlevels 
fld = new((/nsoil,nlat,nlon/),float)
do n=0,nsoil-1
 fld1= fobs->TSOIL_P0_2L106_GLL0(n,nnn:0:1,:)
 print("original stl>>>"+min(fld1)+"  "+max(fld1))
 poisson_grid_fill(fld1, True, 1, 2000, 1.e-2, 0.6, 0)
 if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
 end if
 fld(n,:,:) = linint2_Wrap(fld1&lon_0,fld1&lat_0,fld1,True,lon,lat,0)
 delete(fld@_FillValue)
 print("interpolated soilt>>>"+min(fld(n,:,:))+"  "+max(fld(n,:,:)))
end do
fld!0="zsoil"
fld&zsoil = zsoil
fld&zsoil@units = "m"
fld!1="lat"
fld!2="lon"
fld&lon = lon
fld&lat = lat
soilt = fld
do n=0,nsoil-1
 fld1= fobs->SOILW_P0_2L106_GLL0(n,nnn:0:1,:)
 print("original swvl>>>"+min(fld1)+"  "+max(fld1))
 poisson_grid_fill(fld1, True, 1, 2000, 1.e-2, 0.6, 0)
 if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
 end if
 fld(n,:,:) = linint2_Wrap(fld1&lon_0,fld1&lat_0,fld1,True,lon,lat,0)
 fld(n,:,:) = where(land(:,:).eq.0..or.land(:,:).eq.15.,0.,fld(n,:,:))
 fld(n,:,:) = where(fld(n,:,:).lt.0.,0.,fld(n,:,:))
 delete(fld@_FillValue)
 print("interpolated swvl>>>"+min(fld(n,:,:))+"  "+max(fld(n,:,:)))
end do
fld!0="zsoil"
fld&zsoil = zsoil
fld&zsoil@units = "m"
fld!1="lat"
fld!2="lon"
fld&lon = lon
fld&lat = lat
soilw = fld
delete(fld1)

if(netcdf_out) then
 system("rm -f "+outdir_nc+"/soil_init_"+date+"_"+grid+"_"+dataset+".nc")
 out = addfile(outdir_nc+"/soil_init_"+date+"_"+grid+"_"+dataset+".nc","c")
 out->soilt = soilt
 out->soilw = soilw
end if


fld1 = soilt
fld2 = soilw
print("date:"+date)
print("soilt:"+min(fld1(:,:,:))+"  "+max(fld1(:,:,:)))
print("soilw:"+min(fld2(:,:,:))+"  "+max(fld2(:,:,:)))
fname = "./"+outdir+"/soil_init_"+date+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
fbinrecwrite(fname,-1,nsoil)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,zsoil)
do i=1,nsoil
  fbinrecwrite(fname,-1,fld1(i-1,:,:))
end do
do i=1,nsoil
  fbinrecwrite(fname,-1,fld2(i-1,:,:))
end do



