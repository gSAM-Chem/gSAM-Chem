load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;grid="1536x1536_AMAZ"
;grid="3072x3072_SO"
;grid="2000x500_dyvar"
;grid="4096x2048_IND"
;grid="3072x3072_CA"
;grid="360x180_unif"
;grid="1440x720_dyvar"
;grid="9216x128"
;grid="3840x1920_dyvar"
;grid="2304x2304_dyvar_COMBLE"
;grid="2048x1024_CONUS"
;grid="3072x3072_HELENE"
;grid="2048x2048_ITCZ"
;grid="768x768_HELENE"
;grid="1024x1024_IRMA"
;grid="1024x1024_AGNES"
;grid="2048x2048_AGNES"
;grid="1024x1024_FIONA"
;grid="2048x2048_AGNES"
;grid="1024x1024_MILTON"
;grid="3072x3072_MILTON"
;grid="512x512_MILTON"
;grid="1024x1024_HELENE"
;grid="2560x1280_ATL"
;grid="512x128_TUNE"
;grid="2048x1024_IRMA"
;grid="1024x1024_ISABEL"
;grid="1536x1536_FAXAI"
;grid="2048x2048_HOUSTON"
;grid="1440x720_GLOBAL"
;grid="384x192_USA"
;grid="2048x1024_ITCZ"
;grid="6144x3072_ITCZ"
;grid="6144x3072_USA"
;grid="3072x1536_USA"
;grid="1536x1536_USA"
;grid="2048x1024_USA"
;grid="768x768_USA_SM"
;grid="12800x6400_GLOBAL"
;grid="384x384_EUSA"
;grid="768x768_EUSA"
;grid="2048x2048_HUMBERTO"
grid="3456x3456_LI"

indir = "./NC_D"
outdir = "./BIN_D"
outdir_nc = "./NC_D"

netcdf_out = False

;=====================================================
;=====================================================
;=====================================================

fclay = "./DATA/clay_1x1.nc"
fsand = "./DATA/sand_1x1.nc"
flai = "./DATA/LAI_0.25x0.25.nc"

time_lai = (/15.,46.,74.,105.,135.,166.,196.,227.,258.,288.,319.,349./)

f=addfile(indir+"/landtype_"+grid+".nc","r")
fld1 = f->LAND
fld = fld1*1.
lat = f->lat
lon = f->lon
fld!0="lat"
fld!1="lon"
fld&lat = lat
fld&lon = lon
nlon=dimsizes(lon)
nlat=dimsizes(lat)
print(nlon+"  "+nlat)
land = fld
print("land:"+min(fld)+"  "+max(fld))
system("rm -f "+outdir+"/landtype_"+grid+".bin")
fbinrecwrite(outdir+"/landtype_"+grid+".bin",-1,1)
fbinrecwrite(outdir+"/landtype_"+grid+".bin",-1,nlon)
fbinrecwrite(outdir+"/landtype_"+grid+".bin",-1,nlat)
fbinrecwrite(outdir+"/landtype_"+grid+".bin",-1,0.)
fbinrecwrite(outdir+"/landtype_"+grid+".bin",-1,land)
delete(fld)
delete(fld1)

landmask = land
landmask@long_name="Landmask 0-ocean, 1-land"
landmask = 0.
landmask = where(land.ne.0.,1.,landmask)
system("rm -f  "+outdir_nc+"/landmask_"+grid+".nc")
g=addfile(outdir_nc+"/landmask_"+grid+".nc","c")
g->LANDMASK = landmask
fld = landmask
print("landmask:"+min(fld)+"  "+max(fld))
system("rm -f "+outdir+"/landmask_"+grid+".bin")
fbinrecwrite(outdir+"/landmask_"+grid+".bin",-1,1)
fbinrecwrite(outdir+"/landmask_"+grid+".bin",-1,nlon)
fbinrecwrite(outdir+"/landmask_"+grid+".bin",-1,nlat)
fbinrecwrite(outdir+"/landmask_"+grid+".bin",-1,0.)
fbinrecwrite(outdir+"/landmask_"+grid+".bin",-1,fld)

f=addfile(fclay,"r")
fld1=f->Band1
fld1 = where(fld1.lt.0..or.fld1.gt.100.,0.,fld1)
if(.not.any(land&lon.lt.0.)) then
  fld1=lonFlip(fld1)
end if
fld1@_FillValue = 0.
fld@_FillValue = 0.
poisson_grid_fill(fld1, True, 0, 1500, 1.e-2, 0.6, 0)
system("rm -f "+outdir+"/soil_"+grid+".bin")
fbinrecwrite(outdir+"/soil_"+grid+".bin",-1,nlon)
fbinrecwrite(outdir+"/soil_"+grid+".bin",-1,nlat)
delete(fld1@_FillValue)
fld =linint2_Wrap(fld1&lon,fld1&lat,fld1,True,land&lon,land&lat,0)
if(netcdf_out) then
 system("rm -f "+outdir_nc+"/sand_clay_"+grid+".nc")
 g=addfile(outdir_nc+"/sand_clay_"+grid+".nc","c")
end if
fld = where(landmask.ne.1.,0.,fld)
delete(fld&lon)
delete(fld&lat)
fld!0="lat"
fld!1="lon"
fld&lat = land&lat
fld&lon = land&lon
fld@long_name="Soil Clay percentage"
print("clay:"+min(fld)+"  "+max(fld))
if(netcdf_out) then
 g->CLAY = fld
end if
fbinrecwrite(outdir+"/soil_"+grid+".bin",-1,fld)


f=addfile(fsand,"r")
fld1=f->Band1
fld1 = where(fld1.lt.0..or.fld1.gt.100.,0.,fld1)
if(.not.any(land&lon.lt.0.)) then
  fld1=lonFlip(fld1)
end if
fld1@_FillValue = 0.
poisson_grid_fill(fld1, True, 0, 1500, 1.e-2, 0.6, 0)
delete(fld1@_FillValue)
delete(fld@_FillValue)
fld =linint2_Wrap(fld1&lon,fld1&lat,fld1,True,land&lon,land&lat,0)
fld = where(landmask.ne.1.,0.,fld)
delete(fld&lon)
delete(fld&lat)
fld!0="lat"
fld!1="lon"
fld&lat = land&lat
fld&lon = land&lon
fld@long_name="Soil Sand percentage"
print("sand:"+min(fld)+"  "+max(fld))
if(netcdf_out) then
 g->SAND = fld
end if
fbinrecwrite(outdir+"/soil_"+grid+".bin",-1,fld)
delete(fld)
delete(fld1)

if(netcdf_out) then
 system("rm -f "+outdir_nc+"/LAI_"+grid+".nc")
 g=addfile(outdir_nc+"/LAI_"+grid+".nc","c")
end if
f=addfile(flai,"r")
system("rm -f "+outdir+"/lai_"+grid+".bin")
fbinrecwrite(outdir+"/lai_"+grid+".bin",-1,dimsizes(time_lai))
fbinrecwrite(outdir+"/lai_"+grid+".bin",-1,nlon)
fbinrecwrite(outdir+"/lai_"+grid+".bin",-1,nlat)
fbinrecwrite(outdir+"/lai_"+grid+".bin",-1,time_lai)
do i=1,12
 name = "Band"+i
 fld1=f->$name$
 fld1 = where(fld1.eq.-88.,1.,fld1)
 fld1 = where(fld1.eq.-77.,0.,fld1)
 if(.not.any(land&lon.lt.0.)) then
   fld1=lonFlip(fld1)
 end if
 fld1@_FillValue = 0.
; print("lai in data before poisson:"+min(fld1)+"  "+max(fld1))
 poisson_grid_fill(fld1, True, 0, 1500, 1.e-2, 0.6, 0)
; print("lai in data after poisson:"+min(fld1)+"  "+max(fld1))
 delete(fld1@_FillValue)
; print("fld1&lon:  "+min(fld1&lon)+"   "+max(fld1&lon))   
; print("fld1&lat:  "+min(fld1&lat)+"   "+max(fld1&lat))   
; print("land&lon:  "+min(land&lon)+"   "+max(land&lon))   
; print("land&lat:  "+min(land&lat)+"   "+max(land&lat))   
 fld =linint2_Wrap(fld1&lon,fld1&lat,fld1,True,land&lon,land&lat,0)
 fld = where(landmask.ne.1.,0.,fld)
 fld = where(fld.lt.0.,0.,fld)
 fld = where(land.eq.16,0.,fld)
 delete(fld@_FillValue)
 delete(fld&lon)
 delete(fld&lat)
 fld!0="lat"
 fld!1="lon"
 fld&lat = land&lat
 fld&lon = land&lon
 fld@long_name="LAI Month "+i
 name = "LAI"+i
 if(netcdf_out) then
   g->$name$ = fld
 end if
 print("month: "+i+"  lai:"+min(fld)+"  "+max(fld))
 fbinrecwrite(outdir+"/lai_"+grid+".bin",-1,fld)
end do
delete(fld)
delete(fld1)

delete(lat)
delete(lon)
lat=land&lat
lon=land&lon

f=addfile(outdir_nc+"/terrain_"+grid+".nc","r")
fld = f->terr_smooth
lat = f->lat
lon = f->lon
delete(fld&lon)
delete(fld&lat)
fld!0="lat"
fld!1="lon"
fld&lat = lat
fld&lon = lon
nlon=dimsizes(lon)
nlat=dimsizes(lat)
print(nlon+"  "+nlat)
print("terr:"+min(fld)+"  "+max(fld))
system("rm -f "+outdir+"/terrain_"+grid+".bin")
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,1)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,nlon)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,nlat)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,0.)
fbinrecwrite(outdir+"/terrain_"+grid+".bin",-1,fld)

delete(fld)

