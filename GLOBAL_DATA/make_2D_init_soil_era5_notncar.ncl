load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
;grid = "360x180_unif"
;grid = "1440x720_dyvar"
grid = "3840x1920_dyvar"
;grid = "2304x1152_dyvar"
;grid = "3072x3072_CA"
;grid = "3072x3072_SO"
;grid = "1536x1536_AMAZ"
;grid = "9216x128"

;date=2014060100  ; AMAZONIA 
;date=2021102112  ; CA 
;date=2018021512  ; SOCRATES 
;date=2020052000  ; DYAMOND2 
date=2017090500  ; IRMA 
;date=2022092300  ; IAN 

netcdf_out = True

; ERA5 data file
;filedata = "./ERA5/ERA5_SFC_20-27OCT_2021_CA.nc"
filedata = "./ERA5/ERA5_SOIL_SNOW_Sep5_2017.nc"

indir = "NC_D"
outdir = "BIN_D"
outdir_nc = "NC_D"

;============================================================
;============================================================
;============================================================
; don't need to edit below this line

dataset="era5"
fland=addfile("./"+indir+"/landmask_"+grid+".nc","r")


fstl=addfile(filedata,"r")
fswvl=addfile(filedata,"r")

fobs = fstl
times = calendar_decode2(fobs->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(""+dates)
ntime = dimsizes(fobs->time)
itime = -1
do i=0,ntime-1
 if(dates(i).eq.date) then
     itime = i
 end if
end do
if(itime.eq.-1) then
  print("could not find date "+date)
  exit
end if
time = day_of_year(times(itime,0),times(itime,1),times(itime,2))+times(itime,3)/24.
timed = dates(itime)
print("  time="+time+"  "+timed)

lat = fland->lat
lon = fland->lon
nlon = dimsizes(lon)
nlat = dimsizes(lat)
print("nlon="+nlon)
print("nlat="+nlat)
print("lon: "+min(lon(0))+"  "+max(lon))
print("lat: "+min(lat(0))+"  "+max(lat))
land = fland->LANDMASK(:,:)
print("landmask>>>"+min(land)+"  "+max(land))


nnn=dimsizes(fobs->latitude)-1

nsoil = 4
zsoil = (/0.035,0.175,0.64,1.945/) ; midlevels
fld = new((/nsoil,nlat,nlon/),float)
fobs = fstl
do n=0,nsoil-1
 fld1= short2flt(fobs->$("stl"+(n+1))$(itime,nnn:0:1,:))
 if(max(fld1&longitude).lt.0.) then
  fld1&longitude= fld1&longitude+360.
 end if
 if(n.eq.0) then
   print("ERA5 lon:"+ min(fld1&longitude)+"  "+max(fld1&longitude))
   print("ERA5 lat:"+ min(fld1&latitude)+"  "+max(fld1&latitude))
 end if
 print("original stl>>>"+min(fld1)+"  "+max(fld1))
 delete(fld1@missing_value)
 delete(fld1@_FillValue)
; fld(n,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1,lon,lat,False)
 fld(n,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1,True,lon,lat,0)
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
fobs = fswvl
do n=0,nsoil-1
 fld1= short2flt(fobs->$("swvl"+(n+1))$(itime,nnn:0:1,:))
 if(max(fld1&longitude).lt.0.) then
   fld1&longitude= fld1&longitude+360.
 end if
 print("original swvl>>>"+min(fld1)+"  "+max(fld1))
 delete(fld1@missing_value)
 delete(fld1@_FillValue)
; fld(n,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1,lon,lat,False)
 fld(n,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1,True,lon,lat,0)
 fld(n,:,:) = where(land(:,:).eq.0..or.land(:,:).eq.15.,0.,fld(n,:,:))
 fld(n,:,:) = where(fld(n,:,:).lt.0.,0.,fld(n,:,:))
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



