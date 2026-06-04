load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
;grid = "360x180_unif"
;grid = "1440x720_dyvar"
;grid = "2304x1152_dyvar"
;grid = "3840x1920_dyvar"
;grid = "2000x500_dyvar"
;grid = "4096x2048_IND"
;grid = "3072x3072_CA"
;grid = "3072x3072_SO"
;grid = "1536x1536_AMAZ"
;grid = "9216x128"
;grid="2304x2304_dyvar_COMBLE"
;grid="2048x512_dyvar"
;grid="2048x1024_CONUS"
;grid="3072x3072_HELENE"
;grid="3072x3072_MILTON"
;grid="512x512_MILTON"
;grid="768x768_HELENE"
;grid="1024x1024_HELENE"
;grid="1024x1024_MILTON"
;grid="1024x1024_AGNES"
;grid="2048x2048_AGNES"
;grid="1024x1024_FIONA"
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
;grid="384x384_EUSA"
;grid="768x768_EUSA"
;grid="2048x2048_HUMBERTO"
grid="3456x3456_LI"

;date=2014060100  ; AMAZONIA 
;date=2018021512  ; SOCRATES 
;date=2020052000  ; DYAMOND2 
;date=2022092300  ; IAN 
;date=2021102112  ; CA 
;date=2020050100  ; CLIMO 
;date=2020032800  ;  COMBLE
;date=2024070100  ; USA 
;date=2024092400  ; HELENE 
;date=2024100500  ; MILTON 
;date=2017090400  ; IRMA 
;date=1972061818  ; AGNES 
;date=2022091612  ; FIONA 
;date=2019070100  ; FIONA 
;date=2003091600  ; ISABEL 
;date=2019090500  ; FAXAI 
;date=2022061705  ; HOUSTON 
;date=2020070100  ; USA
;date=2025070300  ;  
;date=2025092700  ;  HUMBERTO
;date=2018070909  ;  LI
date=2025090309  ;  LI

netcdf_out = True
outdir_nc = "NC_D"

indir = "NC_D"
outdir = "BIN_D"

; path to data at NCAR
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc"

; ERA5 land-sea mask
g=addfile("/glade/campaign/collections/rda/data/d633000/e5.oper.invariant/197901/e5.oper.invariant.128_172_lsm.ll025sc.1979010100_1979010100.nc","r")

;========================================================================
;========================================================================
;========================================================================

dataset="era5"

fland=addfile("./"+indir+"/landmask_"+grid+".nc","r")


year = floattointeger(date/1.e6)
month = floattointeger((date-year*1000000)/10000.)
date_obs = (year*1000000+month*10000+100)+"_"+ \
           (year*1000000+month*10000+days_in_month(year,month)*100+23)
year_month = year*100+month

fstl=new((/4/),file)
fstl(0)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_139_stl1.ll025sc."+date_obs+".nc","r")
fstl(1)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_170_stl2.ll025sc."+date_obs+".nc","r")
fstl(2)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_183_stl3.ll025sc."+date_obs+".nc","r")
fstl(3)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_236_stl4.ll025sc."+date_obs+".nc","r")
fswvl=new((/4/),file)
fswvl(0)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_039_swvl1.ll025sc."+date_obs+".nc","r")
fswvl(1)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_040_swvl2.ll025sc."+date_obs+".nc","r")
fswvl(2)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_041_swvl3.ll025sc."+date_obs+".nc","r")
fswvl(3)=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_042_swvl4.ll025sc."+date_obs+".nc","r")

fobs = fstl(0)
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
print("nlat="+nlat)
print("nlon="+nlon)
land = fland->LANDMASK(:,:)
print("landmask>>>"+min(land)+"  "+max(land))

latmin = min(lat)-1.
latmax = max(lat)+1.
lonmin = min(lon)-1.
lonmax = max(lon)+1.

nnn=dimsizes(fobs->latitude)-1

lsm = g->LSM(0,nnn:0:1,:)

nsoil = 4
zsoil = (/0.035,0.175,0.64,1.945/) ; midlevels
fld = new((/nsoil,nlat,nlon/),float)
do n=0,nsoil-1
 fobs = fstl(n)
 fld0= fobs->$("STL"+(n+1))$(itime,nnn:0:1,:)
 fld0 = where(lsm.lt.0.5,-999.,fld0)
 fld1 = fld0({latmin:latmax},{lonmin:lonmax})
 delete(fld1@missing_value)
 delete(fld1@_FillValue)
 fld1@_FillValue = -999.
 poisson_grid_fill(fld1, False, 1, 6000, 1.e-2, 0.6, 0)
 print("original stl>>>"+min(fld1)+"  "+max(fld1))
 delete(fld1@missing_value)
 delete(fld1@_FillValue)
 if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
 end if
; fld(n,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1,lon,lat,False)
 fld(n,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1,True,lon,lat,0)
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
 fobs = fswvl(n)
 fld0= fobs->$("SWVL"+(n+1))$(itime,nnn:0:1,:)
 fld0 = where(lsm.lt.0.5,-999.,fld0)
 fld1 = fld0({latmin:latmax},{lonmin:lonmax})
 delete(fld1@missing_value)
 delete(fld1@_FillValue)
 fld1@_FillValue = -999.
 poisson_grid_fill(fld1, False, 1, 6000, 1.e-2, 0.6, 0)
 print("original swvl>>>"+min(fld1)+"  "+max(fld1))
 delete(fld1@missing_value)
 delete(fld1@_FillValue)
 if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
 end if
; fld(n,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1,lon,lat,False)
 fld(n,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1,True,lon,lat,0)
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

if(netcdf_out) then
 system("rm -f "+outdir_nc+"/soil_init_"+date+"_"+grid+"_"+dataset+".nc")
 out = addfile(outdir_nc+"/soil_init_"+date+"_"+grid+"_"+dataset+".nc","c")
 out->soilt = soilt
 out->soilw = soilw
 out->LSM = lsm({latmin:latmax},{lonmin:lonmax})
 out->fld = fld1({latmin:latmax},{lonmin:lonmax})
 out->fld0 = fld0({latmin:latmax},{lonmin:lonmax})
end if
delete(fld1)


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



