;              W A R N I N G !;
; This script use the convention that day is year*365+calday of that year.
;
;==========================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
;grid = "360x180_unif"
;grid = "1440x720_dyvar"
;grid = "2304x1152_dyvar"
grid = "3840x1920_dyvar"
;grid = "6144x3072_SO"
;grid = "3072x3072_CA"
;grid = "9216x128"

; IAN
;date=2022092300
; IRMA
date=2017090500
; 
; CA
;date=2021102112        

; ERA5 data file
;filedata = "./ERA5/ERA5_SFC_20-27OCT_2021_CA.nc"
filedata = "./ERA5/ERA5_SOIL_SNOW_Sep5_2017.nc"

indir = "./NC_D"
outdir = "BIN_D"

netcdf_out = True
outdir_nc = "NC_D"

;=======================================================
;=======================================================

fland=addfile(indir+"/landmask_"+grid+".nc","r")

dataset="era5"

date_start = date
date_end = date

nskip_sd=1
fsd=addfile(filedata,"r")
fst=addfile(filedata,"r")

fsfc = fsd
times = calendar_decode2(fsfc->time,-5)
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(""+dates)
ntime = dimsizes(fsfc->time)
itime1 = -1
do i=0,ntime-1
 if(dates(i).ge.date_start.and.dates(i).le.date_end) then
  if(itime1.eq.-1) then
     itime1 = i
  end if
  itime2 = i
 end if
end do
if(itime1.eq.-1) then
  print("could not find date "+date_start+" "+date_end)
  exit
end if
ntime = itime2-itime1+1
time = new((/ntime/),float)
timed = new((/ntime/),integer)
do i=itime1,itime2
  time(i-itime1) = times(i,0)*365+day_of_year(times(i,0),times(i,1),times(i,2))+times(i,3)/24.
  timed(i-itime1) = dates(i)
  print("  time="+time(i-itime1)+"  "+timed(i-itime1))
end do

if(ntime.eq.1) then
  date_out = ""+date_start
else
  date_out = ""+date_start+"-"+date_end
end if

lat = fland->lat
lon = fland->lon
nlon = dimsizes(lon)
nlat = dimsizes(lat)
print("nlat="+nlat)
print("nlon="+nlon)
print("lon: "+min(lon(0))+"  "+max(lon))
print("lat: "+min(lat(0))+"  "+max(lat))

land = fland->LANDMASK
print("landmask>>>"+min(land)+"  "+max(land))


nnn=dimsizes(fsfc->latitude)-1

fld = new((/1,nlat,nlon/),float)

fname = "./"+outdir+"/snow_"+date_out+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
mtime=0
timem = time
do i=0,ntime-1,nskip_sd
 timem(mtime) = time(i)
 mtime=mtime+1
end do
fbinrecwrite(fname,-1,mtime)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,timem(0:mtime-1))
do i=0,ntime-1,nskip_sd
  fld1= short2flt(fsd->sd(itime1+i,nnn:0:1,:))
  if(max(fld1&longitude).lt.0.) then
    fld1&longitude = fld1&longitude+360.
  end if
  if(i.eq.0) then
   print("ERA5 lon:"+ min(fld1&longitude)+"  "+max(fld1&longitude))
   print("ERA5 lat:"+ min(fld1&latitude)+"  "+max(fld1&latitude))
  end if
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(timed(i)+" original snow depth>>>"+min(fld1)+"  "+max(fld1))
  fld1 = where(fld1.lt.0.,0.,fld1)
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,lon,lat,0)
  fld(0,:,:) = where(land.ne.1,0.,fld(0,:,:))
  fld!0="time"
  fld&time = timed(i)
  fld&time@units = "day"
  fld!1="lat"
  fld!2="lon"
  fld&lat = lat
  fld&lon = lon
  print("interpolated snow depth>>>"+min(fld)+"  "+max(fld))
  print("sbowd:"+min(fld)+"  "+max(fld))
  fbinrecwrite(fname,-1,fld(0,:,:))
  if(netcdf_out) then
   system("rm -f "+outdir_nc+"/snow_"+timed(i)+"_"+grid+"_"+dataset+".nc")
   out = addfile(outdir_nc+"/snow_"+timed(i)+"_"+grid+"_"+dataset+".nc","c")
   out->snow = fld
  end if
end do

fname = "./"+outdir+"/snowt_"+date_out+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
mtime=0
timem = time
do i=0,ntime-1,nskip_sd
 timem(mtime) = time(i)
 mtime=mtime+1
end do
fbinrecwrite(fname,-1,mtime)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,timem(0:mtime-1))
do i=0,ntime-1,nskip_sd
  fld1= short2flt(fst->tsn(itime1+i,nnn:0:1,:))
  if(max(fld1&longitude).lt.0.) then
    fld1&longitude = fld1&longitude+360.
  end if
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(timed(i)+" original snow temperature>>>"+min(fld1)+"  "+max(fld1))
;  fld1 = where(fld1.lt.0.,0.,fld1)
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,lon,lat,0)
;  fld(0,:,:) = where(land.ne.1,0.,fld(0,:,:))
  fld!0="time"
  fld&time = timed(i)
  fld&time@units = "day"
  fld!1="lat"
  fld!2="lon"
  fld&lat = lat
  fld&lon = lon
  print("interpolated snow temperature>>>"+min(fld)+"  "+max(fld))
  print("sbowt:"+min(fld)+"  "+max(fld))
  fbinrecwrite(fname,-1,fld(0,:,:))
  if(netcdf_out) then
   system("rm -f "+outdir_nc+"/snowt_"+timed(i)+"_"+grid+"_"+dataset+".nc")
   out = addfile(outdir_nc+"/snowt_"+timed(i)+"_"+grid+"_"+dataset+".nc","c")
   out->snowt = fld
  end if
end do



