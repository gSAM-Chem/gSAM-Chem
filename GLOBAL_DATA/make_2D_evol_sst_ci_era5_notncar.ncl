load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
;grid = "2304x1152_dyvar"
;grid = "4608x2304"
;grid = "9216x4608_uniform"
;grid = "9216x4608_dyvar"
;grid = "576x288_dyvar"
;grid = "512x256_dyvar"
;grid = "8192x4096_dyvar"
;grid = "2304x1152_dyvar"
;grid = "1440x720_dyvar"
grid = "3840x1920_dyvar"
;grid = "1536x1536_AMAZ"
;grid = "3072x3072_SO"
;grid = "3072x3072_CA"
;grid = "9216x128"

; IAN
;date_start=2022092300
;date_end=2022093000
; IRMA
date_start=2017090500
date_end=2017091200
; Katrina
;date_start=2005082500
;date_end=2005082500
; SOCRATES
;date_start=2018021512
;date_end=2018022100
; CA
;date_start=2021102112
;date_end=2021102700
; AMZONIA
;date_start=2014060100
;date_end=2014063023

seaice = True ; read the sea ice data or not

netcdf_out = True

; ERA5 data file
;filedata = "./ERA5/ERA5_SFC_20-27OCT_2021_CA.nc"
filedata = "./ERA5/ERA5_SST_CI_Sep5-12_2017_6hourly.nc"

nskip_skt=1
nskip_ci=4
;===============================================================
;===============================================================
;===============================================================

indir = "NC_D"
outdir = "BIN_D"
outdir_nc = "NC_D"

fland=addfile(indir+"/landmask_"+grid+".nc","r")

dataset="era5"
fskt=addfile(filedata,"r")
fci=addfile(filedata,"r")

fsfc = fskt
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
  time(i-itime1) = (times(i,0)-1)*365.+day_of_year(times(i,0),times(i,1),times(i,2))+times(i,3)/24.
  timed(i-itime1) = dates(i)
  print("  time="+todouble(time(i-itime1))+"  "+timed(i-itime1))
end do

if(ntime.eq.1) then
  date = ""+date_start
else
  date = ""+date_start+"-"+date_end
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
fname = "./"+outdir+"/sst_"+date+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
mtime=0
timem = time
do i=0,ntime-1,nskip_skt
 timem(mtime) = time(i)
 mtime=mtime+1
end do
print("mtime="+mtime)
fbinrecwrite(fname,-1,mtime)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,timem(0:mtime-1))
do i=0,ntime-1,nskip_skt
  fld1= short2flt(fskt->skt(itime1+i,nnn:0:1,:))
  if(max(fld1&longitude).lt.0.) then
    fld1&longitude = fld1&longitude+360.
  end if
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
 if(i.eq.0) then
   print("ERA5 lon:"+ min(fld1&longitude)+"  "+max(fld1&longitude))
   print("ERA5 lat:"+ min(fld1&latitude)+"  "+max(fld1&latitude))
 end if
  print(timed(i)+"  original sst>>>"+min(fld1)+"  "+max(fld1))
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,lon,lat,0)
  fld(0,:,:)=where(land.eq.0.and.fld(0,:,:).gt.306.,306.,fld(0,:,:))
fld!0="time"
fld&time = timed(i)
fld&time@units = "day"
fld!1="lat"
fld!2="lon"
fld&lon = lon
fld&lat = lat
print("interpolate sst>>>"+min(fld)+"  "+max(fld))
print("sst:"+min(fld)+"  "+max(fld))
fbinrecwrite(fname,-1,fld(0,:,:))
sst = fld
delete(fld1)
if(netcdf_out) then
 system("rm -f "+outdir_nc+"/sst_"+timed(i)+"_"+grid+"_"+dataset+".nc")
 out = addfile(outdir_nc+"/sst_"+timed(i)+"_"+grid+"_"+dataset+".nc","c")
 out->sst = sst
end if
end do

if(.not.seaice) then
  exit
end if

fname = "./"+outdir+"/icemask_"+date+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
mtime=0
timem = time
do i=0,ntime-1,nskip_ci
 timem(mtime) = time(i)
 mtime=mtime+1
end do
fbinrecwrite(fname,-1,mtime)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,timem(0:mtime-1))
do i=0,ntime-1,nskip_ci
  fld1= short2flt(fci->siconc(itime1+i,nnn:0:1,:))
  if(max(fld1&longitude).lt.0.) then
    fld1&longitude = fld1&longitude+360.
  end if
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  fld1=where(fld1.gt.1.,0.,fld1)
  print(timed(i)+" original sea ice>>>"+min(fld1)+"  "+max(fld1))
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,lon,lat,0)
  fld(0,:,:) = where(land.eq.0.and.fld(0,:,:).gt.0.1,1.,0.)
fld!0="time"
fld&time = timed(i)
fld&time@units = "day"
fld!1="lat"
fld!2="lon"
fld&lat = lat
fld&lon = lon
print("interpolated sea ice>>>"+min(fld)+"  "+max(fld))
fld@long_name = "Ice Mask: 1-ice, 0-no ice"
print("icemask:"+min(fld)+"  "+max(fld))
fbinrecwrite(fname,-1,fld(0,:,:))
if(netcdf_out) then
 system("rm -f "+outdir_nc+"/seaice_"+timed(i)+"_"+grid+"_"+dataset+".nc")
 out = addfile(outdir_nc+"/seaice_"+timed(i)+"_"+grid+"_"+dataset+".nc","c")
 out->icemask = fld
end if
delete(fld1)
end do



