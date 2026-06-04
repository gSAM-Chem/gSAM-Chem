;              W A R N I N G !;
; This script use the convention that day is year*365+calday of that year.
;
;==========================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
;grid = "360x180_unif"
;grid = "1440x720_dyvar"
;grid = "2304x1152_dyvar"
;grid = "3840x1920_dyvar"
;grid = "2000x500_dyvar"
;grid = "6144x3072_SO"
;grid = "3072x3072_CA"
;grid = "4096x2048_IND"
;grid = "9216x128"
;grid = "2304x2304_dyvar_COMBLE"
;grid = "2048x1024_CONUS"
;grid = "1440x720_GLOBAL"
;grid = "3072x1536_USA"
;grid = "1536x1536_USA"
;grid = "2048x1024_USA"
;grid = "768x768_USA_SM"
;grid = "6144x3072_USA"
;grid = "384x384_EUSA"
grid = "768x768_EUSA"

; IAN
;date=2022092300
; IRMA
;date=2017090500
; CA
;date=2021102112
; IND
;date=2020051500
; CLIMO
;date=2023070100
; COMBLE
;date=2020032800
; CONUS
;date=2019070100
; 
date=2019041300
; EUSA
date=2024021618

indir = "./NC_D"
outdir = "BIN_D"

netcdf_out = True
outdir_nc = "NC_D"

; path to data at NCAR
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc"

;=========================================================

fland=addfile(indir+"/landmask_"+grid+".nc","r")

dataset="era5"

date_start = date
date_end = date
nskip_sd = 1

year = floattointeger(date/1.e6)
month = floattointeger((date-year*1000000)/10000.)
date_obs = (year*1000000+month*10000+100)+"_"+ \
           (year*1000000+month*10000+days_in_month(year,month)*100+23)
year_month = year*100+month

fsd=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_141_sd.ll025sc."+date_obs+".nc","r")
fst=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_238_tsn.ll025sc."+date_obs+".nc","r")

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
  time(i-itime1) = times(i,0)*365.+day_of_year(times(i,0),times(i,1),times(i,2))+times(i,3)/24.
  timed(i-itime1) = dates(i)
  print("  time="+todouble(time(i-itime1))+"  "+timed(i-itime1)+" "+times(i,3))
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
  fld1= fsd->SD(itime1+i,nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(timed(i)+" original snow depth>>>"+min(fld1)+"  "+max(fld1))
  fld1 = where(fld1.lt.0.,0.,fld1)
  if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
  end if
;  fld(i,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),land&lon,land&lat,False)
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,land&lon,land&lat,0)
  fld(0,:,:) = where(land.ne.1,0.,fld(0,:,:))
  fld!0="time"
  fld&time = timed(i)
  fld&time@units = "day"
  fld!1="lat"
  fld!2="lon"
  fld&lon = lon
  fld&lat = lat
  delete(fld@_FillValue)
  print("snowd:"+min(fld)+"  "+max(fld))
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
  fld1= fst->TSN(itime1+i,nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(timed(i)+"  "+timem(i)+" original snow temperature>>>"+min(fld1)+"  "+max(fld1))
  if(any(lon.lt.0.)) then
     fld1 = lonFlip(fld1)
  end if
;  fld1 = where(fld1.lt.0.,0.,fld1)
;  fld(i,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),land&lon,land&lat,False)
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,land&lon,land&lat,0)
;  fld(0,:,:) = where(land.ne.1,0.,fld(0,:,:))
  fld!0="time"
  fld&time = timed(i)
  fld&time@units = "day"
  fld!1="lat"
  fld!2="lon"
  fld&lon = lon
  fld&lat = lat
  delete(fld@_FillValue)
  print("sbowt:"+min(fld)+"  "+max(fld))
  fbinrecwrite(fname,-1,fld(0,:,:))
  if(netcdf_out) then
   system("rm -f "+outdir_nc+"/snowt_"+timed(i)+"_"+grid+"_"+dataset+".nc")
   out = addfile(outdir_nc+"/snowt_"+timed(i)+"_"+grid+"_"+dataset+".nc","c")
   out->snowt = fld
  end if
end do



