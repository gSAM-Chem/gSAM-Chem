; interpolate SST and Seaice mask datasets available at NCAR
; This script is for running on NCAR's computers only

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
;grid = "3840x1920_dyvar"
;grid = "2000x500_dyvar"
;grid = "1536x1536_AMAZ"
;grid = "3072x3072_SO"
;grid = "3072x3072_CA"
;grid = "4096x2048_IND"
;grid = "9216x128"
;grid = "2304x2304_dyvar_COMBLE"
;grid = "2048x1024_USA"
;grid = "3072x3072_HELENE"
;grid = "3072x3072_MILTON"
;grid = "512x512_MILTON"
;grid = "2048x2048_ITCZ"
;grid = "768x768_HELENE"
;grid = "1024x1024_IRMA"
;grid = "1024x1024_MILTON"
;grid = "1024x1024_HELENE"
;grid = "512x128_TUNE"
;grid = "2048x1024_IRMA"
;grid = "1024x1024_AGNES"
;grid = "2048x2048_AGNES"
grid = "1024x1024_ISABEL"

; IAN
;date_start=2022092300
;date_end=2022093000
; IRMA
;date_start=2017090400
;date_end=2017091200
; IRMA
;date_start=2017083000
;date_end=2017091200
; Katrina
;date_start=2005082500
;date_end=2005082500
; SOCRATES
;date_start=2018021512
;date_end=2018022100
; CA
;date_start=2020051500
;date_end=2020053100
; CA
;date_start=2021102112
;date_end=2021102700
; AMZONIA
;date_start=2014060100
;date_end=2014063023
; CLIMO
;date_start=2023070100
;date_end=2023090100
; CAMBLE
;date_start=2020032800
;date_end=2020040100
; USA
;date_start=2024070100
;date_end=2024070800
; HELENE
;date_start=2024092100
;date_end=2024092900
; ITCZ
;date_start=2024081700
;date_end=2024092000
; MILTON
;date_start=2024100400
;date_end=2024101200
; AGNES
;date_start=1972061700
;date_end=1972062400
; ISABEL
date_start=2003091500
date_end=2003092000

; skipping record (in hours)
nskip_skt=12

netcdf_out = True

indir = "NC_D"
outdir = "BIN_D"
outdir_nc = "NC_D"

delta_file = "/glade/work/marat/CLIMO_PERT/ISABEL/SKT_PW_2015-2024--1994-2003_Sep_DELTA.nc"
;delta_file = "/glade/work/marat/CLIMO_PERT/SKT_PW_ERA5_June_2020-2024--1968-1972.nc"

dataset="era5_DELTA24-03"

; path to data at NCAR
filepath = "/glade/campaign/collections/rda/data/d633000/e5.oper.an.sfc"

;========================================================================

fland=addfile(indir+"/landmask_"+grid+".nc","r")
lat = fland->lat
lon = fland->lon
nlon = dimsizes(lon)
nlat = dimsizes(lat)
print("lat= "+min(lat)+"  "+max(lat))
print("lon= "+min(lon)+" "+max(lon))
print("nlat="+nlat)
print("nlon="+nlon)
land = fland->LANDMASK
print("landmask>>>"+min(land)+"  "+max(land))

; interpolate DELTA SST to model grid:

f_delta = addfile(delta_file,"r")
fld1 = f_delta->skt(::-1,:)
lon1 = fld1&longitude
lat1 = fld1&latitude
if(any(lon1.lt.0.)) then
 fld1&longitude = fld1&longitude + 360.
end if
lon1 = fld1&longitude
lat1 = fld1&latitude
print("lat1= "+min(lat1)+"  "+max(lat1))
print("lon1= "+min(lon1)+" "+max(lon1))
sst_delta = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,land&lon,land&lat,0)
print("sst_delta: "+min(sst_delta)+"  "+max(sst_delta))
delete(fld1)

;-------------------------------------------------------------
; SST

ntimes = 0
dates = new((/2000/),integer)
times = new((/2000/),float)
date = date_start

; determine dates that go into output:

do while (date.le.date_end)

  ntimes = ntimes+1
  year = floattointeger(date/1.e6)
  month = floattointeger((date-year*1000000)/10000.)
  day = floattointeger((date-year*1000000-month*10000)/100.)
  hour = date-year*1000000-month*10000-day*100
  times(ntimes-1) = (year-1)*365.+day_of_year(year,month,day)+hour/24.
  hour = hour+nskip_skt
  if(hour.gt.23) then
   hour = hour-24
   day = day+1
  end if
  if(day.gt.days_in_month(year,month)) then
   month = month+1
   if(month.eq.13) then
     year = year+1
     month = 1
   end if
   day = 1
   hour = 0
  end if
  dates(ntimes-1) = date
  date = year*1000000+month*10000+day*100+hour

end do
print("dates = "+dates(0:ntimes-1)+"    day-of-year: "+todouble(times(0:ntimes-1)))


if(ntimes.eq.1) then
  dates_out = ""+date_start
else
  dates_out = ""+date_start+"-"+date_end
end if



fld = new((/1,nlat,nlon/),float)

fname = "./"+outdir+"/sst_"+dates_out+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
fbinrecwrite(fname,-1,ntimes)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,times(0:ntimes-1))
month_cur = 0
do i=0,ntimes-1
  
  year = floattointeger(dates(i)/1.e6)
  month = floattointeger((dates(i)-year*1000000)/10000.)
  if(month.ne.month_cur) then
   date_obs = (year*1000000+month*10000+100)+"_"+ \
             (year*1000000+month*10000+days_in_month(year,month)*100+23)
   year_month = year*100+month
   fsfc=addfile(filepath+"/"+year_month+"/e5.oper.an.sfc.128_235_skt.ll025sc."+date_obs+".nc","r")
   nnn=dimsizes(fsfc->latitude)-1
   if(isvar("ttt")) then
    delete(ttt)
    delete(dates_obs)
   end if
   ttt = calendar_decode2(fsfc->time,-5)
   dates_obs = ttt(:,0)*1000000+ttt(:,1)*10000+ttt(:,2)*100+ttt(:,3)
   month_curr = month
  end if
  index = ind(dates_obs.eq.dates(i))

  fld1= fsfc->SKT(index,nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(dates(i)+"  original sst>>>"+min(fld1)+"  "+max(fld1))
  if(any(land&lon.lt.0.)) then
    fld1 = lonFlip(fld1)
  end if
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,land&lon,land&lat,0)
  fld(0,:,:) = fld(0,:,:) + sst_delta(:,:)
  fld(0,:,:)=where(land.eq.0.and.fld(0,:,:).gt.306.,306.,fld(0,:,:))
  delete(fld@_FillValue)
  print("interpolated sst>>>"+min(fld)+"  "+max(fld))
  fbinrecwrite(fname,-1,fld(0,:,:))
  if(netcdf_out) then
   sss = fld
   sss!0="time"
   sss&time = times(i)
   sss&time@units = "day"
   sss!1="lat"
   sss!2="lon"
   sss&lon = lon
   sss&lat = lat
   system("rm -f "+outdir_nc+"/sst_"+dates(i)+"_"+grid+"_"+dataset+".nc")
   out = addfile(outdir_nc+"/sst_"+dates(i)+"_"+grid+"_"+dataset+".nc","c")
   filedimdef(out,"time",-1,True)
   out->sst = sss
   delete(sss)
  end if
  delete(fld1)
end do ; i


