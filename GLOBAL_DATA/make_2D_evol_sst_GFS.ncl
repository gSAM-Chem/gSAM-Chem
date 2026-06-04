; interpolate SST and Seaice mask datasets available at NCAR
; This script is for running on NCAR's computers only

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

grid = "384x384_EUSA"

date_start=2024021618
date_end=2024021718

; skipping record (in hours)
nskip_skt=6

netcdf_out = True

indir = "NC_D"
outdir = "BIN_D"
outdir_nc = "NC_D"

dataset="GFS"


;========================================================================

fland=addfile(indir+"/landmask_"+grid+".nc","r")
lat = fland->lat
lon = fland->lon
nlon = dimsizes(lon)
nlat = dimsizes(lat)
print("nlat="+nlat)
print("nlon="+nlon)
land = fland->LANDMASK
print("landmask>>>"+min(land)+"  "+max(land))

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
;make last time a bit longer (just in case)
times(ntimes-1) = times(ntimes-1)+0.03
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
  
  fobs = addfile("/glade/work/marat/GFS/gfs.0p25."+dates(i)+".f000.grib2","r")  
  nnn=dimsizes(fobs->lat_0)-1
  fld1= fobs->TMP_P0_L1_GLL0(nnn:0:1,:)
  land1= fobs->LAND_P0_L1_GLL0(nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  fld1 = where(land1.gt.0.1,-999.,fld1)
  fld1@_FillValue = -999.
  poisson_grid_fill(fld1, True, 1, 2000, 1.e-2, 0.6, 0)
  print(dates(i)+"  original sst>>>"+min(fld1)+"  "+max(fld1))
  if(any(land&lon.lt.0.)) then
    fld1 = lonFlip(fld1)
  end if
  fld(0,:,:) = linint2_Wrap(fld1&lon_0,fld1&lat_0,fld1(:,:),True,land&lon,land&lat,0)
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


