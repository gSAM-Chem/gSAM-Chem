; make monthly mean ozone dataset from ERA5 dataset at NCAR

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

netcdf_out = False ; output netcdf file as well?

year_start = 2022
month_start = 8
year_end = 2022
month_end = 10 

years_add = 3 ; add number of years to dataset date 
             ; (in case there is no data for year you want)

;lonmin = 264.
;lonmax = 299.
;latmin = 21.
;latmax = 49.
;identifier = "AGNES"

;lonmin = 122.
;lonmax = 171.
;latmin = 9.
;latmax = 51.
;identifier = "FAXAI"

;lonmin = 260.
;lonmax = 272.
;latmin = 25.
;latmax = 35.
;identifier = "HOUSTON"

;lonmin = 216.
;lonmax = 314.
;latmin = 11.
;latmax = 61.
;identifier = "USA"

;lonmin = 287.
;lonmax = 360.
;latmin = -8.
;latmax = 28.
;identifier = "ITCZ"

;lonmin = 252.
;lonmax = 274.
;latmin = 26.
;latmax = 44.
;identifier = "USA_SM"
;lonmin = 241.
;lonmax = 285.
;latmin = 19.
;latmax = 55.
;identifier = "USA_MID"
;lonmin = 219.
;lonmax = 307.
;latmin = 19.
;latmax = 55.
;identifier = "USA"

;lonmin = 264.
;lonmax = 299.
;latmin = 21.
;latmax = 49.
;identifier = "ISABEL"

;lonmin = 284.
;lonmax = 305.
;latmin = 7.
;latmax = 27.
;identifier = "FIONA"

;lonmin = 219.
;lonmax = 302.
;latmin = 20.
;latmax = 54.
;identifier = "CONUS"

;lonmin = 250.
;lonmax = 340.
;latmin = -5.
;latmax = 40.
;lonmin = 270.
;lonmax = 315.
;latmin = 0.
;latmax = 40.
;identifier = "IRMA_LG"

;lonmin = 310.
;lonmax = 340.
;latmin = -5.
;latmax = 25.
;identifier = "ITCZ"

;lonmin = 250.
;lonmax = 300.
;latmin = 10.
;latmax = 50.
;identifier = "HELENE"

;lonmin = 210.
;lonmax = 310.
;latmin = 15.
;latmax = 60.
;identifier = "USA"

;lonmin = 265.
;lonmax = 305.
;latmin = 19.
;latmax = 51.
;identifier = "EUSA"

;lonmin = 289.
;lonmax = 311.
;latmin = 13.
;latmax = 33.
;identifier = "HUMBERTO"

lonmin = 283.5
lonmax = 288.5
latmin = 39.
latmax = 43.
identifier = "LI"

outdir = "BIN_D"
outdir_nc = "NC_D"

;=================================================================
;=================================================================
;=================================================================
;=================================================================

dataset = "era5_monthly"

time_rec = new((/1000/),float)
date_rec = new((/1000/),integer)
ntime_rec = 0
;==================================================================
do pass = 1,2
do year = year_start,year_end

o3_data = "/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl/"+year+"/e5.moda.an.pl.128_203_o3.ll025sc."+year+"010100_"+year+"120100.nc"
z3_data = "/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl/"+year+"/e5.moda.an.pl.128_129_z.ll025sc."+year+"010100_"+year+"120100.nc"
f_o3 = addfile(o3_data,"r")
f_z3 = addfile(z3_data,"r")

times = calendar_decode2(f_o3->time,-5)
times(:,2) = times(:,2) + 16
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(dates)
ntime = dimsizes(f_o3->time)
dayofyear = day_of_year(times(:,0),times(:,1),times(:,2))+times(:,3)/24.
print("  dayofyear="+dayofyear)
time = (times(:,0)-1+years_add)*365+dayofyear
;time = (times(:,0)-1)*365+dayofyear
print("  time="+time)
print("  type of time="+typeof(time))
nzi = dimsizes(f_o3->level)
grid = identifier
print("data grid:"+grid)


do month = 1,12
if(year.eq.year_start.and.month.lt.month_start) then
  continue
end if 
if(year.eq.year_end.and.month.gt.month_end) then
  break
end if 

mmm = month-1

time_rec(ntime_rec) = time(mmm)
date_rec(ntime_rec) = times(mmm,0)*100+times(mmm,1)
ntime_rec = ntime_rec + 1

; determine pressure levels at z:

nlati = dimsizes(f_o3->latitude)
presi = doubletofloat(f_o3->level(nzi-1:0:1))
z3 = f_z3->Z(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 z3=lonFlip(z3)
end if
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
lat!0 = "lat"
lon!0 = "lon"
lat&lat = lat
lon&lon = lon
z33 = z3(:,{latmin:latmax},{lonmin:lonmax})
nlat = dimsizes(z33&latitude)
nlon = dimsizes(z33&longitude)
wgt = cos(z33&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z33,wgt,0,1))/9.80665

if(pass.eq.1) then
 continue
end if

print(date_rec(ntime_rec-1)+"   reading O3...")
fld3=f_o3->O3(mmm,nzi-1:0:1,nlati-1:0:1,:)
if(lonmin.lt.0..or.lonmax.lt.0) then
 fld3=lonFlip(fld3)
end if
print("o3: "+min(fld3)+"  "+max(fld3))
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
print("done")
if(netcdf_out) then
 fnamenc = outdir_nc+"/ozone_"+dataset+"_"+date_rec(ntime_rec-1)+"_"+grid+".nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->O3 = fld3(:,{latmin:latmax},{lonmin:lonmax})
end if
delete(fld3)

end do ; month
end do ; year
if(pass.eq.1) then
  fname = "./"+outdir+"/ozone_"+dataset+"_"+(date_rec(0)+years_add*100)+"-"+(date_rec(ntime_rec-1)+years_add*100)+"_"+grid+".bin"
  system("rm -f "+fname)
  inds = (/nlon,nlat,nzi/)
  fbinrecwrite(fname,-1,inds)
  fbinrecwrite(fname,-1,ntime_rec)
  fbinrecwrite(fname,-1,time_rec(0:ntime_rec-1))
  ntime_rec = 0
end if
end do ; pass
print(""+time_rec(0:ntime_rec-1))
