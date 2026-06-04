; make monthly mean ozone dataset from ERA5 dataset at NCAR

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"

netcdf_out = True ; output netcdf file as well?

; IRMA
year_start = 2017
month_start = 1
year_end = 2017
month_end = 1
fileroot = "./ERA5/ERA5_OZONE_Sep" ; ozone data
identifier = "GLOBAL" ; will be added to output filename

; CA
;year_start = 2021
;month_start = 9
;year_end = 2021
;month_end = 11
;fileroot = "./ERA5/ERA5_PRES_O3_CA" ; ozone data
;identifier = "CA" ; will be added to output filename

outdir = "BIN_D"
outdir_nc = "NC_D"

;===================================================================
;===================================================================
;===================================================================

dataset = "era5_monthly"

time_rec = new((/1000/),float)
date_rec = new((/1000/),integer)
ntime_rec = 0
;==================================================================
do pass = 1,2
do year = year_start,year_end

o3_data = fileroot+"_"+year+".nc"
z3_data = o3_data
f_o3 = addfile(o3_data,"r")
f_z3 = addfile(z3_data,"r")

times = calendar_decode2(f_o3->time,-5)
times(:,2) = times(:,2) + 16
dates = times(:,0)*1000000+times(:,1)*10000+times(:,2)*100+times(:,3)
print(dates)
ntime = dimsizes(f_o3->time)
dayofyear = day_of_year(times(:,0),times(:,1),times(:,2))+times(:,3)/24.
print("  dayofyear="+dayofyear)
time = (times(:,0)-1)*365+dayofyear
print("  time="+time)
print("  type of time="+typeof(time))
nlat = dimsizes(f_o3->latitude)
nlon = dimsizes(f_o3->longitude)
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

if(pass.eq.1) then
 continue
end if

; determine pressure levels at z:

nlati = dimsizes(f_o3->latitude)
presi = doubletofloat(f_o3->level(nzi-1:0:1))
z3 = short2flt(f_z3->z(mmm,nzi-1:0:1,nlati-1:0:1,:))
lat = doubletofloat(z3&latitude)
lon = doubletofloat(z3&longitude)
wgt = cos(z3&latitude*acos(-1.)/180.)
z3m = dim_avg_Wrap(dim_avg_wgt_n_Wrap(z3,wgt,0,1))/9.80665
print("z3m: "+z3m+"  presi: "+presi)

print(date_rec(ntime_rec-1)+"   reading O3...")
fld3=short2flt(f_o3->o3(mmm,nzi-1:0:1,nlati-1:0:1,:))
print("o3: "+min(fld3)+"  "+max(fld3))
print(presi(:)+"  "+fld3(:,nlat/2,nlon/2))
fld3!1="lat"
fld3!2="lon"
if(max(lon).lt.0.) then
  lon = lon+360.
end if
fld3&lat = lat
fld3&lon = lon
fbinrecwrite(fname,-1,lon)
fbinrecwrite(fname,-1,lat)
fbinrecwrite(fname,-1,z3m)
fbinrecwrite(fname,-1,presi)
do k=0,nzi-1
  fbinrecwrite(fname,-1,fld3(k,:,:))
end do
print("done")
if(netcdf_out) then
 fnamenc = outdir_nc+"/ozone_"+dataset+"_"+date_rec(ntime_rec-1)+"_"+grid+".nc"
 system("rm -f "+fnamenc)
 gg=addfile(fnamenc,"c")
 gg->O3 = fld3
end if
delete(fld3)

end do ; month
end do ; year
if(pass.eq.1) then
  fname = "./"+outdir+"/ozone_"+dataset+"_"+date_rec(0)+"-"+date_rec(ntime_rec-1)+"_"+grid+".bin"
  system("rm -f "+fname)
  inds = (/nlon,nlat,nzi/)
  fbinrecwrite(fname,-1,inds)
  fbinrecwrite(fname,-1,ntime_rec)
  fbinrecwrite(fname,-1,time_rec(0:ntime_rec-1))
  ntime_rec = 0
end if
end do ; pass
print(""+time_rec(0:ntime_rec-1))
