load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
grid = "2304x1152_dyvar"
;grid = "4608x2304"
;grid = "9216x4608_dyvar0.5"
;grid = "96x4608_dyvar"
;grid = "576x288_dyvar"
;grid = "512x256_dyvar"
;grid = "8192x4096_dyvar"
;grid = "768x384_dyvar"
;grid = "720x360_dyvar"
;grid = "1440x720"
;;grid = "1440x720"
;grid = "1440x720_dyvar"
;grid = "360x180_unif"
date_start=2005050100
date_end=2005123100

netcdf_out = False

indir = "./DATA"
indir_nc = "./NC_D"
outdir = "./BIN_D"
outdir_nc = "./NC_D"

fland=addfile(indir_nc+"/landmask_"+grid+".nc","r")

dataset="era5_7dayrm"
file_in = indir+"/SFC/sstk_era5_2005_daily_7dayrm.nc"
nskip_sst=1
fsst=addfile(file_in,"r")
file_in = indir+"/SFC/ci_era5_2005_daily_7dayrm.nc"
nskip_ci=1
fci=addfile(file_in,"r")

fsfc = fsst
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
  time(i-itime1) = day_of_year(times(i,0),times(i,1),times(i,2))+times(i,3)/24.
  timed(i-itime1) = dates(i)
  print("  time="+time(i-itime1)+"  "+timed(i-itime1))
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
land = fland->LANDMASK
print("landmask>>>"+min(land)+"  "+max(land))

nnn=dimsizes(fsfc->latitude)-1


fld = new((/1,nlat,nlon/),float)
fname = "./"+outdir+"/sst_"+date+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
mtime=0
timem = time
do i=0,ntime-1,nskip_sst
 timem(mtime) = time(i)
 mtime=mtime+1
end do
fbinrecwrite(fname,-1,mtime)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,timem(0:mtime-1))
do i=0,ntime-1,nskip_sst
  fld1= fsst->sstk(itime1+i,nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(timed(i)+"  original sst>>>"+min(fld1)+"  "+max(fld1))
;  fld(0,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),land&lon,land&lat,False)
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,land&lon,land&lat,0)
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
  fld1= fci->ci(itime1+i,nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  fld1=where(fld1.gt.1.,0.,fld1)
  print(timed(i)+" original sea ice>>>"+min(fld1)+"  "+max(fld1))
;  fld(0,:,:) = area_conserve_remap_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),land&lon,land&lat,False)
  fld(0,:,:) = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,land&lon,land&lat,0)
  fld(0,:,:) = where(land.eq.0.and.fld(0,:,:).gt.0.1,1.,0.)
fld!0="time"
fld&time = timed(i)
fld&time@units = "day"
fld!1="lat"
fld!2="lon"
fld&lon = lon
fld&lat = lat
print("interpolated sea ice>>>"+min(fld)+"  "+max(fld))
fld = where(fld.gt.0.,1.,0.)
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

