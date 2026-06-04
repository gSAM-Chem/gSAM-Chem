load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
date_start=2020052900
date_end=2020070200

grid = "512x512_CP_cart4km"
nlon = 512 
nlat = 512
dx = 4000.
dy = 4000.
;grid = "1944x648_CP_cart4km"
;nlon = 1944
;nlat = 648
;dx = 4000.
;dy = 4000.
;grid = "5184x2592_CP_cart1km"
;nlon = 5184
;nlat = 2592
;dx = 1000.
;dy = 1000.

lon0=215.  ; center of left wall
lat0=6.    ; center of domain

netcdf_out = False

outdir = "BIN_CP"
outdirnc = "NC_CP"

dataset="era5_daily"
nskip_skt=1
fskt=addfile("./BIN_D/SKT/skt_era5_daily_5678.nc","r")

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
  time(i-itime1) = day_of_year(times(i,0),times(i,1),times(i,2))+times(i,3)/24.
  timed(i-itime1) = dates(i)
  print("  time="+time(i-itime1)+"  "+timed(i-itime1))
end do

if(ntime.eq.1) then
  date = ""+date_start
else
  date = ""+date_start+"-"+date_end
end if

pi = acos(-1.)
dlon = dx*cos(lat0*pi/180.)/40000000.*360.
dlat = dy/40000000.*360.
print("dlat= "+dlat)
print("dlon= "+dlon)
lon = new((/nlon/),double)
lat = new((/nlat/),double)
do i=0,nlon-1
 lon(i) = lon0 + i*dlon
end do
do j=0,nlat-1
 lat(j) = lat0 + (j-nlat/2+0.5)*dlat
end do
print("lon: "+min(lon)+" "+max(lon))
print("lat: "+min(lat)+" "+max(lat))
nlon = dimsizes(lon)
nlat = dimsizes(lat)
print("nlat="+nlat)
print("nlon="+nlon)

nnn=dimsizes(fsfc->latitude)-1


fname = "./"+outdir+"/sst_"+date+"_"+grid+"_"+dataset+".bin"
system("rm -f "+fname)
mtime=0
timem = time
do i=0,ntime-1,nskip_skt
 timem(mtime) = time(i)
 mtime=mtime+1
end do
fbinrecwrite(fname,-1,mtime)
fbinrecwrite(fname,-1,nlon)
fbinrecwrite(fname,-1,nlat)
fbinrecwrite(fname,-1,timem(0:mtime-1))
do i=0,ntime-1,nskip_skt
  fld = new((/1,nlat,nlon/),float)
  fld1= fskt->skt(itime1+i,nnn:0:1,:)
  delete(fld1@missing_value)
  delete(fld1@_FillValue)
  print(timed(i)+"  original sst>>>"+min(fld1)+"  "+max(fld1))
  fld = linint2_Wrap(fld1&longitude,fld1&latitude,fld1(:,:),True,lon,lat,0)
  fld =where(fld(0,:,:).gt.306.,306.,fld(0,:,:))
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
  delete(fld)
  delete(fld1)
  if(netcdf_out) then
   system("rm -f "+outdirnc+"/sst_"+timed(i)+"_"+grid+"_"+dataset+".nc")
   out = addfile(outdirnc+"/sst_"+timed(i)+"_"+grid+"_"+dataset+".nc","c")
   out->sst = sst
  end if
end do




