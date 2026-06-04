print("input file: "+filename_in)
print("output file: "+filename_out)

fi = addfile(filename_in,"r")
system("rm -f "+filename_out)
fo = addfile(filename_out,"c")

filedimdef(fo,"time",-1,True)

loni = fi->lon
lati = fi->lat
nloni=dimsizes(fi->lon)
nlati=dimsizes(fi->lat)
print("time="+fi->time)
print("nloni="+nloni)
print("nlati="+nlati)
; uniform grid to interpolate to:

nlat = 180
nlon = 360

dlat = (lati(nlati-1) - lati(0))/nlat
dlon = (loni(nloni-1) - loni(0))/nlon
print("nlon="+nlon)
print("nlat="+nlat)
print("dlon="+dlon)
print("dlat="+dlat)

lon = new((/nlon/),float)
lat = new((/nlat/),float)
do i=0,nlon-1
  lon(i) = loni(0)+dlon*i
end do
do i=0,nlat-1
  lat(i) = lati(0)+dlat*i
end do
lon@units = loni@units
lat@units = lati@units

vNames = getfilevarnames(fi)
wgt = cos(lati*acos(-1.)/180.)

do n=0,dimsizes(vNames)-1
 v = fi->$vNames(n)$
 nd=dimsizes(dimsizes(v))
 if(nd.eq.4) then
   print(""+vNames(n))
   fm = area_hi2lores(loni,lati,v,True,wgt,lon,lat,False)
   fm!0="time"
   fm!1="z"
   fm!2="lat"
   fm!3="lon"
   fm&time = fi->time
   fm&lat=lat
   fm&lon=lon
   fo->$vNames(n)$ = fm(:,:,{-60.:60.},:)
 end if
 delete(v)
end do
fo->p = fi->p


