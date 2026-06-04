filename = "../../OUT_2DZ/MJO_2048x512x38_20s_dyvar_SFLAT300K_THOM_1024_AVG.com2DZ_1.nc"
filename_out = "../../OUT_2DZ/MJO_2048x512x38_20s_dyvar_SFLAT300K_THOM_1024_AVG.com2DZ_1deg.nc"

file_lats = "./../../GRIDS/lat_512x2048.txt"

nlat = 179
lat = new((/nlat/),float)
do i=0,nlat-1
  lat(i) = -89.+1.*i
end do
lat@units ="degrees north"

fi = addfile(filename,"r")
system("rm -f "+filename_out)
fo = addfile(filename_out,"c")
vNames = getfilevarnames(fi)

z = fi->z
lev = fi->lev
nlati = dimsizes(fi->lat)
print("nlat="+nlat)
lati = asciiread(file_lats,(/nlati/),"float")
print("lati="+lati)
nz = dimsizes(z)
fm = new((/1,nz,nlat/),float)
do n=0,dimsizes(vNames)-1
 v = fi->$vNames(n)$
 nd=dimsizes(dimsizes(v))
 if(nd.eq.3) then
   print(""+vNames(n))
   v&lat = lati
   do i=0,nlat-1
     fm(:,:,i) = dim_avg(v(:,:,{lat(i)-0.5:lat(i)+0.5}))
   end do
   fm!0="time"
   fm!1="lev"
   fm!2="lat"
   fm&lev=lev
   fm&lat=lat
   fo->$vNames(n)$ = fm(:,:,:)
 end if
 delete(v)
end do
;fo->p = p

