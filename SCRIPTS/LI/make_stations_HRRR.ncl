names = systemfunc("ls -1 /glade/work/marat/LI_DATA/HRRR/hrrr_201807*/hrrr_ll_*.nc")
doi = 190. ; initial day of year


print(""+names)
nt = dimsizes(names)
print("nt = "+nt)

names_stations = systemfunc("ls -1 /glade/work/marat/LI_STATIONS/isd_bin/*.bin") 


ns = dimsizes(names_stations)
print("ns = "+ns)
number = new((/ns/),integer)
lats = new((/ns/),float)
lons = new((/ns/),float)
do i=0,ns-1
  number(i) = fbinrecread(names_stations(i),0,1,"integer")
  lats(i) = fbinrecread(names_stations(i),1,1,"float")
  lons(i) = fbinrecread(names_stations(i),2,1,"float") 
  print(i+"  "+names_stations(i)+" "+lats(i)+" "+lons(i))
end do

temp = new((/ns,nt/),float)
rh = new((/ns,nt/),float)
wdir = new((/ns,nt/),float)
wspeed = new((/ns,nt/),float)
time = new((/nt/),float)

f=addfiles(names,"r")

do i=0,nt-1
  lat = f[i]->lat
  lon = f[i]->lon
  dlat = 0.5*(lat(1)-lat(0))
  dlon = 0.5*(lon(1)-lon(0))
  time(i) = doi
  doi = doi + 1./24. 
  print("time= "+time(i))
  do j=0,ns-1
   print(j+"  "+lons(j)+"  "+lats(j))
   fld = f[i]->TMP_P0_L103_GLC0({lats(j)-dlat:lats(j)+dlat},{lons(j)-dlon:lons(j)+dlon})
   temp(j,i) = avg(fld)
   fld = f[i]->DPT_P0_L103_GLC0({lats(j)-dlat:lats(j)+dlat},{lons(j)-dlon:lons(j)+dlon})
   dpt = avg(fld)
   rh(j,i) = relhum_ttd(temp(j,i),dpt,0)
   u = 0.
   v = 0.
   wdir(j,i) = wind_direction(u,v,0)
   wspeed(j,i) = wind_speed(u,v)
   print(""+names_stations(j)+" "+number(j)+"  "+(temp(j,i)-273.15)+" "+rh(j,i)+" "+wdir(j,i)+" "+wspeed(j,i))
   delete(fld)
  end do
end do  

do j=0,ns-1
 nameout = "../../OUT_2D/"+sprinti("%03i", j+1)+"_HRRR.bin" 
 print(""+nameout)
 fbinrecwrite(nameout,0,number(j))
 fbinrecwrite(nameout,1,lat(j))
 fbinrecwrite(nameout,2,lon(j))
 fbinrecwrite(nameout,3,nt)
 fbinrecwrite(nameout,4,time(:))
 fbinrecwrite(nameout,5,temp(j,:)-273.15)
 fbinrecwrite(nameout,6,rh(j,:))
 fbinrecwrite(nameout,7,wdir(j,:))
 fbinrecwrite(nameout,8,wspeed(j,:))
end do

