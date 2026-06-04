names = systemfunc("ls -1 ../../OUT_2D/LI_3456x3456x83_100m_3s_SPNUDGE013610_20*.2D_atm.nc")
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
  lons(i) = fbinrecread(names_stations(i),2,1,"float") + 360.
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
  dlat = 1.5*(lat(1)-lat(0))
  dlon = 1.5*(lon(1)-lon(0))
  ttt = f[i]->time
  dates = tointeger(cd_calendar(ttt, 0))
  year = dates(0,0)
  month = dates(0,1)
  day = dates(0,2)
  mins = dates(0,4)
  if(mins.eq.59) then
    dates(0,3) = dates(0,3)+1
    mins = 0
  end if
  hour = dates(0,3)
  time(i) = day+hour/24.+mins/(24.*60.)
  print("time= "+time(i))
  do j=0,ns-1
   land = f[i]->LAND(0,{lats(j)-dlat:lats(j)+dlat},{lons(j)-dlon:lons(j)+dlon})
   print(j+"  "+(lons(j)-360.)+"  "+lats(j))
   print(""+land)
   fld = f[i]->T2m(0,{lats(j)-dlat:lats(j)+dlat},{lons(j)-dlon:lons(j)+dlon})
   fld = where(land.gt.0.,fld,-9999.)
   fld@_FillValue = -9999.
   temp(j,i) = avg(fld)
   fld = f[i]->Q2m(0,{lats(j)-dlat:lats(j)+dlat},{lons(j)-dlon:lons(j)+dlon})
   fld = where(land.gt.0.,fld,-9999.)
   fld@_FillValue = -9999.
   q = 1.e-3*avg(fld)
   rh(j,i) = relhum(temp(j,i),q,100000.)
   u = avg(f[i]->U10m(0,{lats(j)-dlat:lats(j)+dlat},{lons(j)-dlon:lons(j)+dlon}))
   v = avg(f[i]->V10m(0,{lats(j)-dlat:lats(j)+dlat},{lons(j)-dlon:lons(j)+dlon}))
   wdir(j,i) = wind_direction(u,v,0)
   wspeed(j,i) = wind_speed(u,v)
   print(""+names_stations(j)+" "+number(j)+"  "+(temp(j,i)-273.15)+" "+rh(j,i)+" "+wdir(j,i)+" "+wspeed(j,i))
   delete(land)
   delete(fld)
  end do
end do  

do j=0,ns-1
 nameout = "../../OUT_2D/"+sprinti("%03i", j+1)+"_gSAM.bin" 
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

