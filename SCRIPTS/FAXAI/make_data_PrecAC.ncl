files = systemfunc("ls -1 ../../OUT_2D/AGNES_1024x1024x72_3km_10s_CONTROL_TKEDISS_PRECDISS_*1972-06-24-00-00-00_0000047520.2D_atm.nc")
print(""+files)
nfiles = dimsizes(files)

wks = gsn_open_wks("ps","ps")

f=addfiles(files,"r")

lat = f[0]->lat
lon = f[0]->lon
print("lon: "+min(lon)+" "+max(lon))
print("lat: "+min(lat)+" "+max(lat))

res = True
res@cnFillOn              = True
res@cnLinesOn             = False
res@cnLineLabelsOn = False                    ; turn the line labels off
res@cnFillMode            = "RasterFill"

lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)

lonmin = -81.1 + 360.
lonmax = -76.3 + 360.
latmin = 36.6
latmax = 39.1
;lonmin = -83. + 360.
;lonmax = -75. + 360.
;latmin = 36.
;latmax = 40.


do i=0,nfiles-1
  fld = f[i]->Precac(0,{latmin:latmax},{lonmin:lonmax})
  time = f[i]->time
  dates = tointeger(cd_calendar(time, 0))
  year = dates(0,0)
  month = sprinti("%0.2d", dates(0,1))
  day = sprinti("%0.2d", dates(0,2))
  mins = sprinti("%0.2d", dates(0,4))
  if(mins.eq.59) then
   dates(0,3) = dates(0,3)+1
  end if
  hour = sprinti("%0.2d", dates(0,3))
  print(dates+"   "+year + "-" + month + "-" + day + "-" + hour)
  res@gsnRightString = year + "-" + month + "-" + day + "-" + hour
  if(i.eq.0) then
    res@gsnLeftString = "ENS 0"
  else
    res@gsnLeftString = "ENS "+i
  end if
  fld = fld/25.4
  print("Precac: "+min(fld)+"  "+max(fld))
end do

