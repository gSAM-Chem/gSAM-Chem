ensfile = "ens_pmin_vmax.bin"

title = ""

cases = (/" NN " /)

month_obs = 9

p_min = 890.
p_max = 1010.
u_min = 0.
u_max = 90.


lat_min = 10.
lat_max = 50.
lon_min= 122.
lon_max= 171.

leg_lon = 160. ; legend starting lon (left side)
dleg_lon = 2.
leg_lat = 15. ; legend starting lon (left side)
dleg_lat = 1.


day_start = 5.
day_end = 10.

nstart = 47 ; start index of observed track data
nend = 89 ; end of index of observed track data


;========================================
nstart = nstart-1
nend = nend-1


fobs = addfile("/glade/work/marat/IRMA/BestTrack_FAXAI.nc","r")
times_obs = fobs->date_time
print("times_obs= "+times_obs)
lat_obs = tofloat(fobs->lat)
lon_obs = tofloat(fobs->lon)
p_obs = tofloat(fobs->pres)
u_obs = tofloat(fobs->speed)/1.94384
nobs = dimsizes(u_obs)
;nend = nobs-1 ; end of index of observed track data
gust_obs = tofloat(fobs->gust)/1.94384
date_string = chartostring(times_obs)
year  = toint(str_get_field(date_string, 1, "-"))
month = toint(str_get_field(date_string, 2, "-"))
day   = toint(str_get_field(date_string, 3, "-"))
doy = day_of_year(year, month, day)
; Extract time part (HH:MM:SS)
time_string = str_get_field(date_string, 2, " ")
; Extract hour as integer
hour = toint(str_get_field(time_string, 1, ":"))
day_obs = doy+hour/24.
lat_obs0 = day_obs
lat_obs0 = -9999.
lat_obs0@_FillValue = -9999.
lon_obs0 = lat_obs
do k=0,nobs-1
  if(hour(k).eq.0.) then
    lat_obs0(k) = lat_obs(k)
    lon_obs0(k) = lon_obs(k)
  end if
  print((k+1)+"  "+times_obs(k,:)+"  "+day_obs(k)+"  "+lat_obs(k)+"  "+lon_obs(k)+"  "+u_obs(k)+"  "+gust_obs(k)+"  "+p_obs(k))
end do
day0 = day_of_year(year(0), month_obs, 1)
print("day0 = "+day0)
delete(day)

nens = fbinrecread(ensfile,0,1,"integer")
print("nens="+nens)
nrec = 0

do i=0,nens-1
  nrec = nrec+1
  ntime = fbinrecread(ensfile,nrec,1,"integer")
  print("ntime="+ntime)
  if(i.eq.0) then
   pmin = new((/nens,1000/),float)
   pmin = -999.
   pmin@_FillValue = -999.
   vmax = pmin
   day = pmin
   lat = pmin
   lon = pmin
   land = pmin
  end if
  nrec = nrec+1
  day(i,:ntime-1) = fbinrecread(ensfile,nrec,ntime,"float")
  nrec = nrec+1
  pmin(i,:ntime-1) = fbinrecread(ensfile,nrec,ntime,"float")
  nrec = nrec+1
  vmax(i,:ntime-1) = fbinrecread(ensfile,nrec,ntime,"float")
  nrec = nrec+1
  lon(i,:ntime-1) = fbinrecread(ensfile,nrec,ntime,"float")
  nrec = nrec+1
  lat(i,:ntime-1) = fbinrecread(ensfile,nrec,ntime,"float")
  nrec = nrec+1
  land(i,:ntime-1) = fbinrecread(ensfile,nrec,ntime,"float")
  print("day:"+min(day)+"  "+max(day))
  print("pmin:"+min(pmin)+"  "+max(pmin))
  print("vmax:"+min(vmax)+"  "+max(vmax))
  print("lon:"+min(lon)+"  "+max(lon))
  print("lat:"+min(lat)+"  "+max(lat))
  print("land:"+min(land)+"  "+max(land))
end do

 wks  = gsn_open_wks("ps","ps")             ; send graphics to PNG file
 
 res                     = True                   ; plot mods desired
 res@gsnMaximize         = True                   ; maximize plot in frame
 markers          = (/0,2,3,5,6,7,8,9,12,13,14/)                     ; choose type of marker 
 colors = (/"red", "blue", "green", "orange", "purple", "cyan", "magenta", "darkgreen", "brown","gray","black"/)
; colors = (/"green","palegreen","purple3",\
;  "salmon3", "royalblue2","pink3","red2","sienna2","chocolate2","cyan2",\
;  "violet","darkviolet","cadetblue4","blue"/)
 res@xyMarkers          = markers                     ; choose type of marker 
 res@xyMarkerColors     = colors                  ; Marker color
 res@xyMarkerSizeF       = 0.006                  ; Marker size (default 0.01)
 res@xyDashPatterns      = 0                      ; solid line 
 res@xyLineThicknesses   = (/3,3/)                ; set second line to 2
 res@vpWidthF            = 0.5           ; change aspect ratio of plot
 res@vpHeightF           = 0.3

 res@tiMainString  = title

 res@trXMinF = 0.
 res@trXMaxF = 90.
 res@trYMinF = 890.
 res@trYMaxF = 1020.

 res@gsnDraw  = False
 res@gsnFrame = False

day = day - day0 + 1
day_obs = day_obs - day0 + 1


iday = dimsizes(day(0,:))-1
do i=0,dimsizes(day(0,:))-1
 if(day(0,i).le.day_end) then
   print(i+"    "+day(0,i))
   iday = i
 end if
end do

res@xyMarkLineMode     = "Markers" 
plot  = gsn_csm_xy (wks,vmax(:,:iday),pmin(:,:iday),res)        ; create plot
pres                  = True
pres@gsMarkerIndex    = 16
pres@gsMarkerSizeF    = 0.006
pres@gsMarkerColor    = "black"
ppp1 = gsn_add_polymarker(wks, plot, u_obs(nstart:nend), p_obs(nstart:nend), pres)
;pres@gsMarkerColor    = "green"
;ppp2 = gsn_add_polymarker(wks, plot, gust_obs(nstart:nend), p_obs(nstart:nend), pres)
draw(plot)
frame(wks)

 res@trXMinF = day_start
 res@trXMaxF = day_end
 res@trYMinF = 0
 res@trYMaxF = 90


res@xyMarkLineMode     = "Markers"  ; choose which have markers

pres@gsMarkerSizeF    = 0.005
res@gsLineThicknessF = 3.0                ; line thickness
res@xyLineColors = "blue"
nrunning = 5
wgt=new((/nrunning/),float)
wgt = 1./nrunning 
vmax = wgt_runave(vmax,wgt,0)
delete(wgt)
vmax@_FillValue = -9999.
vmax = where(vmax.lt.-9000.,-9999.,vmax)
day = where(vmax.lt.-9000.,-9999.,day)
day@_FillValue = -9999.
plot  = gsn_csm_xy (wks,day(:,:iday),vmax(:,:iday),res)        ; create plot
res@gsLineColor = "black"
gsn_polyline(wks,plot,day_obs(nstart:nend),u_obs(nstart:nend),res)
draw(plot)
frame(wks)

res@trYMinF = 890
res@trYMaxF = 1010

res@gsLineColor = "blue"
nrunning = 5
wgt=new((/nrunning/),float)
wgt = 1./nrunning 
pmin = wgt_runave(pmin,wgt,0)
delete(wgt)
plot  = gsn_csm_xy (wks,day(:,:iday),pmin(:,:iday),res)        ; create plot
res@gsLineColor = "black"
gsn_polyline(wks,plot,day_obs(nstart:nend),p_obs(nstart:nend),res)
draw(plot)
frame(wks)

mres=True
mres@mpMaxLatF  = lat_max                          ; select subregion
mres@mpMinLatF  = lat_min
mres@mpMaxLonF  = lon_max
mres@mpMinLonF  = lon_min
mres@gsnDraw  = False
mres@gsnFrame = False
mres@tiMainString  = title

mres@txFontHeightF = 0.01

mres@mpMaskAreaSpecifiers = (/"Land"/)  ; Try disabling land masking

mres@mpOutlineOn           = True
mres@mpOutlineBoundarySets = "AllBoundaries"
mres@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
mres@mpDataSetName = "Earth..1"

plot = gsn_csm_map_ce(wks,mres)                         ; create map

nrunning = 5
wgt=new((/nrunning/),float)
wgt = 1./nrunning 
lon = wgt_runave(lon,wgt,0)
lat = wgt_runave(lat,wgt,0)
delete(wgt)

lon = lon - 360.
mres@gsLineThicknessF = 2.0                ; line thickness
mres@gsMarkerSizeF = 0.005           ; marker
mres@gsMarkerColor = "blue"           ; marker
legx = new((/2/),float)
legy = new((/2/),float)

do i=0,nens-1
  mres@gsMarkerIndex =  markers(i)
  mres@gsMarkerColor =  colors(i)
  gsn_polymarker(wks,plot,lon(i,:iday:12),lat(i,:iday:12),mres)      ; draw the traj
  mres@gsLineColor = colors(i)
  gsn_polyline(wks,plot,lon(i,:iday),lat(i,:iday),mres)      ; draw the traj
  legx(0) = leg_lon 
  legx(1) = leg_lon+dleg_lon 
  legy(0) = leg_lat-i*dleg_lat 
  legy(1) = leg_lat-i*dleg_lat 
  gsn_polyline(wks,plot,legx,legy,mres)
  gsn_polymarker(wks,plot,legx(0)+0.5*dleg_lon,legy(0),mres)      ; draw the traj
  mres@txJust = 1
  dum1 = gsn_add_text(wks,plot,cases(i),legx(1),legy(1),mres)
end do

mres@gsLineThicknessF = 3.0                ; line thickness
mres@gsLineColor = "black"
gsn_polyline(wks,plot,lon_obs(nstart:nend),lat_obs(nstart:nend),mres)
mres@gsMarkerColor = "black"           ; marker
;gsn_polymarker(wks,plot,lon_obs0(nstart:),lat_obs0(nstart:),mres)      ; draw the traj
i=nens
legx(0) = leg_lon 
legx(1) = leg_lon+dleg_lon 
legy(0) = leg_lat-i*dleg_lat 
legy(1) = leg_lat-i*dleg_lat 
gsn_polyline(wks,plot,legx,legy,mres)
mres@txJust = 1
dum2 = gsn_add_text(wks,plot," BestTrack",legx(1),legy(1),mres)


draw(plot)
frame(wks)

