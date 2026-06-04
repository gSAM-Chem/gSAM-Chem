   obs = (/    \
2017083000.,     16.1,    -26.9,      30.,    1008. ,\
2017083006.,     16.2,    -28.3,      35.,    1007. ,\
2017083012.,     16.3,    -29.7,      45.,    1006. ,\
2017083018.,     16.3,    -30.8,      50.,    1004. ,\
2017083100.,     16.3,    -31.7,      55.,     999. ,\
2017083106.,     16.4,    -32.5,      65.,     994. ,\
2017083112.,     16.7,    -33.4,      80.,     983. ,\
2017083118.,     17.1,    -34.2,      95.,     970. ,\
2017090100.,     17.5,    -35.1,     100.,     967. ,\
2017090106.,     17.9,    -36.1,     100.,     967. ,\
2017090112.,     18.4,    -37.3,     100.,     967. ,\
2017090118.,     18.8,    -38.5,     100.,     967. ,\
2017090200.,     19.1,    -39.7,     100.,     967. ,\
2017090206.,     19.1,    -41.1,     100.,     967. ,\
2017090212.,     18.9,    -42.6,      95.,     973. ,\
2017090218.,     18.7,    -44.1,      95.,     973. ,\
2017090300.,     18.5,    -45.5,      95.,     973. ,\
2017090306.,     18.2,    -46.7,      95.,     973. ,\
2017090312.,     17.9,    -47.9,     100.,     969. ,\
2017090318.,     17.6,    -49.2,     100.,     965. ,\
2017090400.,     17.3,    -50.4,     100.,     959. ,\
2017090406.,     17.0,    -51.5,     105.,     952. ,\
2017090412.,     16.8,    -52.6,     110.,     945. ,\
2017090418.,     16.7,    -53.9,     115.,     944. ,\
2017090500.,     16.6,    -55.1,     125.,     943. ,\
2017090506.,     16.6,    -56.4,     135.,     933. ,\
2017090512.,     16.7,    -57.8,     150.,     929. ,\
2017090518.,     16.9,    -59.2,     155.,     926. ,\
2017090600.,     17.3,    -60.6,     155.,     915. ,\
2017090605.,     17.7,    -61.8,     155.,     914. ,\
2017090606.,     17.7,    -61.9,     155.,     914. ,\
2017090611.,     18.1,    -63.1,     155.,     914. ,\
2017090612.,     18.1,    -63.3,     155.,     915. ,\
2017090616.,     18.5,    -64.4,     155.,     915. ,\
2017090618.,     18.6,    -64.7,     150.,     916. ,\
2017090700.,     19.2,    -66.2,     150.,     916. ,\
2017090706.,     19.7,    -67.6,     145.,     920. ,\
2017090712.,     20.2,    -69.0,     145.,     921. ,\
2017090718.,     20.7,    -70.4,     145.,     922. ,\
2017090800.,     21.1,    -71.8,     140.,     919. ,\
2017090805.,     21.5,    -73.0,     135.,     924. ,\
2017090806.,     21.5,    -73.2,     135.,     925. ,\
2017090812.,     21.8,    -74.7,     135.,     927. ,\
2017090818.,     22.0,    -76.0,     140.,     925. ,\
2017090900.,     22.1,    -77.2,     145.,     924. ,\
2017090903.,     22.3,    -77.9,     145.,     924. ,\
2017090906.,     22.4,    -78.3,     130.,     930. ,\
2017090912.,     22.7,    -79.3,     110.,     941. ,\
2017090918.,     23.1,    -80.2,      95.,     938. ,\
2017091000.,     23.4,    -80.9,     100.,     932. ,\
2017091006.,     23.7,    -81.3,     115.,     930. ,\
2017091012.,     24.5,    -81.5,     115.,     931. ,\
2017091013.,     24.7,    -81.5,     115.,     931. ,\
2017091018.,     25.6,    -81.7,     100.,     936. ,\
2017091019.,     25.9,    -81.7,     100.,     936. ,\
2017091100.,     26.8,    -81.7,      80.,     942. ,\
2017091106.,     28.2,    -82.2,      65.,     961. ,\
2017091112.,     29.6,    -82.7,      50.,     970. ,\
2017091118.,     30.9,    -83.5,      45.,     980. ,\
2017091200.,     31.9,    -84.4,      35.,     986. ,\
2017091206.,     32.9,    -85.6,      25.,     997. ,\
2017091212.,     33.8,    -86.9,      20.,    1000. ,\
2017091218.,     34.8,    -88.1,      15.,    1003. ,\
2017091300.,     35.6,    -88.9,      15.,    1004. ,\
2017091306.,     36.2,    -89.5,      15.,    1004. ,\
2017091312.,     36.8,    -90.1,      15.,    1005. /)

   nobs = dimsizes(obs)/5
   times_obs = new((/nobs/),integer)
   lat_obs = new((/nobs/),float)
   lon_obs = new((/nobs/),float)
   p_obs = new((/nobs/),float)
   u_obs = new((/nobs/),float)
   do k=0,nobs-1
     ttt = obs(5*k+0)
     times_obs(k) = floattoint(ttt)
     lat_obs(k) = obs(5*k+1)
     lon_obs(k) = obs(5*k+2)
     u_obs(k) = obs(5*k+3)
     p_obs(k) = obs(5*k+4)
     print((k+1)+"  "+times_obs(k)+"  "+lat_obs(k)+"  "+lon_obs(k)+"  "+u_obs(k)+"  "+p_obs(k))
   end do
   
   times = times_obs

   colors = (/"red","green","blue","brown","magenta"/)
   wks = gsn_open_wks("ps","ps")               ; send graphics to PNG file
   res            = True                          ; map resources
   res@gsnDraw    = False                         ; don't draw
   res@gsnFrame   = False                         ; don't advance frame
   res@vpWidthF   = 0.80                          ; make map bigger
   res@vpHeightF  = 0.80
   res@mpFillOn              = False


   res@tiMainString    = "GSAM v7.1"  ; title
;   res@mpDataBaseVersion     = "HighRes"               ; better map outlines
   res@pmTickMarkDisplayMode = "Always"                  ; better map tickmarks

   pres                  = True               ; polyline resources
   mres                  = True               ; polyline resources
   pres@gsLineThicknessF = 3.0                ; line thickness

   fileroot = "../../OUT_2D/in_"
   nskip = 4
   nens = 10

   do k=0,nens-1

     f=addfile(fileroot+k+".nc","r") 
     lat = f->lat
     lon = f->lon - 360.
     if(k.eq.0) then
       res@mpMaxLatF  = 58. ;max(lat)                          ; select subregion
       res@mpMinLatF  = 8. ;min(lat) 
       res@mpMaxLonF  = -38.;max(lon)
       res@mpMinLonF  = -110. ;min(lon)
       map = gsn_csm_map_ce(wks,res)                         ; create map
     end if
     psfc = f->PSFC
     psfcm = dim_avg(psfc(lat|:,lon|:,time|:))
     pw = f->PW
     u10m = f->U10m
     v10m = f->V10m
     vel10m = sqrt(u10m^2+v10m^2)
     vel10m!1="lat"
     vel10m!2="lon"
     vel10m&lat = lat
     vel10m&lon = lon
     usfc = f->USFC
     vsfc = f->VSFC
     velsfc = sqrt(usfc^2+vsfc^2)
     velsfc!1="lat"
     velsfc!2="lon"
     velsfc&lat = lat
     velsfc&lon = lon
   
     ntimes = dimsizes(f->time)-nskip
     y = new((/ntimes/),float)
     x = y
     if(k.eq.0) then
      xm = x
      ym = y
      pm = y
      um = y
      xm = 0.
      ym = 0.
      pm = 0.
      um = 0.
      time = new((/nens+1,ntimes/),float)
      p1 = new((/nens+1,ntimes/),float)
      u1 = new((/nens+1,ntimes/),float)
     end if
     pmin = x
     umax = x
     umaxs = x
     time(k,:) = f->time(0:ntimes-1)
     do n=0,ntimes-1
     ;  a1D = ndtooned(pw(n,:,:))
     ;  inds  = ind_resolve(maxind(a1D),dimsizes(pw(n,:,:)))
       a1D = ndtooned(psfc(n,:,:)-psfcm(:,:))
       inds  = ind_resolve(minind(a1D),dimsizes(psfc(n,:,:)))
       pmin(n) = min(a1D)+1000.
       x(n) = lon(inds(0,1))
       y(n) = lat(inds(0,0))
       xm(n) = xm(n) + x(n)
       ym(n) = ym(n) + y(n)
       umax(n) = max(vel10m(n,{y(n)-5.:y(n)+5.},{x(n)-5.:x(n)+5.}))*1.94384
       umaxs(n) = max(velsfc(n,{y(n)-5.:y(n)+5.},{x(n)-5.:x(n)+5.}))*1.94384
       pm(n) = pm(n) + pmin(n)
       um(n) = um(n) + umaxs(n)
     end do
     
     p1(k,:) = pmin(:) 
     u1(k,:) = umaxs(:) 
     print("ens="+k)
     print(time(k,:)+"  "+x+"  "+y+"  "+pmin+"  "+umax+"  "+umaxs)
     pres@gsLineColor = colors(0) 
     gsn_polyline(wks,map,x,y,pres)      ; draw the traj
;     mres@gsMarkerIndex  = k           ; marker
;     gsn_polymarker(wks,map,x(::48),y(::48),mres)      ; draw the traj

     delete(psfc)
     delete(pmin)
     delete(a1D)
     delete(y)
     delete(x)
   pres@gsLineColor = "black"
   pres@gsLineThicknessF = 3.0                ; line thickness
   gsn_polyline(wks,map,lon_obs,lat_obs,pres) 
   end do
   xm = xm/nens
   ym = ym/nens
   pm = pm/nens
   um = um/nens
   print(xm+"  "+ym+"  "+pm+"  "+um)
   pres@gsLineColor = colors(1) 
   gsn_polyline(wks,map,xm,ym,pres)      ; draw the traj
   delete(xm)
   delete(ym)
   delete(pm)
   delete(um)
   delete(lat)
   delete(lon)
   

   f=addfile("/glade/p/unsb0007/SAM7.1_TC/OUT_2D/in_4km_3h.nc","r")
   lat = f->lat
   lon = f->lon - 360.
   psfcd = f->PSFC
   psfcmd = dim_avg(psfcd(lat|:,lon|:,time|:))
     pwd = f->PW
     u10md = f->U10m
     v10md = f->V10m
     vel10md = sqrt(u10md^2+v10md^2)
     vel10md!1="lat"
     vel10md!2="lon"
     vel10md&lat = lat
     vel10md&lon = lon
     usfcd = f->USFC
     vsfcd = f->VSFC
     velsfcd = sqrt(usfcd^2+vsfcd^2)
     velsfcd!1="lat"
     velsfcd!2="lon"
     velsfcd&lat = lat
     velsfcd&lon = lon
     nskip = 3
     ntimes = dimsizes(f->time)-nskip
     y = new((/ntimes/),float)
     x = y
     xm = x
     ym = y
     pm = y
     um = y
     xm = 0.
     ym = 0.
     pm = 0.
     um = 0.
     timed = new((/ntimes/),float)
     p1d = new((/ntimes/),float)
     u1d = new((/ntimes/),float)
     pmind = x
     umaxd = x
     umaxsd = x
     timed(:) = f->time(0:ntimes-1)
     do n=0,ntimes-1
       a1Dd = ndtooned(psfcd(n,:,:)-psfcmd(:,:))
       indsd  = ind_resolve(minind(a1Dd),dimsizes(psfcd(n,:,:)))
       pmind(n) = min(a1Dd)+1000.
       x(n) = lon(indsd(0,1))
       y(n) = lat(indsd(0,0))
       umaxd(n) = max(vel10md(n,{y(n)-5.:y(n)+5.},{x(n)-5.:x(n)+5.}))*1.94384
       umaxsd(n) = max(velsfcd(n,{y(n)-5.:y(n)+5.},{x(n)-5.:x(n)+5.}))*1.94384
     end do
     p1(nens,:) = pmind(:) 
     u1(nens,:) = umaxsd(:) 
     time(nens,:) = timed(:)

     print(timed(:)+"  "+x+"  "+y+"  "+pmind+"  "+umaxd+"  "+umaxsd)
     pres@gsLineColor = colors(2)
     gsn_polyline(wks,map,x,y,pres)      ; draw the traj


   
   pres@gsLineColor = "black"
   pres@gsLineThicknessF = 3.0                ; line thickness
   gsn_polyline(wks,map,lon_obs,lat_obs,pres) 


   draw(map)                                           
   frame(wks)                                         
  
   delete(res)
   res = True
   delete(colors)
   colors=new((/nens+1/),string)
   colors(0:nens-1) = "green"
   colors(nens) = "blue"
   res@xyLineColors = colors
   res@xyDashPattern  = 0
   res@xyLineThicknessF = 3.
   plot  = gsn_xy (wks,time,p1,res)
   plot  = gsn_xy (wks,time,u1,res)


