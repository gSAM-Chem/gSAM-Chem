stations_obs = systemfunc("ls -1 /glade/work/marat/LI_STATIONS/isd_bin/0*.bin")

stations_SAM = systemfunc("ls -1 ../../OUT_2D/0*gSAM.bin")
stations_HRRR = systemfunc("ls -1 ../../OUT_2D/0*HRRR.bin")

ns = dimsizes(stations_obs)
print("ns = "+ns)
ns1 = dimsizes(stations_SAM)
print("ns gSAM = "+ns1)
ns2 = dimsizes(stations_HRRR)
print("ns HRRR = "+ns2)

if(ns.ne. ns1) then
  print("number of stations in obs for gSAM should be the same!"+"")
  exit
end if
if(ns.ne. ns2) then
  print("number of stations in obs for HRRR should be the same!"+"")
  exit
end if

wks = gsn_open_wks("ps","ps")

plot =new((/ns/),graphic)
dumplot =new((/ns/),graphic)
poly = new(ns,graphic)
dum = new(ns,graphic)
plot1 =new((/ns/),graphic)
dumplot1 =new((/ns/),graphic)
poly1 = new(ns,graphic)
dum1 = new(ns,graphic)
plot2 =new((/ns/),graphic)
poly2 = new(ns,graphic)
dum2 = new(ns,graphic)
plot3 =new((/ns/),graphic)
poly3 = new(ns,graphic)
dum3 = new(ns,graphic)

res  = True

res@gsnDraw     = False
res@gsnFrame    = False

res@gsLineThicknessF  = 3.0 
res@xyLineThicknessF  = 3.0 
res@xyLineColors       = (/"red","blue","red","magenta","green"/)

resl = True
resl@gsnMaximize = True

polyres                   = True
polyres@gsMarkerIndex     = 16          ; polymarker style
polyres@gsMarkerSizeF     = 5.          ; polymarker size

  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.04                     ; font smaller. default big

 res@trXMinF = 0.
 res@trXMaxF = 30.

res@tmYLLabelFontHeightF = 0.03
res@tmXBLabelFontHeightF = 0.03

iplot = 0
do i=0,ns-1
  is = fbinrecread(stations_obs(i),0,1,"integer")
  lats = fbinrecread(stations_obs(i),1,1,"float")
  lons = fbinrecread(stations_obs(i),2,1,"float")
  nt_obs = fbinrecread(stations_obs(i),3,1,"integer")
  hour_obs = fbinrecread(stations_obs(i),4,(/nt_obs/),"float")
  hour_obs = (hour_obs - tofloat(tointeger(hour_obs(0))))*24.
  hour_obs = hour_obs - 5. ; LST
  print("hour_obs:"+hour_obs)
  t_obs = fbinrecread(stations_obs(i),5,(/nt_obs/),"float")
  rh_obs = fbinrecread(stations_obs(i),6,(/nt_obs/),"float")
  wdir_obs = fbinrecread(stations_obs(i),7,(/nt_obs/),"float")
  wspeed_obs = fbinrecread(stations_obs(i),8,(/nt_obs/),"float")
  print("wspeed_obs:"+wspeed_obs)
  nt_SAM = fbinrecread(stations_SAM(i),3,1,"integer")
  hour_SAM = fbinrecread(stations_SAM(i),4,(/nt_SAM/),"float")
  hour_SAM = (hour_SAM - tofloat(tointeger(hour_SAM(0))))*24.
  hour_SAM = hour_SAM - 5. ; LST
  print("hour_SAM:"+hour_SAM)
  t_SAM = fbinrecread(stations_SAM(i),5,(/nt_SAM/),"float")
;  t_SAM = 9./5.*t_SAM + 32.
  rh_SAM = fbinrecread(stations_SAM(i),6,(/nt_SAM/),"float")
  wdir_SAM = fbinrecread(stations_SAM(i),7,(/nt_SAM/),"float")
  wspeed_SAM = fbinrecread(stations_SAM(i),8,(/nt_SAM/),"float")
  nt_HRRR = fbinrecread(stations_HRRR(i),3,1,"integer")
  hour_HRRR = fbinrecread(stations_HRRR(i),4,(/nt_HRRR/),"float")
  hour_HRRR = (hour_HRRR - tofloat(tointeger(hour_HRRR(0))))*24.
  hour_HRRR = hour_HRRR - 5. ; LST
  print("hour_HRRR:"+hour_HRRR)
  t_HRRR = fbinrecread(stations_HRRR(i),5,(/nt_HRRR/),"float")
;  t_HRRR = 9./5.*t_HRRR + 32.
  rh_HRRR = fbinrecread(stations_HRRR(i),6,(/nt_HRRR/),"float")
  wdir_HRRR = fbinrecread(stations_HRRR(i),7,(/nt_HRRR/),"float")
  wspeed_HRRR = fbinrecread(stations_HRRR(i),8,(/nt_HRRR/),"float")

  name = ""+is
  
res@tmXBMode   = "Explicit"
res@tmXBValues = (/0., 6., 12., 18., 24.,30./)
res@tmXBLabels = (/"00","06","12","18","00","06"/)
res@tmXMajorGrid        = True
res@tmXMajorGridLineDashPattern = 2    ; dotted
res@tmXMajorGridThicknessF     = 1.0   ; optional
res@tmXMajorGridLineColor      = "gray50"  ; optional
res@tmYMajorGrid        = True
res@tmYMajorGridLineDashPattern = 2    ; dotted
res@tmYMajorGridThicknessF     = 1.0   ; optional
res@tmYMajorGridLineColor      = "gray50"  ; optional

 res@trYMinF = 10.
 res@trYMaxF = 40.
; res@trYMinF = 50.
; res@trYMaxF = 105.
;  if(iplot%5.eq.0) then
    res@tmYLLabelsOn               = True            ; turn off lat labels
;  else
;    res@tmYLLabelsOn               = False            ; turn off lat labels
;  end if
;  if(iplot.ge.ns-5) then
    res@tmXBLabelsOn               = True            ; turn off lon labels
;  else
;    res@tmXBLabelsOn               = False            ; turn off lon labels
;  end if
  plot(iplot)  = gsn_csm_xy (wks,hour_SAM,t_SAM,res)
  res@gsLineColor = "blue"
  dumplot(iplot)  = gsn_add_polyline (wks,plot(iplot),hour_HRRR,t_HRRR,res)
  poly(iplot)  = gsn_add_polymarker(wks,plot(iplot),hour_obs,t_obs,polyres)
  dum(iplot) = gsn_add_text(wks,plot(iplot),name,24.,37.,txres)
  res@trYMinF = 20.
  res@trYMaxF = 115.
  plot1(iplot)  = gsn_csm_xy (wks,hour_SAM,rh_SAM,res)
  dumplot1(iplot)  = gsn_add_polyline (wks,plot1(iplot),hour_HRRR,rh_HRRR,res)
  poly1(iplot)  = gsn_add_polymarker(wks,plot1(iplot),hour_obs,rh_obs,polyres)
  dum1(iplot) = gsn_add_text(wks,plot1(iplot),name,24.,110.,txres)
  res@trYMinF = 0.
  res@trYMaxF = 360.
  plot2(iplot)  = gsn_csm_xy (wks,hour_SAM,wdir_SAM,res)
  poly2(iplot)  = gsn_add_polymarker(wks,plot2(iplot),hour_obs,wdir_obs,polyres)
  dum2(iplot) = gsn_add_text(wks,plot2(iplot),name,24.,110.,txres)
  res@trYMinF = 0.
  res@trYMaxF = 10.
  plot3(iplot)  = gsn_csm_xy (wks,hour_SAM,wspeed_SAM,res)
  poly3(iplot)  = gsn_add_polymarker(wks,plot3(iplot),hour_obs,wspeed_obs,polyres)
  dum3(iplot) = gsn_add_text(wks,plot3(iplot),name,24.,110.,txres)
  iplot = iplot+1
  delete(hour_obs)
  delete(t_obs)
  delete(rh_obs)
  delete(wdir_obs)
  delete(wspeed_obs)
  delete(hour_SAM)
  delete(t_SAM)
  delete(rh_SAM)
  delete(wdir_SAM)
  delete(wspeed_SAM)
end do
gsn_panel(wks,plot,(/8,5/),resl)
gsn_panel(wks,plot1,(/8,5/),resl)
gsn_panel(wks,plot2,(/8,5/),resl)
gsn_panel(wks,plot3,(/8,5/),resl)



