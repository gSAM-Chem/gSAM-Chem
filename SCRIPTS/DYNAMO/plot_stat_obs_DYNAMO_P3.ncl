; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************
 ; set number of files:
 nfiles = 5

 filename = new((/nfiles/),string)

 ; set file names:
;filename(0) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_M2005_CAM.nc"
;filename(1) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_M2005_CAM_SCAU.nc"
;filename(2) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_M2005_CAM_SCAU_ICESNOW250.nc"
;filename(3) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_M2005_CAM_SCAU_ICESNOW250_JOONG.nc"
;filename(0) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_P3_RRTM.nc"
;filename(1) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_P3_RRTM_SCAU.nc"
;filename(2) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_P3_RRTM_SCAU_ICEF8.nc"
;filename(3) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_P3_RRTM_SCAU_ALLICERAD.nc"
;filename(4) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_P3_CAM_SCAU_ALLICERAD.nc"
;filename(5) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_P3_CAM_SCAU.nc"
;filename(0) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK.nc"
;filename(1) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU.nc"
;filename(2) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU_Nc100.nc"
;filename(3) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU_Nc100_icef8.nc"
;filename(4) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU_Nc100_iceifs.nc"
;filename(5) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU_Nc100_iceifs_qi2e-5.nc"
;filename(4) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU_Nc100_iceifs_qi2e-5.nc"
;filename(2) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU_Nc100_iceifs_qi2e-5.nc"
;filename(3) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_KK_SCAU_iceifs_qi2e-5.nc"
;filename(0) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR.nc"
;filename(1) = "./../../OUT_STAT/DYNAMO_32x32x72_4km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR.nc"
;filename(0) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM.nc"
;filename(2) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM.nc"
;filename(3) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR.nc"
;filename(4) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_TKF0.nc"
;filename(5) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1_RRTM_TVELTHR_TKF0_TG49.nc"
;filename(0) = "./../../OUT_STAT/DYNAMO_128x128x64_1km_10s_P3.nc"
;filename(1) = "./../../OUT_STAT/DYNAMO_128x128x64_1km_10s_P3_38C.nc"
;filename(0) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_P3.nc"
;filename(1) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_M2005_CAM.nc"
;filename(2) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM.SAM_CAM_SGSSAM.nc"
;filename(3) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_P3.SAM.nc"
filename(0) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_P3_CAM.nc"
filename(1) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_P3_RRTM.nc"
filename(2) = "./../../OUT_STAT/DYNAMO_128x128x64_1km_10s_P3_RRTM.nc"
filename(3) = "../../../SAM6.11.9/OUT_STAT/DYNAMO_128x128x64_1km_10s_P3.nc"
filename(4) = "../../../SAM6.11.9/OUT_STAT/DYNAMO_128x128x64_1km_10s_P3_RRTM.nc"
;filename(5) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_M2005_CAM_SGSSAM.nc"
;filename(5) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_SAM1MOM_RHHOMO1.nc"
;filename(5) = "./../../OUT_STAT/DYNAMO_128x128x72_1km_10s_M2005_CAM_SCAU.nc"
;filename(1) = "./../../OUT_STAT/DYNAMO_128x128x64_1km_10s_eis.nc"
;filename(3) = "./../../OUT_STAT/DYNAMO_512x512x64_0.25km_10s_P3.nc"
;filename(4) = "./../../OUT_STAT/DYNAMO_512x512x64_0.25m_10s_SAM1MOM.nc"
;filename(5) = "./../../OUT_STAT/DYNAMO_128x128x64_1km_10s_lsf3D.nc"


 Case = new((/nfiles/),string)
 ; Case names:
; Case(0) = "4km M2005 CAM"
; Case(1) = "4km M2005 CAM SCAU"
; Case(2) = "4km M2005 CAM SCAU ICESNOW250"
; Case(3) = "4km M2005 CAM SCAU ICESNOW250 JOONG"
; Case(0) = "4km P3"
; Case(1) = "4km P3 SCAU"
; Case(2) = "4km P3 SCAU ICEF8"
; Case(3) = "4km P3 SCAU ALLICERAD"
; Case(4) = "4km P3 CAM SCAU ALLICERAD"
; Case(5) = "4km P3 CAM SCAU"
; Case(0) = "4km Nc50"
; Case(1) = "4km Nc100"
; Case(2) = "1km Nc100"
; Case(3) = "1km Nc50"
; Case(0) = "SAM1MOM 4km KK"
; Case(1) = "SAM1MOM 4km KK SCAU"
; Case(2) = "SAM1MOM 4km KK SKAU NC100"
; Case(3) = "SAM1MOM 4km KK SKAU NC100 ICEF8"
; Case(4) = "SAM1MOM 4km KK SKAU NC100 ICEIFS"
; Case(5) = "SAM1MOM 4km KK SKAU NC100 ICEIFS QI2E-5"
; Case(4) = "SAM1MOM 1km KK SKAU NC100 ICEIFS QI2E-5"
; Case(0) = "SAM1MOM 1km"
; Case(2) = "SAM1MOM 1km RHHOMO1 RRTM"
; Case(3) = "SAM1MOM 1km RHHOMO1 RRTM TVELTHR"
; Case(4) = "SAM1MOM 1km RHHOMO1 RRTM TVELTHR_TKF0"
; Case(5) = "SAM1MOM 1km RHHOMO1 RRTM TVELTHR_TKF0 TG40"
; Case(0) = "P3 1km 72"
; Case(1) = "M2005 1km 72"
; Case(2) = "SAM1MOM.SAM_CAM_SGSSAM"
; Case(3) = "P3 1km 72 SAM"
 Case(0) = "P3 gSAM 1km 72 CAM"
 Case(1) = "P3 gSAM 1km 72 RRTM"
 Case(2) = "P3 gSAM 1km 64 RRTM"
 Case(3) = "P3 SAM 1km 64 CAM"
 Case(4) = "P3 SAM 1km 64 RRTM"
; Case(5) = "M2005 1km 72 SGSSAM"
; Case(5) = "SAM1MOM 1km RHHOMO1"
; Case(5) = "M2005 1km 72 SCAU"
; Case(6) = "SAM1MOM 1km RHHOMO1 RRTM TVELTHR_TKF0"
; Case(1) = "SAM1MOM EIS 1km"
; Case(3) = "P3 0.25km"
; Case(4) = "SAM1MOM 0.25km"
; Case(5) = "SAM1MOM LSF3D 1km"

 colors = (/"black","red","blue","blue","green","brown","cyan","magenta"/)
 dashes = (/0,0,0,1,0,0,0,0/)

; running mean interval (days)
 time_running = 3.0

 ; Set time interval (in days) for time averaged vertical profile plots:

; time_start = 278.
; time_stop = 285. 
; time_start = 295.
; time_stop = 305. 
 time_start = 280.
 time_stop = 360. 

 lon0 = 76.5
 lat0 = 3.

 fileobs = "./../../CASES/DYNAMO/dynamo.nc"
 day0 = 274.
 plotobscloud = True


 ; define threshold of total water (in g/m3) to defone PBL top:
 qt_threshold = 6.5

 ; Define maximum height (in m) plotted:

; z_max = 18000.
 z_max = 1000.

 ; Time evolution plots: Set min and max values:
 field_names = (/"SST","CLDSHD","CLDHI", "CLDMID","CLDLOW","RH", "LHF","SHF","PW","CWP","IWP","PREC","LWNT","SWNT","TOP"/)
 lows = (/         295., 0.,      0.,      0.,       0.,     50.,    0.,   0.,  20.,  0.,  0.,  0.,    0.,  100. ,-20. /)
 highs = (/        305., 1.2,     1.2,     1.2,      1.2,     100.,  250., 60., 90.,400., 400., 60.,  400., 500., 160./)

 ; Vertical profile plots. Set min and max values:
 field_names_z = (/"TABS","QV",  "QN", "CLD","TKH",  "RELH", "W2", "TLFLUX", "TVFLUX","QTFLUX","PRECIP","RADQR","QC","QI","QR","QN","QP","QP"/)
 lows_z = (/         190.,   0.,   0.,  0.,   0.,   0.,     0.,  -400.,      -20.,     -50.,    0. , -5., 0.,0.,0.,0.,0.,0./)
 highs_z = (/        310.,   20.,  0.2,  1.,  10.,  100.,    0.5,  100.,        50.,     400.,   30., 5., 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 /)
        
;===========================================================================

 fin = new((/nfiles/),file) 
 do i=0,nfiles-1
   fin(i)     = addfile (filename(i) , "r")             ; add file
 end do

fobs = addfile(fileobs,"r")
time_obs = day0+fobs->time_offset/86400.
dtime = 1./24.


wks = gsn_open_wks("ps","ps")
    
;************************************************ 
; plotting parameters
;************************************************

 res                   = True                       ; plot mods desired
 resl                   = True                       ; plot mods desired

; res@tiMainString      = fld@long_name             ; add title
 res@trYReverse        = False                       ; reverse Y-axis
; res@xyDashPatterns    = 1                         ; choose dash patterns
 res@xyLineThicknessF = 2.0
 res@gsnStringFontHeightF = 0.015
 res@tmXBLabelFontHeightF = 0.015                 ; resize tick labels
 res@tmYLLabelFontHeightF = 0.015
 res@tiXAxisFontHeightF   = 0.015
 res@tiMainFontHeightF   = 0.02
 res@tiYAxisFontHeightF   = 0.015
 res@vpWidthF            = 0.8           ; change aspect ratio of plot
 res@vpHeightF           = 0.5
 res@tmXMajorGrid      = True
 res@tmYMajorGrid      = True
 res@tmXMajorGridThicknessF = 0.5
 res@tmYMajorGridThicknessF = 0.5
 res@tmXMinorGridLineDashPattern = 1
 res@tmYMinorGridLineDashPattern = 1


 resl@gsLineThicknessF  = 2.                    ; line thicker
 res@gsnDraw     = False
 res@gsnFrame    = False


 res@xyExplicitLegendLabels = Case

 nploty = 3
 nplotx = 2
 nplotmax = nplotx * nploty ; number of plots per page
 plot =new((/nfiles*nplotmax*2/),graphic)

 do k=0,dimsizes(field_names)-1
   print(""+field_names(k))
   do i=0,nfiles-1
     obs_plot = False
     fi = fin(i)
     res@xyLineColor    =  colors(i+1)
     res@xyDashPattern    =   dashes(i+1)
     resl@gsLineDashPattern = dashes(i+1)                    ; solid line
     resl@gsLineColor       = colors(i+1)                 ; line color
     res@trYMinF = lows(k)
     res@trYMaxF = highs(k)
     res@tiXAxisString = "Days"
   ;  res@pmLegendOrthogonalPosF = -0.8-0.2*i
     if(field_names(k).eq."TKE") then
        z = fi->z
        adz = 0.5*(z(0)+z(1))
        fld = fi->time
        fld = 0.
        do m=1,dimsizes(z)-2
          fld = fld + fi->TKE(:,m)*adz*fi->RHO(:,m)
          adz = 0.5*(z(m+1)-z(m-1))
        end do
        res@tiYAxisString = "kg/m/s2"
        res@tiMainString = "Vertically Integrated TKE"
        delete(z)
     else
     if(field_names(k).eq."Zi") then
        fld = fi->time
        fld = 0.
        z = fi->z
        do n=0,dimsizes(fld)-1
          fld(n) = 0.
          do m=dimsizes(z)-2,1,1
            if(fi->QT(n,m).gt.qt_threshold) then
              fld(n) = z(m-1)-(fi->QT(n,m-1)-fi->QT(n,m))/(z(m-1)-z(m))*(qt_threshold-fi->QT(n,m))
              break
            end if
          end do
        end do
        res@tiYAxisString = "m"
        res@tiMainString = "PBL Top"
        delete(z)
      else
     if(field_names(k).eq."W2MAX") then
        fld = fi->time
        fld = 0.
        fld =  dim_max(fi->W2)
        res@tiYAxisString = "m2/s2"
        res@tiMainString = "Maximum Vertical Velocity Variance"
     else
     if(field_names(k).eq."RH") then
        fld = fi->time
        fld = 0.
        fld =  fi->RELH(:,0)
        res@tiYAxisString = "%"
        res@tiMainString = "Surface Air RH"
     else
     if(field_names(k).eq."TOP") then
        fld = fi->time
        fld = 0.
        fld =  fi->SWNT(:) - fi->LWNT(:)
        res@tiYAxisString = "W/m2"
        res@tiMainString = "NET TOA (ASR - OLR)"
     else
        res@tiYAxisString = fi->$field_names(k)$@units
        res@tiMainString = fi->$field_names(k)$@long_name
        fld = fi->$field_names(k)$(:)
      end if
      end if
      end if
      end if
      end if
   nrunning = round(time_running/(fi->time(1)-fi->time(0)),3)
;   print(time_running+" "+(fi->time(1)-fi->time(0)))
   print("nrunning="+nrunning)
   wgt=new((/nrunning/),float)
   wgt = 1./nrunning
   if(nrunning.gt.1) then
     fld = wgt_runave(fld,wgt,0)
     fld = wgt_runave(fld,wgt,0)
   end if
   delete(wgt)

      if(i.eq.nfiles-1.and.field_names(k).eq."SST") then
        fldobs = fobs->T_skin(:)+273.15
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."LHF") then
        fldobs = fobs->LH(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."SHF") then
        fldobs = fobs->SH(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."CWP") then
        fldobs = fobs->LWP(:)*10000.
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."RH") then
      ;  fldobs = fobs->RH_srf(:)
      ;  timeobs = time_obs
        fldobs = relhum(fi->TABSOBS(:,0),fi->QVOBS(:,0)*1.e-3,fi->Ps(:)*100.)
        timeobs = fi->time
        obs_plot = True
      end if
    if(plotobscloud) then
      if(i.eq.nfiles-1.and.field_names(k).eq."CLDLOW") then
        fldobs = fobs->cld_low(:)*0.01
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."CLDMID") then
        fldobs = fobs->cld_mid(:)*0.01
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."CLDHI") then
        fldobs = fobs->cld_high(:)*0.01
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."CLDSHD") then
        fldobs = fobs->cld_tot(:)*0.01
        timeobs = time_obs
        obs_plot = True
      end if
    end if
      if(i.eq.nfiles-1.and.field_names(k).eq."PREC") then
        fldobs = fobs->prec_srf(:)*24.
        timeobs = time_obs
;        print(timeobs+"   "+fldobs)
        obs_plot = True
        print("total OBS PREC="+sum(fldobs)/(max(timeobs)-min(timeobs)))
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."PW") then
        fldobs = fi->PWOBS(:)
        timeobs = fi->time
;         fldobs = short2flt(fobs1->tcwv(:,{lat0},{lon0}))
;         timeobs = time_obs1
;         print(timeobs+" "+fldobs)
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."LWNT") then
        fldobs = fobs->lw_net_toa(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."SWNT") then
        fldobs = fobs->sw_net_toa(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."TOP") then
        fldobs = fobs->sw_net_toa(:) - fobs->lw_net_toa(:)
        timeobs = time_obs
        obs_plot = True
      end if

;      if(field_names(k).eq."PREC") then
;         print("total PREC="+sum(fld)/(max(fi->time)-min(fi->time)))
;      end if
      print(field_names(k)+"  "+avg(fld({time_start:time_stop})))

      if(i.eq.0) then
         plot(0)  = gsn_csm_xy (wks,fi->time(:),fld,res)
      else
         plot(i)  = gsn_add_polyline (wks,plot(0),fi->time(:),fld,resl)
      end if
      if(obs_plot) then
        ;   nrunning = round(time_running/(fobs->time_offset(1)-fobs->time_offset(0))*86400.,3)
           nrunning = round(time_running/(timeobs(1)-timeobs(0)),3)
;           print("nrunning_ob="+nrunning)
           wgt=new((/nrunning/),float)
           wgt = 1./nrunning
           fldobs = wgt_runave(fldobs,wgt,0)
           delete(wgt)
           res@xyLineColor    =  colors(0)
           res@xyDashPattern    =   dashes(0)
           resl@gsLineDashPattern = dashes(0)                    ; solid line
           resl@gsLineColor       = colors(0)                 ; line color
           resl@gsLineThicknessF  = 3.                    ; line thicker
         plot(i+1)  = gsn_add_polyline (wks,plot(0),timeobs(:),fldobs,resl)
           resl@gsLineThicknessF  = 2.                    ; line thicker
         delete(fldobs)
         delete(timeobs)
      end if
      delete(fld)
   end do
   draw (plot(0))
   yy=0.75
   if(nfiles.gt.1) then
     do i=0,nfiles-1
      res@txFontHeightF   = 0.01
      res@txJust   = "BottomLeft"
      resl@gsLineDashPattern = dashes(i+1)                    ; solid line
      resl@gsLineColor       = colors(i+1)                 ; line color
      gsn_polyline_ndc(wks,(/0.21,0.25/),(/yy-0.02*i,yy-0.02*i/),resl)
      gsn_text_ndc(wks,Case(i),0.26,yy-0.02*i,res)
     end do
   end if
   frame(wks)

 end do

 plots =new(nplotmax,graphic)
 res@vpWidthF            = 0.8           ; change aspect ratio of plot
 res@vpHeightF           = 0.5

 iplot=0
 nplots=0
 do k=0,dimsizes(field_names_z)-1
   print(""+field_names_z(k))
   do i=0,nfiles-1
     fi = fin(i)
     res@xyLineColor    =  colors(i+1)
     res@xyDashPattern    =   dashes(i+1)
     resl@gsLineDashPattern = dashes(i+1)                    ; solid line
     resl@gsLineColor       = colors(i+1)                 ; line color
     res@trYMinF = 0.
     res@trYMaxF = max(fi->z({:z_max}))
     res@trXMinF = lows_z(k)
     res@trXMaxF = highs_z(k)
     res@tiYAxisString = "height, m"
     if(field_names_z(k).eq."RADFLUX") then
        fld = fi->RADLWDN({time_start:time_stop},{:z_max})
        fld = fld - fi->RADLWUP({time_start:time_stop},{:z_max})
        fld = fld + fi->RADSWDN({time_start:time_stop},{:z_max})
        fld = fld - fi->RADSWUP({time_start:time_stop},{:z_max})
        res@tiXAxisString = "W/m2"
        res@tiMainString = "Net Radiative Flux"
     else if(field_names_z(k).eq."RELH") then
        fld = relhum(fi->TABS({time_start:time_stop},{:z_max}), fi->QV({time_start:time_stop},{:z_max})*1.e-3, conform(fi->QV({time_start:time_stop},{:z_max}),fi->p({:z_max}),1)*100.)
      ;  fld = fi->RELH({time_start:time_stop},{:z_max})
        res@tiXAxisString = "RH"
        res@tiMainString = "Relative Hum over Liq.W"
        fld!0 = "time"
        fld!1 = "z"
     else if(field_names_z(k).eq."Ri") then
        fld=fi->TABS({time_start:time_stop},{:z_max})
        fld = 0.
     else
        fld = fi->$field_names_z(k)$({time_start:time_stop},{:z_max})
        res@tiXAxisString = fi->$field_names_z(k)$@units
        res@tiMainString = fi->$field_names_z(k)$@long_name
     end if
     end if
     end if
     print(field_names_z(k)+" "+min(fld)+" "+max(fld))
     if(i.eq.0) then
       plot(iplot)  = gsn_csm_xy (wks,dim_avg(fld(z|:,time|:)),fi->z({:z_max}), res)
       if(field_names_z(k).eq."TABS") then
          fldobs = fi->TABSOBS({time_start:time_stop},{:z_max})
           res@xyLineColor    =  colors(0)
           res@xyDashPattern    =   dashes(0)
           resl@gsLineDashPattern = dashes(0)                    ; solid line
           resl@gsLineColor       = colors(0)                 ; line color
           nplots = nplots+1
           plot(nplotmax-1+nplots) = gsn_add_polyline (wks,plot(iplot),dim_avg(fldobs(z|:,time|:)),fi->z({:z_max}), resl)
;           print(dim_avg(fldobs(z|:,time|:)))
          delete (fldobs)
       end if
       if(field_names_z(k).eq."THETA") then
          fldobs = fi->TABSOBS({time_start:time_stop},{:z_max})
           res@xyLineColor    =  colors(0)
           res@xyDashPattern    =   dashes(0)
           resl@gsLineDashPattern = dashes(0)                    ; solid line
           resl@gsLineColor       = colors(0)                 ; line color
           nplots = nplots+1
           coefz = (1000./(fi->p({:z_max})))^(287./1004.)
           plot(nplotmax-1+nplots) = gsn_add_polyline (wks,plot(iplot),dim_avg(fldobs(z|:,time|:))*coefz,fi->z({:z_max}), resl)
;           print(dim_avg(fldobs(z|:,time|:)))
          delete (fldobs)
       end if
       if(field_names_z(k).eq."Ri") then
           coefz = (fi->p(1)/(fi->p({:z_max})))^(287./1004.)
           fldobs = fi->TABSOBS({time_start:time_stop},{:z_max})
           tobs = dim_avg(fldobs(z|:,time|:))*coefz
           fldobs = fi->QVOBS({time_start:time_stop},{:z_max})
           qobs = dim_avg(fldobs(z|:,time|:))*1.e-3
           tobs = tobs*(1.+0.61*qobs)
           fldobs = fi->UOBS({time_start:time_stop},{:z_max})
           uobs = dim_avg(fldobs(z|:,time|:))
           fldobs = fi->VOBS({time_start:time_stop},{:z_max})
           vobs = dim_avg(fldobs(z|:,time|:))
           velobs = sqrt(uobs^2+vobs^2)
           ts = (avg(fobs->T_skin(:))+273.15)*(1.+0.61*qobs(1)*1.2)
           z = fi->z({:z_max})
           Riobs = coefz
           Riobs(0) = 9.81/tobs(0)*(tobs(0)-ts)*z(0)/velobs(0)^2
           do m=1,dimsizes(coefz)-1
            Riobs(m) = 9.81/tobs(m)*(tobs(m)-tobs(m-1))*(z(m)-z(m-1))/(velobs(m)-velobs(m-1))^2
           end do
           do m=0,dimsizes(coefz)-1
            print(z(m)+"    "+Riobs(m)+"  "+tobs(m)+"   "+(qobs(m)*1.e3)+"  "+velobs(m))
           end do
           res@xyLineColor    =  colors(0)
           res@xyDashPattern    =   dashes(0)
           resl@gsLineDashPattern = dashes(0)                    ; solid line
           resl@gsLineColor       = colors(0)                 ; line color
           nplots = nplots+1
           plot(nplotmax-1+nplots) = gsn_add_polyline (wks,plot(iplot),Riobs,fi->z({:z_max}), resl)
;           print(dim_avg(fldobs(z|:,time|:)))
          delete (fldobs)
       end if
       if(field_names_z(k).eq."QV") then
          fldobs = fi->QVOBS({time_start:time_stop},{:z_max})
           res@xyLineColor    =  colors(0)
           res@xyDashPattern    =   dashes(0)
           resl@gsLineDashPattern = dashes(0)                    ; solid line
           resl@gsLineColor       = colors(0)                 ; line color
           nplots = nplots+1
           plot(nplotmax-1+nplots) = gsn_add_polyline (wks,plot(iplot),dim_avg(fldobs(z|:,time|:)),fi->z({:z_max}), resl)
;           print(dim_avg(fldobs(z|:,time|:)))
          delete (fldobs)
       end if
       if(field_names_z(k).eq."RELH") then
          fldobs = relhum(fi->TABSOBS({time_start:time_stop},{:z_max}), fi->QVOBS({time_start:time_stop},{:z_max})*1.e-3, conform(fi->QVOBS({time_start:time_stop},{:z_max}),fi->p({:z_max}),1)*100.)
           fldobs!0 = "time"
           fldobs!1 = "z"
           res@xyLineColor    =  colors(0)
           res@xyDashPattern    =   dashes(0)
           resl@gsLineDashPattern = dashes(0)                    ; solid line
           resl@gsLineColor       = colors(0)                 ; line color
           nplots = nplots+1
           plot(nplotmax-1+nplots) = gsn_add_polyline (wks,plot(iplot),dim_avg(fldobs(z|:,time|:)),fi->z({:z_max}), resl)
;           print(dim_avg(fldobs(z|:,time|:)))
          delete (fldobs)
       end if
     else
       nplots = nplots+1
       plot(nplotmax-1+nplots)  = gsn_add_polyline (wks,plot(iplot),dim_avg(fld(z|:,time|:)),fi->z({:z_max}), resl)
     end if
     delete (fld)
   end do
   plots(iplot)=plot(iplot)
   iplot=iplot+1
   if(iplot.eq.nplotmax) then
     iplot=0
     nplots=0
     resP = True
     if(nfiles.gt.1) then
      yy = 0.95
      do i=0,nfiles-1
       res@txFontHeightF   = 0.01
       res@txJust   = "BottomLeft"
       resl@gsLineDashPattern = dashes(i+1)                    ; solid line
       resl@gsLineColor       = colors(i+1)                 ; line color
       gsn_polyline_ndc(wks,(/0.2,0.25/),(/yy-0.02*i,yy-0.02*i/),resl)
       gsn_text_ndc(wks,Case(i),0.26,yy-0.02*i,res)
      end do
     end if
     gsn_panel(wks,plots,(/nploty,nplotx/),False)
   end if
 end do




end

