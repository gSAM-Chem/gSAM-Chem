; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************
 ; set number of files:
 nfiles = 3

 filename = new((/nfiles/),string)

 ; set file names:
 filename(0) = "./../../OUT_STAT/OCEAN_128x128x83_100m_3s.nc"
 filename(1) = "./../../OUT_STAT/OCEAN_128x128x83_100m_3s_LSF.nc"
 filename(2) = "./../../OUT_STAT/OCEAN_128x128x83_100m_3s_NUDGEUVTQ.nc"

 Case = new((/nfiles/),string)
 ; Case names:
 Case(0) = "CTRL"
 Case(1) = "LHF"
 Case(2) = "NUDGEUVTQ"

; running mean interval (days)
; time_running = 1.
 time_running = 1./24.

 ; Set time interval (in days) for time averaged vertical profile plots:

 time_start = 190.5
 time_stop = 191.


 ; define threshold of total water (in g/m3) to defone PBL top:
 qt_threshold = 6.5

 ; Define maximum height (in m) plotted:

 z_max = 3000.
; z_max = 19000.

 ; Time evolution plots: Set min and max values:
 field_names = (/"SST","T2m","VEL10m","CLDSHD","CLDHI", "CLDLOW","RH", "LHF","SHF","TAU","PW","CWP","WMAX","PREC","TKE","LWNT","SWNT","SFC"/)
 lows = (/        290., 290.,  -20.,    0.,     0.,       0.,     40.,   0.,   0.,   0.,  0.,  0.,    0.,   0.,    0.,       0.,    0.,  -100./)
 highs = (/       300., 300.,   20.,    1.,     1.,       1.,    100.,  300., 100.,  0.1, 70., 800.,  70., 100., 100000.,    400.,   400., 200./)

 ; Vertical profile plots. Set min and max values:
 field_names_z = (/"TABS","QV",  "QN", "CLD","TKH",  "RELH", "W2", "TLFLUX", "TVFLUX","QTFLUX","PRECIP","RADQR","QC","QI","QR","QN","QP","QP"/)
 lows_z = (/         190.,   0.,   0.,  0.,   0.,   0.,     0.,  -400.,      -20.,     -50.,    0. , -5., 0.,0.,0.,0.,0.,0./)
 highs_z = (/        310.,   20.,  0.2,  1.,  10.,  100.,    0.5,  100.,        50.,     400.,   30., 5., 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 /)
        
 colors = (/"black","red","blue","green","magenta","brown"/)
 dashes = (/0,0,0,0,0,0,0/)
;===========================================================================

 fin = new((/nfiles/),file) 
 do i=0,nfiles-1
   fin(i)     = addfile (filename(i) , "r")             ; add file
 end do


wks = gsn_open_wks("ps","ps")
    
;************************************************ 
; plotting parameters
;************************************************

 res                   = True                       ; plot mods desired
 resl                   = True                       ; plot mods desired

; res@tiMainString      = fld@long_name             ; add title
 res@trYReverse        = False                       ; reverse Y-axis
; res@xyDashPatterns    = 1                         ; choose dash patterns
 res@xyLineThicknessF = 4.0
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


 resl@gsLineThicknessF  = 3.                    ; line thicker
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
     res@pmLegendOrthogonalPosF = -0.4-0.2*i
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
     if(field_names(k).eq."TAU") then
        fld = fi->UW(:,0)
        fld = sqrt(fld^2 + fi->VW(:,0)^2)
        fld = fld * fi->RHO(:,0)
        res@tiYAxisString = "N/m2"
        res@tiMainString = "Surface Stress"
      else
     if(field_names(k).eq."VEL10m") then
        fld = fi->time
        fld = 0.
        fld =  fi->U10m
        fld =  sqrt(fld^2 + fi->V10m^2)
        res@tiYAxisString = "m/s"
        res@tiMainString = "10-m Wind"
     else
     if(field_names(k).eq."RH") then
        fld = fi->time
        fld = 0.
        fld =  fi->RELH(:,0)
        res@tiYAxisString = "%"
        res@tiMainString = "Surface Air RH"
     else
     if(field_names(k).eq."SFC") then
        fld = fi->time
        fld = 0.
        fld =  fi->SWNS
        fld =  fld - fi->LWNS
        fld =  fld - fi->LHF
        fld =  fld - fi->SHF
        res@tiYAxisString = "W/m2"
        res@tiMainString = "SFC SW - LW - LHF - SHF"
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
   print(time_running+" "+(fi->time(1)-fi->time(0)))
   print("nrunning="+nrunning)
   wgt=new((/nrunning/),float)
   wgt = 1./nrunning
   fld = wgt_runave(fld,wgt,0)
   delete(wgt)

      if(i.eq.nfiles-1.and.field_names(k).eq."T2m") then
        fldobs = fi->T2m_OBS(:)
        timeobs = fi->time
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."VEL10m") then
        fldobs = fi->VEL10m_OBS(:)
        timeobs = fi->time
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."SHF") then
        fldobs = fi->SHFOBS(:)
        timeobs = fi->time
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."LHF") then
        fldobs = fi->LHFOBS(:)
        timeobs = fi->time
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."TAU") then
        fldobs = fi->TAU_OBS(:)
        timeobs = fi->time
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."PW") then
        fldobs = fi->PWOBS(:)
        timeobs = fi->time
        obs_plot = True
      end if


      if(i.eq.0) then
         plot(0)  = gsn_csm_xy (wks,fi->time(:),fld,res)
      else
         plot(i)  = gsn_add_polyline (wks,plot(0),fi->time(:),fld,resl)
      end if
      if(obs_plot) then
        ;   nrunning = round(time_running/(fobs->time_offset(1)-fobs->time_offset(0))*86400.,3)
           nrunning = round(time_running/(timeobs(1)-timeobs(0)),3)
           print("nrunning_ob="+nrunning)
           wgt=new((/nrunning/),float)
           wgt = 1./nrunning
           fldobs = wgt_runave(fldobs,wgt,0)
           delete(wgt)
           res@xyLineColor    =  colors(0)
           res@xyDashPattern    =   dashes(0)
           resl@gsLineDashPattern = dashes(0)                    ; solid line
           resl@gsLineColor       = colors(0)                 ; line color
         plot(i+1)  = gsn_add_polyline (wks,plot(0),timeobs(:),fldobs,resl)
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
      gsn_polyline_ndc(wks,(/0.8,0.85/),(/yy-0.02*i,yy-0.02*i/),resl)
      gsn_text_ndc(wks,Case(i),0.86,yy-0.02*i,res)
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
       gsn_polyline_ndc(wks,(/0.8,0.85/),(/yy-0.02*i,yy-0.02*i/),resl)
       gsn_text_ndc(wks,Case(i),0.86,yy-0.02*i,res)
      end do
     end if
     gsn_panel(wks,plots,(/nploty,nplotx/),False)
   end if
 end do




end

