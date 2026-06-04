; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************
 ; set number of files:
 nfiles = 1

 filename = new((/nfiles/),string)

 ; set file names:
 filename(0) = "./../../OUT_2DL/DIURNAL_72x72x70_3km_land12_36_AVG.2D_lnd.nc"

 finv = addfile("./../../OUT_INV/DIURNAL_72x72x70_3km_land12_36_INV.2D_lnd.nc","r")

 Case = new((/nfiles/),string)
 ; Case names:
 Case(0) = "LAND12"

; running mean interval (days)
 time_running = 1./24.


 fileobs = "/glade/work/marat/ERA5/ERA5_July2020_-90_37.nc"
 fobs = addfile(fileobs,"r")
 time_obs = fobs->time


 ; Set time interval (in days) for time averaged vertical profile plots:

 ; Time evolution plots: Set min and max values:
 field_names = (/"LHFC","LHFCW","LHFS","SHFC","SHFS","TCAS","TCNP", "TS1","WS1","TS6","WS6","QCAS","GRFLUX0","Rc","Rb","TOTSOIL","LAI","USTAR","MW","PCNP","DCNP"/)
 lows = (/        -100., -100.,  -100., -100., -100.,   290.,  290.,   290., 0.,  290.,  0.,   0.,    -200.,    0.,  0.,    0.,     0.,    0.,    0.,   0.,     0./)
 highs = (/        300.,  300.,  300.,  300.,  300.,   340.,  340.,   340.,  1., 340.,  1.,   25.,    200. , 100., 100.,   0.5,    10.,   0.5,   1.,  1.,   1./)

print(">>>"+dimsizes(field_names)+"  "+dimsizes(lows)+"  "+dimsizes(highs))
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
     res@tiYAxisString = fi->$field_names(k)$@units
     res@tiMainString = fi->$field_names(k)$@long_name
     if(field_names(k).eq."TOTSOIL") then
      fld = fi->WS1(:)*finv->S_DZ1(0,1,1)*finv->PORO1(0,1,1)
      do m=2,9
        fld = fld+fi->$("WS"+m)$(:) * finv->$("S_DZ"+m)$(0,1,1)*finv->$("PORO"+m)$(0,1,1)
      end do 
      print(min(fld)+" "+max(fld))
     else
       fld = fi->$field_names(k)$(:)
     end if

   nrunning = round(time_running/(fi->time(1)-fi->time(0)),3)
   print(time_running+" "+(fi->time(1)-fi->time(0)))
   print("nrunning="+nrunning)
   wgt=new((/nrunning/),float)
   wgt = 1./nrunning
   fld = wgt_runave(fld,wgt,0)
   delete(wgt)

      if(i.eq.nfiles-1.and.field_names(k).eq."TS1") then
        fldobs = fobs->stl1(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."WS1") then
        fldobs = fobs->swvl1(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."TS6") then
        fldobs = fobs->stl2(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."WS6") then
        fldobs = fobs->swvl2(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."LAI") then
        fldobs = fobs->lai_hv(:)
        fldobs = fldobs+ fobs->lai_lv(:)
        timeobs = time_obs
        obs_plot = True
      end if
      if(i.eq.nfiles-1.and.field_names(k).eq."USTAR") then
        fldobs = fobs->zust(:)
        timeobs = time_obs
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


end

