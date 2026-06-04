; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************
 ; set number of files:
 nfiles = 6

 filename = new((/nfiles/),string)

 ; set file names:

 filename(0) = "../../OUT_STAT/PBL_DRY_288x288x150_20m_20m.nc"
 filename(1) = "../../OUT_STAT/PBL_DRY_64x64x57_100m.nc"
 filename(2) = "../../OUT_STAT/PBL_DRY_64x64x57_100m_nudgeuv_3600s.nc"
 filename(3) = "../../OUT_STAT/PBL_DRY_64x64x57_100m_nudgeuv_1800s.nc"
 filename(4) = "../../OUT_STAT/PBL_DRY_64x64x57_100m_nudgeuv_900s.nc"
 filename(5) = "../../OUT_STAT/PBL_DRY_64x64x57_100m_nudgeuv_450s.nc"
; filename(1) = "../../OUT_STAT/PBL_DRY_64x64x57_100m_Cs0.19.nc"
; filename(2) = "../../OUT_STAT/PBL_DRY_64x64x57_100m_Cs0.19_Pr2.nc"
; filename(3) = "../../OUT_STAT/PBL_DRY_64x64x57_100m_Cs0.19_Pr3.nc"

; filename(3) = "./../../../SAM6.11.7/OUT_STAT/PBL_DRY_40x40x40_160m_60m_1.5TKE.nc"
 Case = new((/nfiles/),string)

 ; Case names:
 Case(0) = "20m"
 Case(1) = "100m LI"
 Case(2) = "100m LI_3600"
 Case(3) = "100m LI_1800"
 Case(4) = "100m LI_900"
 Case(5) = "100m LI_900"
; Case(1) = "100m LI_Cs0.19"
; Case(2) = "100m LI_Cs0.19 Pr2"
; Case(3) = "100m LI_Cs0.19 Pr3"

 ; Set time interval (in days) for time averaged vertical profile plots:

time_start=1.+10800./86400.
time_stop=1.+14400./86400.
z_max = 3000.

 ; Vertical profile plots. Set min and max values:
 field_names_z = (/"THETA", "W2", "TK", "U2", "TLFLUX", "TKE",  "TLFLUXS","TKES"/)
 lows_z = (/         300.,   0.,   0.,    0.,  -20.,   0.,        -20.,    0. /)
 highs_z = (/        301.,   2.,   50.,    2.,  100.,   2.,       100.,    2. /)
        
 colors = (/"black","blue","green","red","blue","green"/)
 dashes = (/0,0,0,0,1,1/)
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
 res@xyLineThicknessF = 3.0
 res@gsnStringFontHeightF = 0.015
 res@tmXBLabelFontHeightF = 0.015                 ; resize tick labels
 res@tmYLLabelFontHeightF = 0.015
 res@tiXAxisFontHeightF   = 0.015
 res@tiMainFontHeightF   = 0.02
 res@tiYAxisFontHeightF   = 0.015
 res@vpWidthF            = 0.4           ; change aspect ratio of plot
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

 nplots = 4
 plot =new((/nfiles*nplots/),graphic)
 plots =new(nplots,graphic)

 iplot=0
 do k=0,dimsizes(field_names_z)-1
   print(""+field_names_z(k))
   do i=0,nfiles-1
     fi = fin(i)
     res@xyLineColor    =  colors(i)
     res@xyDashPattern    =   dashes(i)
     resl@gsLineDashPattern = dashes(i)                    ; solid line
     resl@gsLineColor       = colors(i)                 ; line color
     res@trYMinF = 0.
     res@trYMaxF = max(fi->z({:z_max}))
     res@trXMinF = lows_z(k)
     res@trXMaxF = highs_z(k)
     res@tiYAxisString = "height, m"
     if(field_names_z(k).eq."TKE") then
        fld = fi->TKE({time_start:time_stop},{:z_max})
        fld = fld + fi->TKES({time_start:time_stop},{:z_max})
        res@tiXAxisString = "m2/s2"
        res@tiMainString = "Total TKE"
     else
        fld = fi->$field_names_z(k)$({time_start:time_stop},{:z_max})
        res@tiXAxisString = fi->$field_names_z(k)$@units
        res@tiMainString = fi->$field_names_z(k)$@long_name
     end if
     if(i.eq.0) then
       plot(iplot)  = gsn_csm_xy (wks,dim_avg(fld(z|:,time|:)),fi->z({:z_max}), res)
     else
       print(iplot+"   "+i+"   "+(6+i+iplot*(nfiles-1)))
       plot(nplots-1+i+iplot*(nfiles-1))  = gsn_add_polyline (wks,plot(iplot),dim_avg(fld(z|:,time|:)),fi->z({:z_max}), resl)
     end if
     delete (fld)
   end do
   plots(iplot)=plot(iplot)
   iplot=iplot+1
   if(iplot.eq.nplots) then
     iplot=0
     resP = True
     if(nfiles.gt.1) then
      yy = 0.95
      do i=0,nfiles-1
       res@txFontHeightF   = 0.01
       res@txJust   = "BottomLeft"
       resl@gsLineDashPattern = dashes(i)                    ; solid line
       resl@gsLineColor       = colors(i)                 ; line color
       gsn_polyline_ndc(wks,(/0.7,0.75/),(/yy-0.02*i,yy-0.02*i/),resl)
       gsn_text_ndc(wks,Case(i),0.76,yy-0.02*i,res)
      end do
     end if
     gsn_panel(wks,plots,(/nplots/2,2/),False)
   end if
 end do




end

