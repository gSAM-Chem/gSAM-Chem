; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
; read in data
;************************************************
 filename = "./../../OUT_STAT/RCE_SAM_72x72x64_27km_cup_epmax0.95_rad_minorig5.nc"
; filename = "./../../OUT_STAT/RCE_SAM_144x144x64_1km_SAM1MOM_rad.nc"

 ; Set time interval (in hours) for time averaged vertical profile plots:

time_start = 0.
time_end =  5.

; Define maximum height (in m) plotted:

zmax = 1000.

names=(/" CUP"," Trans"," ToRain"," RainEv"," Residual"/)
colors = (/"blue","red","green","magenta","cyan","brown","blueviolet"/)
dashes = (/0,1,2,3,4,0,0,0/)

wks = gsn_open_wks("ps","ps")
    
;************************************************ 
; plotting parameters
;************************************************

 res                   = True                       ; plot mods desired

 res@trYReverse        = False                       ; reverse Y-axis
 res@xyLineThicknessF = 4.0
 res@gsnStringFontHeightF = 0.015
 res@tmXBLabelFontHeightF = 0.015                 ; resize tick labels
 res@tmYLLabelFontHeightF = 0.015
 res@tiXAxisFontHeightF   = 0.015
 res@tiMainFontHeightF   = 0.02
 res@tiYAxisFontHeightF   = 0.015
 res@vpWidthF            = 0.6           ; change aspect ratio of plot
 res@vpHeightF           = 0.75
 res@tmXMajorGrid      = True
 res@tmYMajorGrid      = True
 res@tmXMajorGridThicknessF = 0.5
 res@tmYMajorGridThicknessF = 0.5
 res@tmXMinorGridLineDashPattern = 1
 res@tmYMinorGridLineDashPattern = 1


 f = addfile(filename,"r")
 nz = dimsizes(f->z)

 fld = new((/5,nz/),float)
 time0 = f->time(0)

 fld(0,:) = dim_avg_n_Wrap(f->QNUDGE({time_start+time0:time_end+time0},:),0)
; fld(0,:) = fld(0,:) /(12.*12.)
 print(fld(0,:))
 fld(1,:) = dim_avg_n_Wrap(f->QTADV({time_start+time0:time_end+time0},:),0)
 fld(1,:) = fld(1,:) + dim_avg_n_Wrap(f->QTDIFF({time_start+time0:time_end+time0},:),0)
 fld(2,:) = dim_avg_n_Wrap(f->QTSINK({time_start+time0:time_end+time0},:),0)
 fld(3,:) = dim_avg_n_Wrap(f->QTSRC({time_start+time0:time_end+time0},:),0)
 fld(4,:) = fld(0,:)+fld(1,:)+fld(2,:)+fld(3,:)

 res@xyLineColors    =  colors
 res@xyDashPatterns    =   dashes
 res@trXMinF = -10.
 res@trXMaxF = 10.
 res@pmLegendDisplayMode    = "Always"
 res@pmLegendSide           = "Top"                 ; Change location of
 res@pmLegendParallelPosF   = .70                   ; move units right
 res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down

 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.20                  ; height of legend.
 res@lgLabelFontHeightF     = .02                   ; change font height
 res@lgPerimOn              = False                 ; no box around

 res@xyExplicitLegendLabels = names
 plot  = gsn_csm_xy (wks,fld(:,{0.:zmax}),f->z({0.:zmax}),res)
