; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;************************************************
begin
;************************************************
; read in data
;************************************************

zmax = 2000.
wstar = 1.46
zstar = 1350.
fluxt0 = 60.

 f = new((/10/),file)

; all data should have identical vertical grid

 nfiles = 1
 f(0) = addfile("../../OUT_STAT/PBL_DRY_128x128x64_50m.nc","r")

time1=9900./86400.
time2=10800./86400.
zmax = 3000.
ff=f(0)
time = ff->time
print(""+time)

wks = gsn_open_wks("ps","ps")

;************************************************
; plotting parameters
;************************************************

 res                   = True                       ; plot mods desired

; res@tiMainString      = fld@long_name             ; add title
 res@trYReverse        = False                       ; reverse Y-axis
; res@xyDashPatterns    = 1                         ; choose dash patterns
 res@pmLegendDisplayMode    = "Always"
 res@pmLegendSide           = "Top"                 ; Change location of 
 res@pmLegendParallelPosF   = .70                   ; move units right
 res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down
 
 res@pmLegendWidthF         = 0.12                  ; Change width and
 res@pmLegendHeightF        = 0.20                  ; height of legend.
 res@lgLabelFontHeightF     = .02                   ; change font height
 res@lgPerimOn              = False                 ; no box around
 res@xyLineThicknessF = 4.0 

 res@xyExplicitLegendLabels = (/"144x144x64 50m"/)
 ff = f(0)
 z = ff->z({0.:zmax}) ;/zstar
 data = new((/nfiles,dimsizes(z)/),float)

 do i=0,nfiles-1
   ff=f(i)
   fld = ff->TL({time1:time2},{0.:zmax})
   data(i,:) = dim_avg(fld(z|:,time|:))
 end do
 data@long_name = fld@long_name
 res@xyLineColors    =   (/1,2,3,4,6,7,8/) 
 plot  = gsn_csm_xy (wks,data,z,res) ; create plot

 do i=0,nfiles-1
   ff=f(i)
   fld = ff->TLFLUX({time1:time2},{0.:zmax})
   fld = fld/fluxt0
   data(i,:) = dim_avg(fld(z|:,time|:))
 end do
 data@long_name = fld@long_name
 res@xyLineColors    =   (/1,2,3,4,6,7,8/)
 plot  = gsn_csm_xy (wks,data,z,res) ; create plot

 do i=0,nfiles-1
   ff=f(i)
   fld = ff->W2({time1:time2},{0.:zmax})
   fld = fld/wstar^2
   data(i,:) = dim_avg(fld(z|:,time|:))
 end do
 data@long_name = fld@long_name
 res@xyLineColors    =   (/1,2,3,4,6,7,8/)
 plot  = gsn_csm_xy (wks,data,z,res) ; create plot

 do i=0,nfiles-1
   ff=f(i)
   fld = ff->W3({time1:time2},{0.:zmax})
   fld1 = ff->W2({time1:time2},{0.:zmax})
   fld1 = fld1^1.5
   fld = fld/fld1
   data(i,:) = dim_avg(fld(z|:,time|:))
 end do
 data@long_name = fld@long_name
 res@xyLineColors    =   (/1,2,3,4,6,7,8/)
 plot  = gsn_csm_xy (wks,data,z,res) ; create plot





end

