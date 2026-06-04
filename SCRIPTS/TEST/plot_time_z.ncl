
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;=============================================
; data processing
;=============================================
 filename = "../../OUT_STAT/TEST_32x32x83_TKE_NORAD.nc"
; filename = "../../OUT_STAT/TEST_32x32x83.nc"

height_max=15000.

; gsn_define_colormap(wks,"rainbow+white+gray")        ; choose color map

 f     = addfile (filename , "r")             ; add file

; ===========================
 wks = gsn_open_wks("ps","ps")
 res                 = True                    ; plot mods desired

 res@cnFillOn        = True                    ; turn on color
 res@cnLinesOn       = False
 res@gsnSpreadColors = True                    ; use full colormap

 res@lbLabelStride         = 2
 res@lbBoxLinesOn   = False
 res@lbOrientation   = "vertical"
; res@vpWidthF            = 1.0           ; change aspect ratio of plot
; res@vpHeightF           = 0.5


 res@gsnXAxisIrregular2Linear = True
 res@gsnYAxisIrregular2Linear = True
; res@gsnMaximize = True                    ; fill the page
; res@gsnPaperOrientation = "auto"

 res@trYMinF = 0.
 res@trYMaxF = height_max*0.001

 fld = f->QV(:,:)
 print("fld="+min(fld)+" "+max(fld))
 fld&z = fld&z * 0.001
 fld&time = (fld&time - fld&time(0)) * 86400.
 res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
; res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
; res@cnMinLevelValF          =  0.             ; min level
; res@cnMaxLevelValF          =  12.             ; max level
; res@cnLevelSpacingF         =  1. ; interval
 plot = gsn_csm_contour(wks, fld(z|:,time|:), res) 

 fld = f->TKE(:,:)
 print("fld="+min(fld)+" "+max(fld))
 fld&z = fld&z * 0.001
 fld&time = (fld&time - fld&time(0)) * 86400.
 res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
; res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
; res@cnMinLevelValF          =  0.             ; min level
; res@cnMaxLevelValF          =  3.             ; max level
; res@cnLevelSpacingF         =  0.1 ; interval
 plot = gsn_csm_contour(wks, fld(z|:,time|:), res)  

 fld = f->W2(:,:)
 print("fld="+min(fld)+" "+max(fld))
 fld&z = fld&z * 0.001
 fld&time = (fld&time - fld&time(0)) * 86400.
 res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
; res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
; res@cnMinLevelValF          =  0.             ; min level
; res@cnMaxLevelValF          =  1.             ; max level
; res@cnLevelSpacingF         =  0.1 ; interval
 plot = gsn_csm_contour(wks, fld(z|:,time|:), res)

 fld = f->TABS(:,:)
 print("fld="+min(fld)+" "+max(fld))
 fld&z = fld&z * 0.001
 fld&time = (fld&time - fld&time(0)) * 86400.
 res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
; res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
; res@cnMinLevelValF          =  288.             ; min level
; res@cnMaxLevelValF          =  294.             ; max level
; res@cnLevelSpacingF         =  0.5 ; interval
 plot = gsn_csm_contour(wks, fld(z|:,time|:), res)
