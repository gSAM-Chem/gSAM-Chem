
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;=============================================
; data processing
;=============================================
 filename = "../../OUT_3D/CLOUD_72x72x74_P3_3D.nc"
; filename = "../../OUT_3D/CLOUD_72x72x74_M2005_3D.nc"

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

 f = addfile(filename,"r")

 nx = dimsizes(f->lon)
 ny = dimsizes(f->lat)

 fld = f->QV(:,{0:height_max},ny/2+1,nx/2+1)
 fld&z = fld&z * 0.001
 fld&time = (fld&time - fld&time(0)) * 86400.
; res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
 res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
 res@cnMinLevelValF          =  0.             ; min level
 res@cnMaxLevelValF          =  16.             ; max level
 res@cnLevelSpacingF         =  2. ; interval
 plot = gsn_csm_contour(wks, fld(z|:,time|:), res) 

 fld = f->QC(:,{0:height_max},ny/2+1,nx/2+1)
 fld&z = fld&z * 0.001
 fld&time = (fld&time - fld&time(0)) * 86400.
; res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
 res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
 res@cnMinLevelValF          =  0.             ; min level
 res@cnMaxLevelValF          =  3.             ; max level
 res@cnLevelSpacingF         =  0.1 ; interval
 plot = gsn_csm_contour(wks, fld(z|:,time|:), res)  

 fld = f->QR(:,{0:height_max},ny/2+1,nx/2+1)
 fld&z = fld&z * 0.001
 fld&time = (fld&time - fld&time(0)) * 86400.
; res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
 res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
 res@cnMinLevelValF          =  0.             ; min level
 res@cnMaxLevelValF          =  3.             ; max level
 res@cnLevelSpacingF         =  0.1 ; interval
 plot = gsn_csm_contour(wks, fld(z|:,time|:), res)
