
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;=============================================
; data processing
;=============================================
 filename = "./out.nc"


 f     = addfile (filename , "r")             ; add file

; ===========================
 wks = gsn_open_wks("ps","ps")
 res                 = True                    ; plot mods desired

 res@cnFillOn        = True                    ; turn on color
 res@cnLinesOn       = False
 res@cnLineLabelsOn = False                    ; turn the line labels off

 res@gsnSpreadColors = True                    ; use full colormap
 res@cnFillMode            = "RasterFill"

 res@lbBoxLinesOn   = False
 res@lbOrientation   = "vertical"
 res@vpWidthF            = 1.0           ; change aspect ratio of plot
 res@vpHeightF           = 0.1


 res@gsnXAxisIrregular2Linear = True
 res@gsnYAxisIrregular2Linear = True
 res@gsnMaximize = True                    ; fill the page
; res@gsnPaperOrientation = "auto"

 fld = f->W(:,:)
 print("fld="+min(fld)+" "+max(fld))
 fld&z = fld&z * 0.001
; print(""+fld&z)
; print(""+fld&time)
 time = fld&time
 do i=1,dimsizes(time)-1
  if(time(i).lt.time(i-1)) then
   print("non-monotonic:"+time(i)+"  "+time(i-1)) 
  end if
 end do
 
; res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
 res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
 res@cnMinLevelValF          =  -5.             ; min level
 res@cnMaxLevelValF          =  5.             ; max level
 res@cnLevelSpacingF         =  0.1 ; interval
 res@lbLabelStride         = 5
 gsn_define_colormap(wks,"NCV_jaisnd")        ; choose color map

 plot = gsn_csm_contour(wks, fld(z|:,time|:), res) 


 fld = f->QV(:,:)
 print("fld="+min(fld)+" "+max(fld))
 fld&z = fld&z * 0.001
; res@cnLevelSelectionMode    = "AutomaticLevels"   ; manual levels
 res@cnLevelSelectionMode    = "ManualLevels"   ; manual levels
 res@cnMinLevelValF          =  0.             ; min level
 res@cnMaxLevelValF          =  11.             ; max level
 res@cnLevelSpacingF         =  0.1 ; interval
 res@lbLabelStride         = 5
 gsn_define_colormap(wks,"rainbow+white+gray")        ; choose color map


 plot = gsn_csm_contour(wks, fld(z|:,time|:), res)

