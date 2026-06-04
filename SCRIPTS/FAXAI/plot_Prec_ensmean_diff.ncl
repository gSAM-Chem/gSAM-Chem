

file1 = "../../OUT_2D/ISABEL_PRECAC_MEAN.nc"
file2 = "../../OUT_2D/ISABEL_PRECAC_MEAN_DELTA.nc"


 f1 = addfile(file1,"r")
 f2 = addfile(file2,"r")

lon = f1->lon
lat = f1->lat

wks = gsn_open_wks("ps","ps")

res = True
res@cnFillOn              = True
res@cnLinesOn             = False
res@cnLineLabelsOn = False                    ; turn the line labels off
res@cnFillMode            = "RasterFill"

res@gsnAddCyclic          = False
res@mpLimitMode            = "LatLon"

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = True           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
;  resP@gsnPanelLabelBar                 = True                    ; add common colorbar
;  resP@gsnPanelFigureStringsFontHeightF = 0.015                   
;  resP@gsnMaximize = True  

  plot    = new((/4/),graphic)  


  resP@gsnPanelMainString  = "Precipitation Accumulation"

  res@mpOutlineOn           = True
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1" 

do mmm = 0,1

if(mmm.eq.0) then

lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)

else

lonmin = -83. + 360.
lonmax = -75. + 360.
latmin = 36.
latmax = 40.

end if
res@mpMinLatF              = latmin
res@mpMaxLatF              = latmax
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax

  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillPalette)
  res@cnLevelSelectionMode =  "ManualLevels"   
  colors = (/ (/255,255,255/), (/248,185,125/), \
              (/243,120,41/), (/240,47,29/), (/120,24,15/), \
              (/126,58,226/), (/24,68,124/)/) / 255.
  res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  res@cnLevels             = (/1.,3.,5.,7.,10.,15./)
  res@cnFillPalette        = colors

  fld1 = f1->Precac(0,:,:)
  fld1 = fld1/25.4
  print("Precac1: "+min(fld1)+"  "+max(fld1))
  res@gsnRightString = "Inches"
  res@gsnLeftString = "Past"
  res@cnLevelSpacingF      = 1.0
  res@cnMinLevelValF       = 0.
  res@cnMaxLevelValF       = 15.
  plot(0) = gsn_csm_contour_map(wks,fld1,res) 

  fld2 = f2->Precac(0,:,:)
  fld2 = fld2/25.4
  res@gsnLeftString = "Present"
  print("Precac2: "+min(fld2)+"  "+max(fld2))
  plot(1) = gsn_csm_contour_map(wks,fld2,res) 
  
  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillPalette)
  res@cnFillPalette        = "GHRSST_anomaly"
  fld = fld2
  fld = fld - fld1
  print("diffPrecac: "+min(fld)+"  "+max(fld))
  res@gsnLeftString = "Present - Past"
  res@cnLevelSelectionMode =  "ManualLevels"   
  res@cnLevelSpacingF      = 1.0
  res@cnMinLevelValF       = -6.
  res@cnMaxLevelValF       = 6.
  plot(2) = gsn_csm_contour_map(wks,fld,res) 

  fld1 = where(fld1.lt.0.1,0.1,fld1)
  fld = fld/fld1*100.
  fld = where(fld2.lt.1.,0.,fld)
  print("diffPrecac: "+min(fld)+"  "+max(fld))
  res@gsnLeftString = "Present - Past : Relative Change"
  res@gsnRightString = "%"
  res@cnLevelSpacingF      = 10.0
  res@cnMinLevelValF       = -100.
  res@cnMaxLevelValF       = 100.
  plot(3) = gsn_csm_contour_map(wks,fld,res) 

gsn_panel(wks,plot,(/2,2/),resP)                     ; now draw as one plot

end do ; mmm
