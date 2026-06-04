files = systemfunc("ls -1 ../../OUT_2DL/AGNES_1024x1024x72_3km_10s_ICEF6_NONUDGE_19*.2D_lnd.nc")
;files = systemfunc("ls -1 ../../OUT_2DL/AGNES_1024x1024x72_3km_10s_ICEF6_NUDGEUVWT_19*.2D_lnd.nc")
nfiles = dimsizes(files)

wks = gsn_open_wks("png","image")
;wks = gsn_open_wks("ps","ps")

f=addfiles(files,"r")

lat = f[0]->lat
lon = f[0]->lon
print("lon: "+min(lon)+" "+max(lon))
print("lat: "+min(lat)+" "+max(lat))

res = True
res@cnFillOn              = True
res@cnLinesOn             = False
res@cnLineLabelsOn = False                    ; turn the line labels off
res@cnFillMode            = "RasterFill"

lonmin = min(lon)
lonmax = max(lon)
latmin = min(lat)
latmax = max(lat)

res@gsnAddCyclic          = False
res@mpLimitMode            = "LatLon"
res@mpMinLatF              = min(lat)  
res@mpMaxLatF              = max(lat) 
res@mpMinLonF              = min(lon)
res@mpMaxLonF              = max(lon) 

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = False           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
  resP@gsnPanelLabelBar                 = True                    ; add common colorbar
  resP@gsnPanelFigureStringsFontHeightF = 0.015                   

  plot    = new(1,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   

  resP@gsnPanelMainString  = "On-grownd Water Depth"
  res@cnFillPalette       = "MPL_jet"

  res@mpOutlineOn           = True
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1" 


do i=0,nfiles-1
  fobs = f[i]->MWS(0,:,:)
  print("Precac: "+min(fobs)+"  "+max(fobs))
  res@gsnLeftString = " "
  res@cnMinLevelValF       = 0.
  res@cnMaxLevelValF       = 30.
  res@cnLevelSpacingF      = 1. 
  plot(0) = gsn_csm_contour_map(wks,fobs,res) 

  gsn_panel(wks,plot,(/1,1/),resP)                     ; now draw as one plot
end do
