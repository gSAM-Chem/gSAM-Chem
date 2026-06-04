filename = "/glade/work/marat/HRRR/DATA/hrrr2D_20180709_18z_0.06x0.06_deg_sfc.nc"

;wks = gsn_open_wks("png","image")
wks = gsn_open_wks("ps","ps")

f=addfile(filename,"r")

lat = f->lat
lon = f->lon
print("lon: "+min(lon)+" "+max(lon))
print("lat: "+min(lat)+" "+max(lat))

res = True
res@cnFillOn              = True
res@cnLinesOn             = False
res@cnLineLabelsOn = False                    ; turn the line labels off
res@cnFillMode            = "RasterFill"

lonmin = 283.9
lonmax = 288.06
latmin = 39.44
latmax = 42.55

res@gsnAddCyclic          = False
res@mpLimitMode            = "LatLon"
res@mpMinLatF              = latmin  
res@mpMaxLatF              = latmax 
res@mpMinLonF              = lonmin
res@mpMaxLonF              = lonmax 

res@gsnDraw             = False           ; don't draw
res@gsnFrame            = False           ; don't advance frame
res@lbLabelBarOn        = True           ; turn off individual cb's

;************************************************
  resP                                  = True                    ; modify the panel plot
  resP@gsnPanelLabelBar                 = False                    ; add common colorbar
;  resP@gsnPanelFigureStrings            = (/"gSAM","ERA5"/)           ; add strings to panel
;  resP@gsnPanelFigureStringsFontHeightF = 0.015                   

  plot    = new(1,graphic)  

  res@cnLevelSelectionMode =  "ManualLevels"   

  resP@gsnPanelMainString  = "2-M Temperature and RH"
  res@cnFillPalette       = "MPL_jet"

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"


  fld = f->T2m(:,:)
  fld = (fld-32.)/(9./5.)
  fld1 = dim_avg_Wrap(fld)
  lat = f->lat
  w = cos(lat * 0.0174533)
  w2D = conform(fld1, w, 0)
  fld_mean = dim_avg_wgt_n(fld1,w2D,0,0)
  print("mean="+fld_mean)
  res@gsnLeftString = "T2m "+(fld_mean-273.15)
  fld2 = fld
  res@cnMinLevelValF       = 12.
  res@cnMaxLevelValF       = 35.
  res@cnLevelSpacingF      = 1. 
  plot(0) = gsn_csm_contour_map(wks,fld2,res) 

  gsn_panel(wks,plot,(/1,1/),resP)                     ; now draw as one plot
