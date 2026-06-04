
begin
;
; Open file and get variables.
;
  a = addfile("../../OUT_2D/LI_3456x3456x83_100m_3s_SPNUDGE013610_2018-07-10-09-00-00_0000028800.2D_atm.nc","r")

  lat = a->lat
  lon = a->lon
;
; Terrain
;
  ter = a->ZSFC(0,:,:)
  land = a->LAND(0,:,:)
  ter@description = "Terrain Height"
  ter@units       = "m"
  ter = where(ter.eq.0.,-250.,ter)
;
; Get non-convective, convective
; Calculate total precipitation
;
  wks = gsn_open_wks("ps","ps")    ; Send graphics to PNG file

;
; Set up resource list that will be shared between the
; two wrf_contour calls.
;
  res                      = True

;  res@gsnMaximize           = True  
;  res@vpXF                 = 0.1   ; Set X,Y location of upper left
;  res@vpYF                 = 0.98   ; of plot
  res@vpWidthF             = 0.7   ; Set width
  res@vpHeightF             = 0.7   ; Set height
  res@mpShapeMode = "FreeAspect"

  res@tmXBTickSpacingF = 1              ; change maj lat tm spacing
  res@tmYLTickSpacingF = 1              ; change maj lon tm spacing
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01


  res@gsnAddCyclic         = False    ; Don't add lon cyclic point

  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnRightString       = ""
  res@gsnLeftString        = ""

  res@cnFillOn             = True
  res@cnLinesOn            = False
  res@cnFillMode            = "RasterFill"

  res@mpMinLonF            = min(lon)               ; select a subregion
  res@mpMaxLonF            = max(lon)
  res@mpMinLatF            = min(lat)
  res@mpMaxLatF            = max(lat)

  res@pmTickMarkDisplayMode = "Always"     ; Turn on map tickmarks.
  res@mpOutlineOn = True

;
; Generate plot of terrain for background.
;
; First set up resource list specific to terrain plot.
;
  opts_ter                        = res
  opts_ter@lbTitleString          = "Terrain (m)"
  opts_ter@lbTitlePosition        = "Bottom"  
  opts_ter@lbTitleFontHeightF     = 0.02
  opts_ter@lbLabelFontHeightF     = 0.02
  opts_ter@pmLabelBarWidthF       = 0.6
  opts_ter@pmLabelBarHeightF      = 0.1

;  opts_ter = wrf_map_resources(a,opts_ter)

;---If you want to see what wrf_map_resources does, then uncomment the next line.
;  print(opts_ter)

 opts_ter@mpOutlineBoundarySets       = "AllBoundaries"
 opts_ter@mpGeophysicalLineColor      = "Black"
 opts_ter@mpNationalLineColor         = "Black"
 opts_ter@mpUSStateLineColor          = "Black"
 opts_ter@mpGridLineColor             = "Black"
 opts_ter@mpLimbLineColor             = "Black"
 opts_ter@mpPerimLineColor            = "Black"

;  cmap1              = read_colormap_file("MPL_gist_earth")
  cmap1              = read_colormap_file("MPL_terrain")
;  cmap1              = read_colormap_file("cmocean_turbid")
;  cmap1(0,0:2)         = ((/ 0.,  0., 1./))
  opts_ter@cnLevelSelectionMode = "ManualLevels"
  opts_ter@cnMinLevelValF          =  0.             ; min level
  opts_ter@cnMaxLevelValF          = 1200.              ; max level
  opts_ter@cnLevelSpacingF         = 50.
  opts_ter@lbLabelStride    = 4


  opts_ter@mpOutlineBoundarySets = "AllBoundaries"
  opts_ter@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  opts_ter@mpDataSetName = "Earth..2"

  opts_ter@cnFillPalette  = cmap1

  map_terrain = gsn_csm_contour_map(wks,ter,opts_ter)

  ; ============================================================
  ; Read station lon/lat list
  ; File format: lon lat (one per line)
  ; ============================================================
  data = asciiread("/glade/work/marat/LI_STATIONS/stations_xy.txt", -1, "float")

  npts = dimsizes(data) / 2
  lons  = data(0:2*npts-2:2)+360.
  lats  = data(1:2*npts-1:2)

  print("Stations plotted: " + npts)

  ; Station index labels: 1..npts
  labels = new(npts, string)
  do i = 0, npts-1
    labels(i) = tostring(i+1)
  end do

  ; ============================================================
  ; Plot station markers
  ; ============================================================
  mkres = True
  mkres@gsMarkerIndex = 16      ; filled circle
  mkres@gsMarkerSizeF = 0.007
  mkres@gsMarkerColor = "red"

  draw(map_terrain)
  gsn_polymarker(wks, map_terrain, lons, lats, mkres)

; ============================================================
  ; Plot station index numbers (robust method)
  ; ============================================================
  txres = True
  txres@txFontHeightF = 0.01   ; make it visible
  txres@txFontColor   = "red"
  txres@txJust        = "CenterLeft"

  ; offset so numbers don't overlap dots
  dlon = -0.03
  dlat = -0.05

  do i = 0, npts-1
    gsn_text(wks, map_terrain, labels(i), lons(i)+dlon, lats(i)+dlat, txres)
  end do



  frame(wks)


info = (/                             \ ;Majority_Land_Cover_Type_1
        " 0 water"                     ,\ 
        " 1 evergreen needleleaf"      ,\ ;_forest
        " 2 evergreen broadleaf"       ,\ ;_forest     
        " 3 deciduous needleleaf"      ,\ ;_forest
        " 4 deciduous broadleaf"       ,\ ;_forest     
        " 5 mixed forests"             ,\
        " 6 closed shrubland"          ,\
        " 7 open shrublands"           ,\
        " 8 woody savannas"            ,\
        " 9 savannas"                  ,\
        " 10 grasslands"               ,\
        " 11 permanent wetlands"       ,\
        " 12 croplands"                ,\
        " 13 urban and built_up"       ,\
        " 14 crop nat veg_mosaic"      ,\ ; cropland_natural_vegetation_mosaic
        " 15 snow and ice"             ,\
        " 16 barren or sparse"         /) ; barren_or_sparsely_vegetated       

  colorscheme = (/"blue4", "darkgreen", "chartreuse", "olivedrab3","palegreen",\
		"darkseagreen","mediumpurple1","mistyrose","burlywood3","gold",\
		"tan1","royalblue","lightgoldenrod1","red","lightgoldenrod4","white",\
		"brown"/) 
  ninfo = dimsizes(info)        ; # of classifications

  res@cnLevelSelectionMode = "ExplicitLevels" ; set explict contour levels
  res@cnLevels         = tobyte( ispan(1,ninfo-1,1) )
  res@cnFillPalette    = colorscheme        ; distinct colors for categories
  res@gsnSpreadColors  = False              ; use each color sequentially


  res@lbLabelStrings   = ispan(0,ninfo,1)
  res@lbLabelPosition  = "Center"           ; label position
  res@lbLabelAlignment = "BoxCenters"       ; label orientation
  res@lbLabelStrings   = ispan(0,ninfo,1)  
  res@lbLabelStride    = 1
  res@lbLabelAutoStride = False
  res@lbLabelFontHeightF = 0.012

  res@pmLabelBarHeightF        =  0.075
  res@pmLabelBarWidthF         =  0.60      ; default is 0.6
  res@pmLabelBarOrthogonalPosF =  0.05 
  res@pmLabelBarParallelPosF = 0.5

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..2"


  plot = gsn_csm_contour_map(wks, land, res)

  rtxt = True
  rtxt@txJust        = "TopLeft"
  rtxt@txFontHeightF = 0.01

  nrow = 4       ; # rows
  ncol = 4       ; # columns

  n  = 0        ; counter
  xx = 0.1     ; iterate to get right  
  do nc=0,ncol-1
     yy = 0.060   ; 0.065 
    do nr=0,nrow-1
       n = n+1
       if (n.le.(ninfo-1)) then    ; handle less than ninfo
           gsn_text_ndc (wks,info(n),xx,yy,rtxt)
           yy = yy - 1.5*rtxt@txFontHeightF
       end if
    end do
     xx = xx + 0.200    ; iterate to get right 
  end do
  draw(plot)
  frame(wks)

end
