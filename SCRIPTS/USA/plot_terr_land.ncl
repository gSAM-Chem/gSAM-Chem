
begin
;
; Open file and get variables.
;
  a = addfile("../../OUT_2D/USA_3072x1536x72_2.5km_10s_APR13_2019-04-14-00-00-00_0000008640.2D_atm.nc","r")

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

  res@vpXF                 = 0.1   ; Set X,Y location of upper left
  res@vpYF                 = 0.98   ; of plot
  res@vpWidthF             = 0.80   ; Set width
  res@vpHeightF             = 0.40   ; Set height

  res@tmXBTickSpacingF = 10              ; change maj lat tm spacing
  res@tmYLTickSpacingF = 10              ; change maj lon tm spacing
  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015


  res@gsnAddCyclic         = False    ; Don't add lon cyclic point

  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnRightString       = ""
  res@gsnLeftString        = ""

  res@cnFillOn             = True
  res@cnLinesOn            = False
  res@cnFillMode            = "RasterFill"

  res@trMinLonF            = 285.8               ; select a subregion
  res@trMaxLonF            = 287.8
  res@trMinLatF            = 40.
  res@trMaxLatF            = 42.

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

;---wrf_map_resources sets these to gray.
;  opts_ter@mpOutlineBoundarySets       = "AllBoundaries"
; opts_ter@mpGeophysicalLineColor      = "Black"
; opts_ter@mpNationalLineColor         = "Black"
;  opts_ter@mpUSStateLineColor          = "Black"
;  opts_ter@mpGridLineColor             = "Black"
;  opts_ter@mpLimbLineColor             = "Black"
;  opts_ter@mpPerimLineColor            = "Black"

;  cmap1              = read_colormap_file("MPL_gist_earth")
  cmap1              = read_colormap_file("MPL_terrain")
;  cmap1              = read_colormap_file("cmocean_turbid")
;  cmap1(0,0:2)         = ((/ 0.,  0., 1./))
  opts_ter@cnLevelSelectionMode = "ManualLevels"
  opts_ter@cnMinLevelValF          =  0.             ; min level
  opts_ter@cnMaxLevelValF          = 4000.              ; max level
  opts_ter@cnLevelSpacingF         = 200.
  opts_ter@lbLabelStride    = 4

  opts_ter@mpOutlineBoundarySets = "AllBoundaries"
  opts_ter@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  opts_ter@mpDataSetName = "Earth..1"


  opts_ter@cnFillPalette  = cmap1

  map_terrain = gsn_csm_contour(wks,ter,opts_ter)

  draw(map_terrain)
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
  res@lbLabelFontHeightF = 0.01

  res@pmLabelBarHeightF        =  0.075
  res@pmLabelBarWidthF         =  0.60      ; default is 0.6
  res@pmLabelBarOrthogonalPosF =  0.03; move up smidge

  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataBaseVersion     = "MediumRes"   ; necessary for mpDataSetName to be effective
  res@mpDataSetName = "Earth..1"


  plot = gsn_csm_contour(wks, land, res)

  rtxt = True
  rtxt@txJust        = "TopLeft"
  rtxt@txFontHeightF = 0.015

  nrow = 4       ; # rows
  ncol = 4       ; # columns

  n  = 0        ; counter
  xx = 0.01     ; iterate to get right  
  do nc=0,ncol-1
     yy = 0.120   ; 0.065 
    do nr=0,nrow-1
       n = n+1
       if (n.le.(ninfo-1)) then    ; handle less than ninfo
           gsn_text_ndc (wks,info(n),xx,yy,rtxt)
           yy = yy - 2*rtxt@txFontHeightF
       end if
    end do
     xx = xx + 0.2500    ; iterate to get right 
  end do
  draw(plot)
  frame(wks)

end
