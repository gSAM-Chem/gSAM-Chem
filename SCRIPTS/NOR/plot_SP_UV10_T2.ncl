;----------------------------------------------------------------------
; wrf_title_6.ncl
;----------------------------------------------------------------------
; Concepts illustrated:
;   - Annotating titles in plots created by wrf_xxxx functions
;----------------------------------------------------------------------
; This script shows how to overlay line contours and wind barbs on
; top of a filled contour plot over a map.
;----------------------------------------------------------------------
; These files are loaded by default in NCL V6.4.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;----------------------------------------------------------------------

;----------------------------------------------------------------------
; This procedure annotates various WRF plot titles with boxes and
; explanatory text.
;----------------------------------------------------------------------
procedure annotate_wrf_plot(wks)
local lnres, txres, xbox, ybox
begin
;  drawNDCGrid(wks)   ; Draw a labeled unit square grid to help determine
                      ; NDC values for boxes and text.
;---Line options
  lnres                       = True
  lnres@gsLineThicknessF      = 3.0
  lnres@gsLineColor           = "brown"

;---Text options
  txres                       = True
  txres@txFontHeightF         = 0.010
  txres@txFont                = "Helvetica-bold"
  txres@txPerimOn             = True
  txres@txBackgroundFillColor = "darkseagreen1" ; "yellow"

;---Annotate MainTitle area
  txres@txJust = "CenterLeft" 
  xbox = (/0.19,0.37,0.37,0.19,0.19/)
  ybox = (/0.99,0.99,0.96,0.96,0.99/)
  gsn_polyline_ndc(wks,xbox,ybox,lnres)
  gsn_text_ndc(wks,"MainTitle",0.38,0.975,txres)

;---Annotate top right field titles area
  xbox = (/0.20,0.38,0.38,0.20,0.20/)
  ybox = (/0.93,0.93,0.87,0.87,0.93/)
  gsn_polyline_ndc(wks,xbox,ybox,lnres)
  gsn_text_ndc(wks,"FieldTitle (UnitsLabel)",0.39,0.92,txres)
  gsn_text_ndc(wks,"FieldTitle (UnitsLabel)",0.39,0.90,txres)
  gsn_text_ndc(wks,"FieldTitle (UnitsLabel)",0.39,0.88,txres)

;---Annotate Labelbar title area
  xbox = (/0.46,0.60,0.60,0.46,0.46/)
  ybox = (/0.24,0.24,0.215,0.215,0.24/)
  gsn_polyline_ndc(wks,xbox,ybox,lnres)
  gsn_text_ndc(wks,"FieldTitle (UnitsLabel)",0.61,0.228,txres)

;---Annotate info label area
  xbox = (/0.62,0.87,0.87,0.62,0.62/)
  ybox = (/0.29,0.29,0.26,0.26,0.29/)
  gsn_polyline_ndc(wks,xbox,ybox,lnres)
  txres@txJust = "CenterRight" 
  gsn_text_ndc(wks,"FieldTitle (UnitsLabel)",0.61,0.275,txres)

;---Annotate init title area
  xbox = (/0.69,0.87,0.87,0.69,0.69/)
  ybox = (/0.99,0.99,0.97,0.97,0.99/)
  gsn_polyline_ndc(wks,xbox,ybox,lnres)
  gsn_text_ndc(wks,"Init time title not modifiable~C~Turn off with InitTime=False",0.865,0.93,txres)

;---Annotate valid title area
  xbox = (/0.68,0.87,0.87,0.68,0.68/)
  ybox = (/0.97,0.97,0.95,0.95,0.97/)
  gsn_polyline_ndc(wks,xbox,ybox,lnres)
  txres@txJust = "TopRight" 
  gsn_text_ndc(wks,"ValidTime=False by default~C~Change with TimeLabel",0.865,0.90,txres)

;---Annotate footer area
  xbox = (/0.20,0.69,0.69,0.20,0.20/)
  ybox = (/0.042,0.042,0.01,0.01,0.042/)
  gsn_polyline_ndc(wks,xbox,ybox,lnres)
  txres@txJust = "BottomLeft"
  gsn_text_ndc(wks,"Bottom title is not modifiable~C~" +\
               "Turn off with Footer=False or NoHeaderFooter=False",0.21,0.053,txres)
end

;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
begin
;---Open a WRF output file and calculate several diagnostics
  f = addfile("../../OUT_2D/NOR_1536x1536x72_2.5km_10s_JAN3_SNOW_2018-01-04-09-00-00_0000011880.2D_atm.nc","r")
  slp   = f->PSFC(0,:,:)
  t2    = f->T2m(0,:,:)
  u10   = f->U10m(0,:,:)
  v10   = f->V10m(0,:,:)
  times = f->time  ; get all times in the file

;---Indicate where you want to set the graphical output. "svg" was used for the web.
  wks = gsn_open_wks("png","wrf_title")

;---Set resources common to all plots
  res           = True
  res@MainTitle = "REAL-TIME WRF"
  res@TimeLabel = times(0)   ; Set Valid time to use on plots

;--- Line contours
  os                  = res
  os@cnLineColor      = "NavyBlue"
  os@cnLineThicknessF = 2.0
  os@FieldTitle       = "SLP"
  plot_slp = wrf_contour(f,wks,slp,os)

;--- Filled contours
  ot             = res
  ot@cnFillOn    = True
; ot@FieldTitle = "T2 contours"
  plot_tc = wrf_contour(f,wks,t2,ot)

;--- Vectors
  ov             = res
  ov@NumVectors  = 47
; ov@FieldTitle = "U/V vectors"
  plot_vec = wrf_vector(f,wks,u10,v10,ov)

;--- Overlay vectors, line contours, and filled contours on a map
  pltres = True
  pltres@PanelPlot = True
  plot  = wrf_map_overlays(f,wks,(/plot_tc,plot_slp,plot_vec/),pltres,False)
  draw(plot)
  annotate_wrf_plot(wks)
  frame(wks)
end

