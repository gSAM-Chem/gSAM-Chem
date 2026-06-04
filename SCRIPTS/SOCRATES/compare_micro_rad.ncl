; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

nfiles = 6
f=new((/nfiles/),file)
f(0) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_M2005_576_2018-02-17.nc","r")
f(1) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_M2005_noradsnow_576_2018-02-17.nc","r")
f(2) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_THOM_576_2018-02-17.nc","r")
;f(3) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_THOM_noradsnow_576_2018-02-17.nc","r")
;f(4) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_P3_576_2018-02-17.nc","r")
;f(5) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_P3_radallice_576_2018-02-17.nc","r")
f(3) = addfile("../../../SAM7.4_FORK/OUT_2D/SOCRATES_576x288x74_10s_x16_SAM1MOM_576_17Feb.nc","r")
f(4) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_SAM1MOM_576_2018-02-17.nc","r")
f(5) = addfile("../../OUT_2D/SOCRATES_576x288x74_10s_x16_SAM1MOM_oldflx_576_2018-02-17.nc","r")
plot = new(nfiles,graphic)

names = (/"M2005","M2005 NoRadSnow","THOM","SAM1MOM SAM7.4","SAM1MOM","SAM1MOM OLDFLX"/)
;;names = (/"M2005","M2005 NoRadSnow","THOM","THOM NoRadSnow","P3","SAM1MOM"/)
;names = (/"M2005","M2005 NoRadSnow","THOM","THOM NoRadSnow","P3","P3 All Ice"/)

; produce area weights to compute zonal mean

fi = f(0)
www = fi->LWNT(0,:,:)
lat = fi->lat
lon = fi->lon
nlat = dimsizes(lat)
pi = acos(-1.)
a = sum(cos(lat*pi/180.))
do j=0,nlat-1
  www(j,:) = cos(lat(j)*pi/180.)
end do

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
res@cnLinesOn            = False
res@gsnDraw         = False               ; Do not draw plot
res@gsnFrame        = False               ; Do not advance frame
res@cnLevelSelectionMode = "ManualLevels"
;res@cnLevelSelectionMode = "AutomaticLevels"
res@cnFillPalette        = "MPL_rainbow"
;res@cnFillPalette        = "BlueWhiteOrangeRed"
res@cnRasterModeOn = True
;==================================================================
; LWNT
do k=0,nfiles-1
  fi = f(k)
  fld=fi->LWNT(0,:,:)
  res@cnMinLevelValF       =  120.                 ; set min contour level
  res@cnMaxLevelValF       =  320.                  ; set max contour level
  res@cnLevelSpacingF      =   10.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do 
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "LW NET TOA (W/m2)"     
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
;==================================================================
; SWNT
do k=0,nfiles-1
  fi = f(k)
  fld=fi->SWNT(0,:,:)
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  400.                  ; set max contour level
  res@cnLevelSpacingF      =   20.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "SW NET TOA (W/m2)" 
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
; LWCF
do k=0,nfiles-1
  fi = f(k)
  fld=fi->LWNTC(0,:,:) - fi->LWNT(0,:,:) 
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  100.                  ; set max contour level
  res@cnLevelSpacingF      =   5.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "LWCF (W/m2)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
;==================================================================
; SWCF
do k=0,nfiles-1
  fi = f(k)
  fld=fi->SWNT(0,:,:) - fi->SWNTC(0,:,:)
  res@cnMinLevelValF       =  -200.                 ; set min contour level
  res@cnMaxLevelValF       =  0.                  ; set max contour level
  res@cnLevelSpacingF      =   10.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "SWCF (W/m2)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
;==================================================================
; LWP 
do k=0,nfiles-1
  fi = f(k)
  fld=fi->CWP(0,:,:)+fi->RWP(0,:,:)
  fld = fld*1000.
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  200.                  ; set max contour level
  res@cnLevelSpacingF      =  10.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "LWP (ALL) g/kg" 
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
;==================================================================
; CWP
do k=0,nfiles-1
  fi = f(k)
  fld=fi->CWP(0,:,:)
  fld = fld*1000.
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  200.                  ; set max contour level
  res@cnLevelSpacingF      =  10.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "CWP  g/kg"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
; RWP
do k=0,nfiles-1
  fi = f(k)
  fld=fi->RWP(0,:,:)
  fld = fld*1000.
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  200.                  ; set max contour level
  res@cnLevelSpacingF      =  10.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "RWP  g/kg"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
; IWP
do k=0,nfiles-1
  fi = f(k)
  fld=fi->IWP(0,:,:)+fi->SWP(0,:,:)+fi->GWP(0,:,:)
  fld = fld*1000.
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  300.                  ; set max contour level
  res@cnLevelSpacingF      =   30.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "IWP (ALL)  (g/kg)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
;==================================================================
; IWP
do k=0,nfiles-1
  fi = f(k)
  fld=fi->IWP(0,:,:)
  fld = fld*1000.
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  300.                  ; set max contour level
  res@cnLevelSpacingF      =   30.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "IWP (g/kg)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================
; GWP
do k=0,nfiles-1
  fi = f(k)
  fld=fi->GWP(0,:,:)
  fld = fld*1000.
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  300.                  ; set max contour level
  res@cnLevelSpacingF      =   30.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "GWP (g/kg)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down

;==================================================================
; SWP
do k=0,nfiles-1
  fi = f(k)
  fld=fi->SWP(0,:,:)
  fld = fld*1000.
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  300.                  ; set max contour level
  res@cnLevelSpacingF      =   30.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "SWP (g/kg)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down

;==================================================================
; PW
do k=0,nfiles-1
  fi = f(k)
  fld=fi->PW(0,:,:)
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  60.                  ; set max contour level
  res@cnLevelSpacingF      =  4.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "Precipitable Water (mm)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================


;==================================================================
; PREC
do k=0,nfiles-1
  fi = f(k)
  fld=fi->Preca(0,:,:)
  res@cnMinLevelValF       =  0.                 ; set min contour level
  res@cnMaxLevelValF       =  20.                  ; set max contour level
  res@cnLevelSpacingF      =  1.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "Precipitation (mm/day)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================

;==================================================================
; LHF
do k=0,nfiles-1
  fi = f(k)
  fld=fi->LHF(0,:,:)
  res@cnMinLevelValF       =  -100.                 ; set min contour level
  res@cnMaxLevelValF       =  400.                  ; set max contour level
  res@cnLevelSpacingF      =  20.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "Latent Heat Flux (W/m2)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================

;==================================================================
; LHF
do k=0,nfiles-1
  fi = f(k)
  fld=fi->SHF(0,:,:)
  res@cnMinLevelValF       =  -100.                 ; set min contour level
  res@cnMaxLevelValF       =  100.                  ; set max contour level
  res@cnLevelSpacingF      =  10.                  ; set contour spacing
  fld@units=dim_avg_Wrap(dim_avg_Wrap(fld*www))
  fld@long_name = names(k)
  plot(k) = gsn_csm_contour_map(wks,fld,res)        ; create plot
end do
resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = "Sensible Heat Flux (W/m2)"
gsn_panel(wks,plot,(/nfiles/2,2/),resp)           ; draw: 2-across, 2-down
;==================================================================

