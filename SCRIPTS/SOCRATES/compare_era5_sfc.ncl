; ***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

f=addfile("../../OUT_2D/SOCRATES_9216x4608x74_THOM_4608_20180219020000_0000026640.2D_1.nc","r")
fobs=addfile("../../GLOBAL_DATA/NC/era5_sfc_576x288_dyvar_2018021600.nc","r")

nobs = 0

wks = gsn_open_wks("ps","ps")                ; send graphics to PNG file
res                      = True                  ; plot mods desired
res@cnFillOn             = True                  ; turn on color fill
res@cnLinesOn            = False
res@gsnDraw         = False               ; Do not draw plot
res@gsnFrame        = False               ; Do not advance frame
res@cnLevelSelectionMode = "ManualLevels"
plot = new(3,graphic)
res@cnFillPalette        = "BlueWhiteOrangeRed"
res@cnRasterModeOn = True
;==================================================================
; SKIN TEMPERATURE

fld=f->SKT(0,:,:)
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->SKT(:,:)
fld_obs = fld_obs-273.15
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -20.                 ; set min contour level
res@cnMaxLevelValF       =  20.                  ; set max contour level
res@cnLevelSpacingF      =   2.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Skin Temparature, C"     
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; T2M

fld=f->T2m(0,:,:)
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->T2m(:,:)
fld_obs = fld_obs-273.15
res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -20.                 ; set min contour level
res@cnMaxLevelValF       =  20.                  ; set max contour level
res@cnLevelSpacingF      =   2.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " 2m Temparature, C"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down
                                                                   
;==================================================================
; U10M 

fld=f->U10m(0,:,:)
res@cnMinLevelValF       =  -20.                 ; set min contour level
res@cnMaxLevelValF       =  20.                  ; set max contour level
res@cnLevelSpacingF      =   2.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->U10m(:,:)
res@cnMinLevelValF       =  -20.                 ; set min contour level
res@cnMaxLevelValF       =  20.                  ; set max contour level
res@cnLevelSpacingF      =   2.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " 10m Zonal Wind, m/s"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; V10M

fld=f->V10m(0,:,:)
res@cnMinLevelValF       =  -20.                 ; set min contour level
res@cnMaxLevelValF       =  20.                  ; set max contour level
res@cnLevelSpacingF      =   2.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->V10n(:,:)
res@cnMinLevelValF       =  -20.                 ; set min contour level
res@cnMaxLevelValF       =  20.                  ; set max contour level
res@cnLevelSpacingF      =   2.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -10.                 ; set min contour level
res@cnMaxLevelValF       =  10.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " 10m Meridional Wind, m/s"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; SHF

fld=f->SHF(0,:,:)
res@cnMinLevelValF       =  -200.                 ; set min contour level
res@cnMaxLevelValF       =  200.                  ; set max contour level
res@cnLevelSpacingF      =   20.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->SHF(:,:)
fld_obs = -fld_obs/3600.
res@cnMinLevelValF       =  -200.                 ; set min contour level
res@cnMaxLevelValF       =  200.                  ; set max contour level
res@cnLevelSpacingF      =   20.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -100.                 ; set min contour level
res@cnMaxLevelValF       =  100.                  ; set max contour level
res@cnLevelSpacingF      =   10.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " SHF, W/m2"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; LHF

fld=f->LHF(0,:,:)
res@cnMinLevelValF       =  -200.                 ; set min contour level
res@cnMaxLevelValF       =  200.                  ; set max contour level
res@cnLevelSpacingF      =   20.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->LHF(:,:)
fld_obs = -fld_obs/3600.
res@cnMinLevelValF       =  -200.                 ; set min contour level
res@cnMaxLevelValF       =  200.                  ; set max contour level
res@cnLevelSpacingF      =   20.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -100.                 ; set min contour level
res@cnMaxLevelValF       =  100.                  ; set max contour level
res@cnLevelSpacingF      =   10.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " LHF, W/m2"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; SOLIN

fld=f->SOLINA(0,:,:)
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =   50.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->tisr(:,:)
fld_obs = fld_obs/3600.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =   50.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -100.                 ; set min contour level
res@cnMaxLevelValF       =  100.                  ; set max contour level
res@cnLevelSpacingF      =   10.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " SOLIN, W/m2"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; SWNT

fld=f->SWNTA(0,:,:)
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =   50.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->tsr(:,:)
fld_obs = fld_obs/3600.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  1000.                  ; set max contour level
res@cnLevelSpacingF      =   50.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -100.                 ; set min contour level
res@cnMaxLevelValF       =  100.                  ; set max contour level
res@cnLevelSpacingF      =   10.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " SWNT, W/m2"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; LWNT

fld=f->LWNTA(0,:,:)
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  350.                  ; set max contour level
res@cnLevelSpacingF      =   10.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->ttr(:,:)
fld_obs = -fld_obs/3600.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  350.                  ; set max contour level
res@cnLevelSpacingF      =   10.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -50.                 ; set min contour level
res@cnMaxLevelValF       =  50.                  ; set max contour level
res@cnLevelSpacingF      =   5.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " LWNT, W/m2"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

;==================================================================
; PREC

res@cnFillPalette        = "precip2_17lev"

fld=f->Preca(0,:,:)
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  17.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
fld@units="SAM"
plot(0) = gsn_csm_contour_map(wks,fld,res)        ; create plot

fld_obs=fobs->tp(:,:)
fld_obs = fld_obs*1000./3600.*86400.
res@cnMinLevelValF       =  0.                 ; set min contour level
res@cnMaxLevelValF       =  17.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
fld_obs@units="ERA5"
plot(1) = gsn_csm_contour_map(wks,fld_obs,res)        ; create plot

res@cnMinLevelValF       =  -5.                 ; set min contour level
res@cnMaxLevelValF       =  5.                  ; set max contour level
res@cnLevelSpacingF      =   1.                  ; set contour spacing
fld@units="SAM-ERA5"
fld = fld-fld_obs
plot(2) = gsn_csm_contour_map(wks,fld,res)        ; create plot

resp             = True                    ; panel mods desired
resp@gsnMaximize = True                    ; fill the page
resp@txString   = " Precipitation, mm/day"
gsn_panel(wks,plot,(/3,1/),resp)           ; draw: 2-across, 2-down

