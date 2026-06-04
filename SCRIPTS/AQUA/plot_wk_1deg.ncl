load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

;==================================================================
;                     USER INTERFACE 
;                   MAIN DRIVER SCRIPT
;==================================================================
  case     = "AQUA 300K SAM1MOM 1x1"
  diro     = "./"     ; output dir: location of plots  [wkdir]
  diri     = "../../OUT_2D/"
;  fili     = "AQUA_2048x512x36_300K_P3_LWNT.2D_atm.nc"
  fili     = "AQUA_2048x512x36_300K_SAM1MOM_LWNT.2D_1deg.nc"
;  fili     = "AQUA_2048x512x36_QOBS_SAM1MOM_LWNT.2D_atm.nc"
;  fili     = "AQUA_2048x512x36_QOBS_P3_LWNT.2D_atm.nc"

  latBound = 10
  level    = -999     ; only used for 4D [set -999 otherwise]
;  VAR      = "Prec"    
  VAR      = "LWNT"    
  spd      = 1        ;     
;  nDayWin  = 60      ; Wheeler-Kiladis [WK] temporal window length (days)
  nDayWin  = 96      ; Wheeler-Kiladis [WK] temporal window length (days)
;  nDayWin  = 140      ; Wheeler-Kiladis [WK] temporal window length (days)
  nDaySkip = -20     ; time (days) between temporal windows [segments]
                     ; negative means there will be overlapping temporal segments
  opt      = True
;  opt@pltColorMap = "WhiteBlueGreenYellowRed"
  opt@cnLinesOn=False
;  opt@Fig_1 =(/-0.5,-0.25,0.,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.0,2.25,2.5/)
;  opt@Fig_2 =(/-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4/)
;  opt@Fig_3b =(/0.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11./)
;  opt@Fig_3a =(/0.,1.,1.5,2.,2.5,3.,3.5,4.,4.5/)

  wkSpaceTime_cam( diri        \ ; input directory
                 , fili        \ ; input file
                 , diro        \ ; output directory [plots]
                 , case        \ ; case name
                 , VAR         \ ; variable name
                 , latBound    \ ; NH lat bound [positive]
                 , spd         \ ; samples per day
                 , level       \ ; if 4D ... what level? set -999 if 3D
                 , nDayWin     \ ; temporal window length [days]
                 , nDaySkip    \ ; days between time windows 
                 , opt         \ ; options .... 
                 )
